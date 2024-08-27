#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>

#include <mpi.h>

typedef double real_t;
typedef int64_t int_t;
#define INT_FMT "lld"
#define REAL_FMT "lf"


/* Simulation parameters. Default values can be overridden by options(). */
int_t
    height    = 400,        /* Vertical coordinates (index i)   */
    width     = 600,        /* Horizontal coordinates (index j) */
    snap_freq = 100,        /* Number of iterations between saved states */
    max_iter  = 40001;      /* Total number of iterations */
bool
    quiet = false;  /* Suppresses stdio and file I/O when true */


/* Map of the simulated domain, this defines the problem geometry.
 * initialize_domain() sets up a set of SOLID/FLUID points, and
 * then calculates which of the SOLID points are actually WALL,
 * on account of having at least one FLUID neighbor.
 * SOLID points are omitted by the compute kernels, as they have no
 * direct impact on fluid flow.
 */
typedef enum { FLUID, WALL, SOLID } node_type_t;
node_type_t *map = NULL;
#define MAP(i,j) map[(i)*width+(j)]


/* Lattice geometry and physical data */
real_t
        /* Density distributions for present and next iteration */
    *density[2] = { NULL, NULL },

        /* 2D velocity vectors derived from the present density distribution */
    *velocity   = NULL,

        /* Magnitude of velocity vectors, for visualization purposes */
    *abs_velocity = NULL,

        /* Sine and cosine values of the 6 lattice directions */
    c[6][2],

        /* Weight coefficient of relaxation operator */
    lambda = -1.0,

        /* Artificial external force of inflow */
    force[2] = { 0.0, 1e-2 },

        /* Radius of the bluff obstruction.
         * This value should be positive, its default value is computed by
         * init_domain() unless it is overridden by a command line argument.
         */
    radius = 0.0;

/* y and x coordinate shifts to neighboring lattice nodes */
int_t offsets[2][6][2] = {
    { {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0} }, /* Even rows */
    { {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1} }  /* Odd rows */
};

/* MPI Rank and size */
int rank, size;

float *buf = NULL;
#define Buf(i,j) buf[(i)*(width+1)+(j)]


#define D_now(i,j,d) density[0][(d)*width*height/size+(i)*width+(j)]
#define D_nxt(i,j,d) density[1][(d)*width*height/size+(i)*width+(j)]
#define V(i,j,x) velocity[2*((i)*width+(j))+(x)]
#define V_abs(i,j) abs_velocity[(i)*width+(j)]


/* Functions used in the time integration loop */
    /* Configure problem geometry: periodic horizontal channel with a
     * bluff (circular) obstruction in the middle of the flow.
     */
void initialize_domain ( void );
    /* Relax the densities at all lattice sites toward equilibrium. */
void collide ( void );
    /* Propagate densities from each lattice site to its neighbors */
void stream ( void );
    /* Save simulation state in gnuplot binary matrix format.
     * File names are numbered by iter / snap_freq, and stored in 'data/'
     */
void store ( int_t iter );


/* Command line option parsing */
void options ( int argc, char **argv );


/* Main program stages:
 *  - argument parsing
 *  - memory allocation
 *  - initialization
 *  - time integration loop
 *  - deallocation
 */
int
main ( int argc, char **argv )
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        options ( argc, argv );
        for (int_t i = 1; i < size; ++i) {
            MPI_Send(&height, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&snap_freq, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&max_iter, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&snap_freq, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_iter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    /* Make all required dynamic memory allocations */
    density[0] = malloc ( 6 * width * height / size * sizeof(real_t) );
    density[1] = malloc ( 6 * width * height / size * sizeof(real_t) );
    velocity = malloc ( 2 * width * height / size * sizeof(real_t) );
    abs_velocity = malloc ( width * height / size * sizeof(real_t) );
    map = malloc ( width * height / size * sizeof(node_type_t) );
    buf = malloc( (width+1) * (height/size+1) * sizeof(float));

    /* Initialize density distribution to equilibrium, i.e. assume
     * a density of 1 per lattice point, and distribute it evenly
     * in 6 directions.
     */
    for ( int_t i=0; i<height/size; i++ )
        for ( int_t j=0; j<width; j++ )
            for ( int_t d=0; d<6; d++ )
                D_now(i,j,d) = D_nxt(i,j,d) = 1.0 / 6.0;

    /* Initialize direction vectors for the 6 directions at each point
     * to match a hexagonal lattice.
     */
    for( int_t d=0; d<6; d++ )
    {
        c[d][0] = sin ( M_PI * d / 3.0 );
        c[d][1] = cos ( M_PI * d / 3.0 );
    }

    /* Set up the simulated geometry */
    initialize_domain();

    /* Time integration loop */
    for ( int_t iter=0; iter<max_iter; iter++ )
    {
        collide();
        stream();

        /* Dump velocity field for visualization every snap_freq steps */
        if ( quiet==false && (iter % snap_freq) == 0 )
        {
            store ( iter / snap_freq );
            if (quiet == false && rank == 0)
                printf ( "Iteration %" INT_FMT " / %" INT_FMT " (%.f %%) \n",
                    iter, max_iter, 100.0*iter/(float)max_iter
                );
        }
    }

    /* Store only last state if I/O is disabled */
    if ( quiet == true )
        store ( max_iter / snap_freq );

    MPI_Finalize();

    /* Release the memory and halt */
    free ( density[0] );
    free ( density[1] );
    free ( velocity );
    free ( abs_velocity );
    free ( map );
    free ( buf );
    exit ( EXIT_SUCCESS );
}


/* Collision/relaxation operator:
 * redistribute the 6 densities at each lattice point so that they
 * approach equilibrium.
 */
void
collide ( void )
{
    for ( int_t i=0; i<height/size; i++ )
    {
        for ( int_t j=0; j<width; j++ )
        {
            /* Disregard solid points */
            if ( MAP(i,j) == FLUID || MAP(i,j) == WALL )
            {
                real_t rho = 0.0, uc = 0.0;
                V(i,j,0) = V(i,j,1) = 0.0;

                /* Calculate the current velocity field in the fluid */
                if ( MAP(i,j) == FLUID )
                {
                    for ( int_t d=0; d<6; d++ )
                    {
                        rho += D_now(i,j,d);
                        V(i,j,0) += c[d][0] * D_now(i,j,d);
                        V(i,j,1) += c[d][1] * D_now(i,j,d);
                    }
                    V(i,j,0) /= rho;
                    V(i,j,1) /= rho;
                }

                /* Calculate the D2Q6 lattice collision operator */
                for ( int_t d=0; d<6; d++ )
                {
                    real_t qi_uaub, N_eq, delta_N;
                    qi_uaub =
                        (c[d][1] * c[d][1] - 0.5) * V(i,j,1) * V(i,j,1)
                      + (c[d][1] * c[d][0])       * V(i,j,1) * V(i,j,0)
                      + (c[d][0] * c[d][1])       * V(i,j,0) * V(i,j,1)
                      + (c[d][0] * c[d][0] - 0.5) * V(i,j,0) * V(i,j,0)
                    ;
                    uc = c[d][0] * V(i,j,0) + c[d][1] * V(i,j,1);
                    N_eq = ( rho / 6.0 ) * ( 1.0 + 2.0 * uc + 4.0 * qi_uaub );
                    delta_N = lambda * ( D_now(i,j,d) - N_eq );

                    /* External force at j=1 */
                    if ( j==1 && MAP(i,j) == FLUID )
                        delta_N +=
                            (1.0/3.0) * (c[d][0]*force[0]+c[d][1]*force[1]);

                    switch ( MAP(i,j) )
                    {
                        /* Redistribute fluid according to density/velocity */
                        case FLUID:
                            D_nxt(i,j,d) = D_now(i,j,d) + delta_N;
                            break;
                        /* Walls reflect incoming mass in opposite direction */
                        case WALL:
                            D_nxt(i,j,(d+3)%6) = D_now(i,j,d);
                            break;
                        /* No work to do on solid points */
                        case SOLID:
                            break;
                    }
                }
            }
        }
    }
}


void
stream ( void )
{
    for ( int_t i=0; i<height/size; i++ )
    {
        for ( int_t j=0; j<width; j++ )
        {
            for ( int_t d=0; d<6; d++ )
            {
                /* Compute neighbor indices, wrap around edges */
                int_t
                    ni = ( i + offsets[(i%2)][d][0] + (height/size) ) % (height/size),
                    nj = ( j + offsets[(i%2)][d][1] + width ) % width;

                /* Propagate present fluid density to neighbor */
                D_now(ni,nj,d) = D_nxt(i,j,d);
            }
        }
    }
}


/* Save a snapshot of the present velocity field in 'data/xxxxx.dat',
 * in gnuplot binary matrix format. Indices are given by iteration
 * count / snapshot frequency.
 */
void
store ( int_t iter )
{
    for ( int_t i=0; i<height/size; i++ )
        for ( int_t j=0; j<width; j++ )
            V_abs(i,j) = sqrt ( V(i,j,0)*V(i,j,0) + V(i,j,1)*V(i,j,1) );

    char filename[256];
    memset ( filename, 0, 256 );
    sprintf ( filename, "data/%.5" INT_FMT ".dat", iter );

    MPI_File out = NULL;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);

    // Write data to buffer for writing
    if (rank == 0) {
        Buf(0, 0) = (float) width;
        for ( int_t j=0; j<width; j++ ) {
            Buf(0, j+1) = (float) j;
        }
        MPI_File_write_at(out, 0, buf, width+1, MPI_FLOAT, MPI_STATUS_IGNORE);
    }

    int_t row_offset = rank * height / size;
    for (int_t i = 0; i<height/size; i++) {
        Buf(i+1, 0) = (float) i + row_offset;
        for (int j = 0; j < width; ++j) {
            Buf(i, j+1) = (float) V_abs(i, j);
        }
    }

    int_t buf_len = (width+1) * height/size;
    MPI_Offset offset = (rank * buf_len + width+1) * sizeof(float);

    MPI_File_write_at_all(out, offset, &Buf(1, 0), buf_len, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&out);
}


/* Configure the domain as a channel with a solid circular obstruction
 * positioned at 1/4 of its length.
 */
void
initialize_domain ( void )
{
    int_t center[2] = { (real_t)(height/2.0), (real_t)(width/4.0) };

    /* If no valid value was set from cmd. line, compute obstruction radius */
    if ( radius <= 0.0 )
        radius = height / 20.0;

    /* Bottom wall */
    if (rank == 0) {
        for ( int_t j=0; j<width; j++ ) {
            MAP(0, j) = WALL;
        }
    }
    /* Top wall */
    if (rank == size - 1) {
        for ( int_t j=0; j<width; j++ ) {
            MAP(height/size-1,j) = WALL;
        }
    }

    /* Solid cylinder */
    int_t grid_offset = rank * height / size;
    for ( int_t i=0; i<height/size; i++ )
    {
        for ( int_t j=0; j<width; j++ )
        {
            if ( radius >
                sqrt( (i+grid_offset-center[0])*(i+grid_offset-center[0])+(j-center[1])*(j-center[1]) )
            )
                MAP(i,j) = SOLID;
            else
                MAP(i,j) = FLUID;
        }
    }

    /* Check all solid points for fluid neighbors, and classify them as
     * WALL if they have any.
     */
    for ( int_t i=0; i<height/size; i++ )
    {
        for ( int_t j=0; j<width; j++ )
        {
            if ( MAP(i,j) == SOLID )
            {
                for ( int_t d=0; d<6; d++ )
                {
                    int_t
                        ni = ( i + offsets[(i%2)][d][0] + height/size ) % (height/size),
                        nj = ( j + offsets[(i%2)][d][1] + width ) % width;

                    if ( MAP(ni,nj) == FLUID )
                        MAP(i,j) = WALL;
                }
            }
        }
    }

}


/* Command line options and usage information */
static const char *usage =
"    -h           : Usage information\n"
"    -q           : Suppress input/output\n"
"    -s <integer> : Save state every <integer> steps (default 100)\n"
"    -H <integer> : Set domain height (default 400)\n"
"    -W <integer> : Set domain width  (default 600)\n"
"    -R <integer> : Set obstruction radius (default height/20)\n"
"    -F <real>    : Set horizontal force magnitude (default 1e-2)\n"
"    -I <integer> : Set iteration count (default 40001)\n"
;

void
options ( int argc, char **argv )
{
    int o;
    while ( (o = getopt(argc, argv, "qhs:H:W:R:F:I:")) != -1 )
    {
        switch ( o )
        {
            case 'q': quiet = true; break;
            case 's': snap_freq = strtol ( optarg, NULL, 10 ); break;
            case 'H': height = strtol ( optarg, NULL, 10 ); break;
            case 'W': width  = strtol ( optarg, NULL, 10 ); break;
            case 'R': radius = (real_t)strtol ( optarg, NULL, 10 ); break;
            case 'F': force[1] = (real_t) strtod ( optarg, NULL ); break;
            case 'I': max_iter = strtol ( optarg, NULL, 10 ); break;
            case 'h':
                fprintf ( stderr, "%s\n%s", argv[0], usage );
                exit ( EXIT_FAILURE );
                break;
            default:
                exit ( EXIT_FAILURE );
                break;
        }
    }
}
