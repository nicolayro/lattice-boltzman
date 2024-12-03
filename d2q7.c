#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#define DIRECTIONS 7
#define ALPHA 0.5
#define TAU 1.0

typedef int64_t int_t;

typedef enum {
    SOLID,
    WALL,
    FLUID
} domain_t;

int OFFSETS[2][DIRECTIONS][2] = {
    { {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1}, {0,0} },  /* Odd rows */
    { {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0}, {0,0} } /* Even rows */
};

void init_domain(void);     // Initialize domain geometry from input file
void init_cart_grid(void);       // Initialize MPI Cartesian grid
void init_types(void);      // Initialize MPI Custom datatypes
void scatter_domain(void);
void collide(void);         // Collision step
void border_exchange(void); // MPI border exchange
void stream(void);          // Streaming step

void save(int iteration);

void options(int argc, char **argv);

char *input = NULL; // Input file (in raw .ppm format)
domain_t *domain = NULL;

int_t W, H;      // Width and height of domain
int_t timesteps; // Number of timesteps in simulation

domain_t *lattice = NULL; // Domain geometry
double *densities[2] = {
    NULL,                 // Densities in current timestep
    NULL                  // Densities in next timestep
};
double *v = NULL;         // Velocities
double e[DIRECTIONS][2];           // Directinal vectors

double force[2] = {
    0.00, // External force in y direction
    0.01  // External force in x direction
};

float *outbuf = NULL; // Output buffer (Note that this is a float)

#define LATTICE(i,j) lattice[(i)*(local_W+2)+(j)]

#define D_now(i,j,d) densities[0][(d)*(local_W+2)*(local_H+2)+(i)*(local_W+2)+(j)]
#define D_nxt(i,j,d) densities[1][(d)*(local_W+2)*(local_H+2)+(i)*(local_W+2)+(j)]

#define V_y(i,j) v[2*((i)*(local_W+2)+(j))]
#define V_x(i,j) v[2*((i)*(local_W+2)+(j))+1]

#define OUTBUF(i,j) outbuf[(i)*(local_W)+(j)]

/* MPI */

int rank; // MPI rank
int comm_size; // Total number of ranks

typedef enum { NORTH, EAST, SOUTH, WEST } Direction;
MPI_Comm comm_cart; // Cartesian communicator
int dims[2];        // Dimensions of cartesian grid [y, x]
int cart_pos[2];    // Position in cartesian grid   [y, x]
int cart_nbo[4];    // Neighbors in grid            [N, E, S, W]
MPI_Datatype subgrid;       // Datatype for local grid in global grid
MPI_Datatype column, row;   // Column and row in subgrid (including halo)

int local_H, local_W, local_x_offsett;


#define MPI_RANK_ROOT 0

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (rank == MPI_RANK_ROOT) {
        options(argc, argv);
        printf("Initializing domain\n");
        init_domain();
    }

    MPI_Bcast(&W, 1, MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&H, 1, MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&timesteps, 1, MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);

    init_cart_grid();

    lattice = malloc((local_W+2) * (local_H+2) * sizeof(domain_t));
    densities[0] = malloc(7 * (local_W+2) * (local_H+2) * sizeof(double));
    densities[1] = malloc(7 * (local_W+2) * (local_H+2) * sizeof(double));
    v = malloc(2 * (local_H+2) * (local_W+2) * sizeof(double));
    outbuf = malloc((local_H) * (local_W) * sizeof(float));

    scatter_domain();

    for (int i = 0; i < local_H+2; i++) {
        for (int j = 0; j < local_W+2; j++) {
            for (int d = 0; d < DIRECTIONS; d++) {
                D_nxt(i,j,d) = D_now(i,j,d) = 1.0 / 7.0;
            }
        }
    }

    for(int_t d=0; d<6; d++) {
        e[d][0] = sin(M_PI * d / 3.0); // y
        e[d][1] = cos(M_PI * d / 3.0); // x
    }
    e[6][0] = 0.0;
    e[6][1] = 0.0;


    init_types();

    for (int_t i = 0; i < timesteps; i++) {
        collide();
        border_exchange();
        stream();

        if (i % 100 == 0) {
            if (rank == MPI_RANK_ROOT)
                printf("Iteration %lld/%lld\n", i, timesteps);
            save(i/100);
        }
    }

    MPI_Type_free(&column);
    MPI_Type_free(&row);
    MPI_Type_free(&subgrid);

    free(lattice);
    free(densities[0]);
    free(densities[1]);
    free(v);
    free(outbuf);

    MPI_Finalize();

    return EXIT_SUCCESS;
}

void init_cart_grid(void)
{
    int periods[2] = { 1, 1 };

    MPI_Dims_create(comm_size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);
    MPI_Cart_coords(comm_cart, rank, 2, cart_pos);
    MPI_Cart_shift(comm_cart, 0, 1, &cart_nbo[NORTH], &cart_nbo[SOUTH]);
    MPI_Cart_shift(comm_cart, 1, 1, &cart_nbo[WEST], &cart_nbo[EAST]);

    local_H = H / dims[0];
    local_W = W / dims[1];
}

void init_types(void)
{
    int start[2] = { cart_pos[0]*local_H, cart_pos[1]*local_W };
    int subgrid_size[2] = { local_H, local_W };
    int grid_size[2] = { H, W };

    MPI_Type_create_subarray(2, grid_size, subgrid_size, start, MPI_ORDER_C, MPI_FLOAT, &subgrid);
    MPI_Type_commit(&subgrid);

    MPI_Type_vector(local_H+2, 1, local_W+2, MPI_DOUBLE, &column);
    MPI_Type_vector(1, local_W+2, local_W+2, MPI_DOUBLE, &row);

    MPI_Type_commit (&column);
    MPI_Type_commit (&row);

}

void init_domain(void)
{
    FILE *file = fopen(input, "r");
    if (!file) {
        fprintf(stderr, "ERROR: Unable to open file %s.\n", input);
        exit(EXIT_FAILURE);
    }

    char magic[4];
    if (!fgets(magic, sizeof(magic), file)) {
        fprintf(stderr, "ERROR: Unable to read magic.\n");
        exit(EXIT_FAILURE);
    }
    if (magic[0] != 'P' || magic[1] != '6') {
        fprintf(stderr, "ERROR: Invalid magic '%.*s'. Only P6 is supported.\n", 3, magic);
        exit(EXIT_FAILURE);
    }

    // Comments
    int c;
    while ((c = fgetc(file)) == '#') {
        while ((c = fgetc(file)) != '\n');
    }
    ungetc(c, file);

    // Data size
    int height, width;
    if (fscanf(file, "%d %d", &width, &height) != 2) {
        fprintf(stderr, "ERROR: Invalid format for width and height.\n");
        exit(EXIT_FAILURE);
    }
    if (height < 0 || height > 4096) {
        fprintf(stderr, "ERROR: Unsupported height %d\n", height);
        exit(EXIT_FAILURE);
    }
    if (width < 0 || width > 4096) {
        fprintf(stderr, "ERROR: Unsupported width %d\n", width);
        exit(EXIT_FAILURE);
    }

    // Pixel size / depth
    int pixel_size;
    if (fscanf(file, "%d", &pixel_size) != 1) {
        fprintf(stderr, "ERROR: Invalid format for pixel size.\n");
        exit(EXIT_FAILURE);
    }
    if (pixel_size != 255) {
        fprintf(stderr, "ERROR: Unsupported pixel size %d\n", pixel_size);
        exit(EXIT_FAILURE);
    }

    uint8_t *geometry = malloc(height * width * 3 * sizeof(uint8_t));
    if (!geometry) {
        fprintf(stderr, "ERROR: Not enough memory...\n");
        exit(EXIT_FAILURE);
    }

    if ((c = fgetc(file)) != '\n') {
        ungetc(c, file);
    }

    int read = fread(geometry, sizeof(uint8_t), height * width * 3, file);
    if (read != height*width*3) {
        fprintf(stderr, "ERROR: Unable to read file.\n"
                "only %d/%d objects were read.\n", read, height*width*3);
        exit(EXIT_FAILURE);
    }

    W = width;
    H = height;

    domain = malloc(H * W * sizeof(domain_t));
    if (!domain) {
        fprintf(stderr, "ERROR: Not enough memory...\n");
        exit(EXIT_FAILURE);
    }

    for (int_t i = 0; i < H; i++) {
        for (int_t j = 0; j < W; j++) {
            int16_t value = geometry[3*(i*W+j)] + geometry[3*(i*W+j)+1]
                + geometry[3*(i*W+j)+2];

            domain[i*W+j] = value > 0 ? FLUID : SOLID;
        }
    }

    // All SOLID points that are next to FLUID points are categorized as WALL
    for (int i = 0; i < H; i++ ) {
        for (int j = 0; j < W; j++) {
            if (domain[i*W+j] != SOLID)
                continue;

            for (int d = 0; d < DIRECTIONS-1; d++) {
                int_t ni = (i + OFFSETS[i%2][d][0]+H)%H;
                int_t nj = (j + OFFSETS[i%2][d][1]+W)%W;

                if (domain[ni*W+nj] == FLUID)
                    domain[i*W+j] = WALL;
            }
        }
    }

    /* Bottom wall */
    for (int j = 0; j < W; j++) {
        domain[j] = WALL;
    }
    /* Top wall */
    for (int j = 0; j < W; j++) {
        domain[(H-1)*W+j] = WALL;
    }

    free(geometry);
    fclose(file);
}

void scatter_domain(void)
{
    if (rank == MPI_RANK_ROOT) {
        // Scatter to others
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                if (i*dims[1]+j==MPI_RANK_ROOT)
                    continue;

                MPI_Datatype map_subgrid;

                int start[2] = { i*local_H, j*local_W };
                int subgrid_size[2] = { local_H, local_W };
                int grid_size[2] = { H, W };

                MPI_Type_create_subarray(2, grid_size, subgrid_size, start, MPI_ORDER_C, MPI_INT32_T, &map_subgrid);
                MPI_Type_commit(&map_subgrid);

                MPI_Send(domain, 1, map_subgrid, i*dims[1]+j, 1, comm_cart);

                MPI_Type_free(&map_subgrid);
            }
        }

        for (int i = 0; i < local_H; i++) {
            for (int j = 0; j < local_W; j++) {
                LATTICE(i+1,j+1) = domain[i*W+j];
            }
        }

        free(domain);
    } else {
        MPI_Datatype local_grid;
        MPI_Type_vector(local_H, local_W, local_W+2, MPI_INT32_T, &local_grid);
        MPI_Type_commit(&local_grid);
        MPI_Recv(&LATTICE(1,1), 1, local_grid, MPI_RANK_ROOT, 1, comm_cart, MPI_STATUS_IGNORE);
        MPI_Type_free(&local_grid);
    }

    MPI_Datatype map_column, map_row;
    MPI_Type_vector(local_H+2, 1, local_W + 2, MPI_INT32_T, &map_column);
    MPI_Type_vector(1, local_W+2, local_W + 2, MPI_INT32_T, &map_row);

    MPI_Type_commit (&map_column);
    MPI_Type_commit (&map_row);

    // Send north
    MPI_Sendrecv(&LATTICE(1, 0), 1, map_row, cart_nbo[NORTH], 1,
                 &LATTICE(local_H+1, 0), 1, map_row, cart_nbo[SOUTH], 1,
                 comm_cart, MPI_STATUS_IGNORE);

    // Send south
    MPI_Sendrecv(&LATTICE(local_H, 0), 1, map_row, cart_nbo[SOUTH], 1,
                 &LATTICE(0, 0), 1, map_row, cart_nbo[NORTH], 1,
                 comm_cart, MPI_STATUS_IGNORE);

    // Send west
    MPI_Sendrecv(&LATTICE(0, 1), 1, map_column, cart_nbo[WEST], 1,
                 &LATTICE(0, local_W+1), 1, map_column, cart_nbo[EAST], 1,
                 comm_cart, MPI_STATUS_IGNORE);

    // Send east
    MPI_Sendrecv(&LATTICE(0, local_W), 1, map_column, cart_nbo[EAST], 1,
                 &LATTICE(0, 0), 1, map_column, cart_nbo[WEST], 1,
                 comm_cart, MPI_STATUS_IGNORE);

    MPI_Type_free(&map_column);
    MPI_Type_free(&map_row);
}


void border_exchange(void) {
    // Send north
    for (int_t d = 0; d < 6; ++d) {
        MPI_Sendrecv(&D_nxt(1, 0, d), 1, row, cart_nbo[NORTH], d,
                     &D_nxt(local_H+1, 0, d), 1, row, cart_nbo[SOUTH], d,
                     comm_cart, MPI_STATUS_IGNORE);
    }

    // Send south
    for (int_t d = 0; d < 6; ++d) {
        MPI_Sendrecv(&D_nxt(local_H, 0, d), 1, row, cart_nbo[SOUTH], d+6,
                     &D_nxt(0, 0, d), 1, row, cart_nbo[NORTH], d+6,
                     comm_cart, MPI_STATUS_IGNORE);
    }

    // Send west
    for (int_t d = 0; d < 6; ++d) {
        MPI_Sendrecv(&D_nxt(0, 1, d), 1, column, cart_nbo[WEST], d+12,
                     &D_nxt(0, local_W+1, d), 1, column, cart_nbo[EAST], d+12,
                     comm_cart, MPI_STATUS_IGNORE);
    }

    // Send east
    for (int_t d = 0; d < 6; ++d) {
        MPI_Sendrecv(&D_nxt(0, local_W, d), 1, column, cart_nbo[EAST], d+18,
                     &D_nxt(0, 0, d), 1, column, cart_nbo[WEST], d+18,
                     comm_cart, MPI_STATUS_IGNORE);
    }

}

void collide(void)
{
     double rho      = 0.0;  // Density
     double ev       = 0.0;  // Dot product of e and v;
     double N_eq     = 0.0;  // Equilibrium at i
     double delta_N  = 0.0;  // Change

    for (int i = 1; i <= local_H; i++) {
        for (int j = 1; j <= local_W; j++) {
            if (LATTICE(i,j) != WALL && LATTICE(i,j) != SOLID && LATTICE(i,j) != FLUID) {
                printf("(%d, %d)\n", i, j);
                exit(EXIT_FAILURE);
            }
            assert(LATTICE(i,j) == WALL || LATTICE(i,j) == SOLID || LATTICE(i,j) == FLUID);

            // Ignore solid sites
            if (LATTICE(i,j) == SOLID) {
                continue;
            }

            rho = 0.0;
            V_x(i,j) = V_y(i,j) = 0.0;
            if (LATTICE(i,j) == FLUID) {
                for (int d = 0; d < DIRECTIONS; d++) {
                    rho += D_now(i,j,d);
                    V_y(i,j) += e[d][0] * D_now(i,j,d);
                    V_x(i,j) += e[d][1] * D_now(i,j,d);
                }
                assert(rho != 0.0);
                V_y(i,j) /= rho;
                V_x(i,j) /= rho;
            }

            for (int d = 0; d < DIRECTIONS; d++) {
                // Boundary condition: Reflect of walls
                if (LATTICE(i,j) == WALL) {
                    if (d != 6) {
                        D_nxt(i,j,(d+3)%6) = D_now(i,j,d);
                    }
                    continue;
                }

                ev = e[d][1] * V_x(i,j) + e[d][0] * V_y(i,j);
                if (d == 6) {
                    // Rest particle
                    N_eq = ALPHA*rho - rho*(V_x(i,j)*V_x(i,j) + V_y(i,j)*V_y(i,j));
                } else {
                    // Outgoing vectors
                    N_eq =
                        // F_eq_i
                        rho*(1.0-ALPHA)/6.0
                        + rho/3.0*ev
                        + (2.0*rho/3.0)*ev*ev
                        - rho/6.0*(V_x(i,j)*V_x(i,j)+V_y(i,j)*V_y(i,j));
                }

                delta_N = -(D_now(i,j,d)-N_eq)/TAU;

                if (cart_pos[1] * local_W + j == 1)
                    delta_N += (1.0/3.0) * force[1] * e[d][1];

                D_nxt(i,j,d) = D_now(i,j,d) + delta_N;
            }

        }
    }
}

void stream(void)
{
    for (int i = 0; i < local_H+2; i++) {
        for (int j = 0; j < local_W+2; j++) {
            for (int d = 0; d < DIRECTIONS; d++) {
                int ni = (i + OFFSETS[i%2][d][0]);
                int nj = (j + OFFSETS[i%2][d][1]);

                if (ni < 0 || ni >= local_H+2 || nj < 0 || nj >= local_W+2)
                    continue;

                D_now(ni,nj,d) = D_nxt(i,j,d);
            }
        }
    }
}

void save(int iteration)
{
    // Caculate absolute velocity (without halo)
    for (int i = 1; i <= local_H; i++) {
        for (int j = 1; j <= local_W; j++) {
            OUTBUF((i-1),(j-1)) = (float) sqrt(V_y(i,j)*V_y(i,j) + V_x(i,j)*V_x(i,j));
        }
    }

    char filename[256];
    memset(filename, 0, 256);
    sprintf(filename, "data/%05d.dat", iteration);

    MPI_File output = NULL;
    MPI_File_open(comm_cart, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &output);

    if (!output) {
        fprintf(stderr, "WARNING: Unable to open file '%s'. Did not save iteration %d\n", filename, iteration);
        exit(EXIT_FAILURE);
    }

    MPI_File_set_view(output, 0, MPI_FLOAT, subgrid, "native", MPI_INFO_NULL);
    MPI_File_write_all(output, outbuf, local_H*local_W, MPI_FLOAT, MPI_STATUS_IGNORE);

    MPI_File_close(&output);
}

void options(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "i:h")) != -1 ) {
        switch (c) {
            case 'i':
                timesteps = strtol(optarg, NULL, 10);
                break;
            case 'h':
                printf("Usage: d2q7 [-i iter] file\n");
                printf("  file       set domain from P6 ppm file\n");
                printf("  options\n");
                printf("    -i iter  number of iterations (default 1000)\n");
                printf("    -h       display this message\n");
                exit(EXIT_SUCCESS);
                break;
            default:
                opterr = 0;
                fprintf(stderr, "ERROR: Illegal argument %c\n", c);
                exit(EXIT_FAILURE);
                break;
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "ERROR: Missing file.\n");
        exit(EXIT_FAILURE);
    }

    input = argv[optind];
}
