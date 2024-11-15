#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#define DIRECTIONS 6
#define ALPHA 0.142
#define TAU 1.0

typedef int64_t int_t;

typedef enum {
    SOLID,
    WALL,
    FLUID
} domain_t;

int OFFSETS[2][6][2] = {
    { {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0} }, /* Even rows */
    { {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1} }  /* Odd rows */
};

void init_domain(void);     // Initialize domain geometry from input file
void init_cart_grid(void);  // Initialize MPI cartesian grid
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
double e[6][2];           // Directinal vectors

double force[2] = {
    0.00, // External force in y direction
    0.05  // External force in x direction
};

float *outbuf = NULL; // Output buffer (Note that this is a float)

#define LATTICE(i,j) lattice[(i)*W+(j)]

#define D_now(i,j,d) densities[0][(d)*W*H+(i)*W+(j)]
#define D_nxt(i,j,d) densities[1][(d)*W*H+(i)*W+(j)]

#define V_y(i,j) v[2*((i)*W+(j))]
#define V_x(i,j) v[2*((i)*W+(j))+1]

#define OUTBUF(i,j) outbuf[(i)*W+(j)]


/* MPI */

int rank; // MPI rank
int comm_size; // Total number of ranks

typedef enum { NORTH, EAST, SOUTH, WEST } Direction;
MPI_Comm comm_cart; // Cartesian communicator
int dims[2];        // Dimensions of cartesian grid [y, x]
int cart_pos[2];    // Position in cartesian grid   [y, x]
int cart_nbo[4];    // Neighbors in grid            [N, E, S, W]
MPI_Datatype subgrid;

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

    printf("[%d] Broadcasting\n", rank);
    MPI_Bcast(&W, 1, MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&H, 1, MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&timesteps, 1, MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);

    printf("[%d] Creating cartesian grid\n", rank);
    int periods[2] = { 0, 0 };
    MPI_Dims_create(comm_size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);
    MPI_Cart_coords(comm_cart, rank, 2, cart_pos);
    MPI_Cart_shift(comm_cart, 1, 1, &cart_nbo[WEST], &cart_nbo[EAST]);
    MPI_Cart_shift(comm_cart, 0, 1, &cart_nbo[SOUTH], &cart_nbo[NORTH]);
    local_H = H / dims[0]; // TODO: This or opposite?
    local_W = W / dims[1];

    printf("[%d] Creating datatypes\n", rank);

    int start[2] = { 0, 0 };
    int subgrid_size[2] = { local_H, local_W };
    int grid_size[2] = { H, W };

    MPI_Type_create_subarray(2, grid_size, subgrid_size, start, MPI_ORDER_C, MPI_INT, &subgrid);
    MPI_Type_commit(&subgrid);

    if (rank == MPI_RANK_ROOT)
        free(domain);

    printf("[%d] Allocating\n", rank);
    lattice = malloc((local_W+2) * (local_H+2) * sizeof(domain_t));
    densities[0] = malloc(7 * (local_W+2) * (local_H+2) * sizeof(double));
    densities[1] = malloc(7 * (local_W+2) * (local_H+2) * sizeof(double));
    v = malloc(2 * local_H * local_W * sizeof(double));
    outbuf = malloc((local_H+1) * (local_W+1) * sizeof(float));

    printf("[%d] Scattering\n", rank);

    printf("[%d] Init densities\n", rank);
    for (int i = 0; i < local_H+2; i++) {
        for (int j = 0; j < local_W+2; j++) {
            for (int d = 0; d < 7; d++) {
                D_nxt(i,j,d) = D_now(i,j,d) = 1.0 / 7.0;
            }
        }
    }

    printf("[%d] Init directional vectors\n", rank);
    for(int_t d=0; d<6; d++) {
        e[d][0] = sin(M_PI * d / 3.0); // y
        e[d][1] = cos(M_PI * d / 3.0); // x
    }

    for (int_t i = 0; i < timesteps; i++) {
        printf("[%d] Colliding\n", rank);
        collide();
        border_exchange();
        printf("[%d] Streaming\n", rank);
        stream();

        if (i % 1 == 0) {
            printf("Iteration %lld/%lld\n", i, timesteps);
            save(i/1);
        }
    }

    free(lattice);
    free(densities[0]);
    free(densities[1]);
    free(v);
    free(outbuf);

    MPI_Finalize();

    return EXIT_SUCCESS;
}

void init_cart_grid(void) {
    // TODO: Offsets?
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

            for (int d = 0; d < DIRECTIONS; d++) {
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

void border_exchange(void) {
    // Setup column and row datatypes for easier border exchange
    MPI_Datatype column, row;
    MPI_Type_vector(local_H, 1, local_W + 2, MPI_DOUBLE, &column);
    MPI_Type_vector(local_W, 1, 1, MPI_DOUBLE, &row);

    MPI_Type_commit (&column);
    MPI_Type_commit (&row);

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
        MPI_Sendrecv(&D_nxt(0, 1, d), 1, column, cart_nbo[SOUTH], d+6,
                     &D_nxt(0, 0, d), 1, column, cart_nbo[NORTH], d+6,
                     comm_cart, MPI_STATUS_IGNORE);
    }

    // West
    for (int_t d = 0; d < 6; ++d) {
        MPI_Sendrecv(&D_nxt(local_W, 0, d), 1, column, cart_nbo[SOUTH], d+6,
                     &D_nxt(0, 0, d), 1, column, cart_nbo[NORTH], d+6,
                     comm_cart, MPI_STATUS_IGNORE);
    }
}

void collide(void)
{
     double rho      = 0.0;  // Density
     double ev       = 0.0;  // Dot product of e and v;
     double N_eq     = 0.0;  // Equilibrium at i
     double delta_N  = 0.0;  // Change

    for (int_t i = 1; i <= local_H; i++) {
        for (int_t j = 1; j <= local_W; j++) {
            if (LATTICE(i,j) != WALL && LATTICE(i,j) != SOLID && LATTICE(i,j) != FLUID) {
                printf("[%d] EEEEEEH (%lld, %lld): %d\n", rank, i, j, LATTICE(i,j));
            }
            assert(LATTICE(i,j) == WALL || LATTICE(i,j) == SOLID || LATTICE(i,j) == FLUID);

            // Ignore solid sites
            if (LATTICE(i,j) == SOLID)
                continue;

            rho = D_now(i,j,6);
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
                    D_nxt(i,j,(d+3)%6) = D_now(i,j,d);
                    continue;
                }

                ev = e[d][1] * V_x(i,j) + e[d][0] * V_y(i,j);
                N_eq =
                    // F_eq_i
                    rho*(1.0-ALPHA)/6.0
                    + rho/3.0*ev
                    + (2.0*rho/3.0)*ev*ev
                    - rho/6.0*(V_x(i,j)*V_x(i,j)+V_y(i,j)*V_y(i,j))
                    // F_eq_0
                    + (ALPHA*rho/6.0 - rho/6.0*(V_x(i,j)*V_x(i,j) + V_y(i,j)*V_y(i,j)));

                delta_N = -(D_now(i,j,d)-N_eq)/TAU;
                // Add external force
                if (j == 1)
                    delta_N += (1.0/3.0) * force[1] * e[d][1];

                D_nxt(i,j,d) = D_now(i,j,d) + delta_N;
            }

            if (LATTICE(i,j) == FLUID) {
                N_eq = ALPHA*rho - rho*(V_x(i,j)*V_x(i,j) + V_y(i,j)*V_y(i,j));
                D_nxt(i,j,6) = D_now(i,j,6) - (D_now(i,j,6)-N_eq)/TAU;
            }
        }
    }
}

void stream(void)
{
    for (int i = 0; i < local_H+2; i++) {
        for (int j = 0; j < local_W+2; j++) {
            for (int d = 0; d < DIRECTIONS; d++) {
                int ni = (i + OFFSETS[i%2][d][0]+H)%H;
                int nj = (j + OFFSETS[i%2][d][1]+W)%W;

                D_now(ni,nj,d) = D_nxt(i,j,d);
            }
        }
    }
}

void save(int iteration)
{
    // Write data to output buffer
    if (rank == MPI_RANK_ROOT) {
        OUTBUF(0, 0) = (float) W;
        for (int_t j = 0; j < W; j++) {
            OUTBUF(0, j+1) = (float) j;
        }
    }
    for (int_t i = 1; i<=local_H; i++) {
        OUTBUF(i, 0) = (float)(i - 1 + rank*local_H);
        for (int j = 1; j <= local_W; ++j) {
            OUTBUF(i, j) = (float) sqrt (V_y(i,j)*V_y(i,j) + V_x(i,j)*V_x(i,j));
        }
    }

    char filename[256];
    memset(filename, 0, 256);
    sprintf(filename, "data/%05d.dat", iteration);

    MPI_File output;
    MPI_File_open(comm_cart, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &output);
    if (!output) {
        fprintf(stderr, "WARNING: Unable to open file '%s'. Did not save iteration %d\n", filename, iteration);
        exit(EXIT_FAILURE);
    }

    if (rank == MPI_RANK_ROOT) {
        MPI_File_write_at(output, 0, outbuf, W+1, MPI_FLOAT, MPI_STATUS_IGNORE);
    }

    MPI_Offset offset = (rank*local_H+1) * (rank*local_W+1) * sizeof(float);

    MPI_File_set_view(output, offset, MPI_FLOAT, subgrid, "native", MPI_INFO_NULL);
    MPI_File_write_all(output, outbuf, 1, subgrid, MPI_STATUS_IGNORE);

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
