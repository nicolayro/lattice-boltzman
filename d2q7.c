#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

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


void init_domain(void);
void collide(void);
void stream(void);

void save(int iteration);

void options(int argc, char **argv);

int_t W, H;
int_t timesteps;
char *input = NULL;

domain_t *lattice = NULL;
double *densities[2] = { NULL, NULL };
double *v = NULL;
double *v_abs = NULL;
double e[6][2];

double force[2] = {
    0.00, // y
    0.05  // x
};

#define LATTICE(i,j) lattice[(i)*W+(j)]

#define D_now(i,j,d) densities[0][(d)*W*H+(i)*W+(j)]
#define D_nxt(i,j,d) densities[1][(d)*W*H+(i)*W+(j)]

#define V_y(i,j) v[2*((i)*W+(j))]
#define V_x(i,j) v[2*((i)*W+(j))+1]
#define V_abs(i,j) v_abs[(i)*W+(j)]

int main(int argc, char **argv)
{
    options(argc, argv);

    init_domain();

    densities[0] = malloc((DIRECTIONS+1) * W * H * sizeof(double));
    densities[1] = malloc((DIRECTIONS+1) * W * H * sizeof(double));
    v = malloc(2 * H * W * sizeof(double));
    v_abs = malloc(H * W * sizeof(double));

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            for (int d = 0; d < DIRECTIONS+1; d++) {
                D_nxt(i,j,d) = D_now(i,j,d) = 1.0 / 7.0;
            }
        }
    }

    for( int_t d=0; d<DIRECTIONS; d++ )
    {
        e[d][0] = sin ( M_PI * d / 3.0 ); // y
        e[d][1] = cos ( M_PI * d / 3.0 ); // x
    }

    for (int_t i = 0; i < timesteps; i++) {
        collide();
        stream();
        if (i % 100 == 0) {
            printf("Iteration %lld/%lld\n", i, timesteps);
            save(i/100);
        }
    }

    free(lattice);
    free(densities[0]);
    free(densities[1]);
    free(v);
    free(v_abs);

    return EXIT_SUCCESS;
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

    lattice = malloc(H * W * sizeof(domain_t));
    if (!lattice) {
        fprintf(stderr, "ERROR: Not enough memory...\n");
        exit(EXIT_FAILURE);
    }

    for (int_t i = 0; i < H; i++) {
        for (int_t j = 0; j < W; j++) {
            int16_t value = geometry[3*(i*W+j)] + geometry[3*(i*W+j)+1]
                + geometry[3*(i*W+j)+2];

            LATTICE(i,j) = value > 0 ? FLUID : SOLID;
        }
    }

    // All SOLID points that are next to FLUID points are categorized as WALL
    for (int i = 0; i < H; i++ ) {
        for (int j = 0; j < W; j++) {
            if (LATTICE(i,j) != SOLID)
                continue;

            for (int d = 0; d < DIRECTIONS; d++) {
                int_t ni = (i + OFFSETS[i%2][d][0]+H)%H;
                int_t nj = (j + OFFSETS[i%2][d][1]+W)%W;

                if (LATTICE(ni,nj) == FLUID)
                    LATTICE(i,j) = WALL;
            }
        }
    }

    /* Bottom wall */
    for (int j = 0; j < W; j++) {
        LATTICE(0,j) = WALL;
    }
    /* Top wall */
    for (int j = 0; j < W; j++) {
        LATTICE((H-1),j) = WALL;
    }

    free(geometry);
    fclose(file);
}

void collide(void)
{
     double rho      = 0.0;  // Density
     double ev       = 0.0;  // Dot product of e and v;
     double N_eq     = 0.0;  // Equilibrium at i
     double delta_N  = 0.0;  // Change

    for (int_t i = 0; i < H; i++) {
        for (int_t j = 0; j < W; j++) {
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
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
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
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            V_abs(i,j) = sqrt(V_x(i,j)*V_x(i,j) + V_y(i,j)*V_y(i,j));
        }
    }

    FILE *file;
    char filename[256];
    memset(filename, 0, 256);
    sprintf(filename, "data/%05d.dat", iteration);

    file = fopen(filename, "wb");
    if (!file) {
        fprintf(stderr, "WARNING: Unable to open file '%s'. Did not save iteration %d\n", filename, iteration);
        exit(EXIT_FAILURE);
    }

    float f;
    f = W;
    fwrite(&f, sizeof(float), 1, file);
    for (int x = 0; x < W; x++) {
        f = x;
        fwrite(&f, sizeof(float), 1, file);
    }
    for (int i = 0; i < H; i++) {
        f = i;
        fwrite(&f, sizeof(float), 1, file);
        for (int j = 0; j < W; j++) {
            f = V_abs(i,j);
            fwrite(&f, sizeof(float), 1, file);
        }
    }

    fclose(file);

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
