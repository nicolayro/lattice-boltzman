#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#define DIRECTIONS 6
#define ALPHA 0.17
#define TAU 1.0

typedef int64_t int_t;

typedef enum {
    SOLID,
    WALL,
    FLUID
} domain_t;

typedef struct {
    double x, y;
} vec2;

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

domain_t *sites = NULL;
double *densities[2] = { NULL, NULL };
vec2 *v = NULL;
double *v_abs = NULL;
vec2 e[6];

vec2 force = {
    .x = 0.01,
    .y = 0.00
};

int main(int argc, char **argv)
{
    options(argc, argv);

    init_domain();

    printf("=== Domain ===\n");
    printf("  Width  (N)    %lld\n", W);
    printf("  Height (M)    %lld\n", H);

    densities[0] = malloc((DIRECTIONS+1) * W * H * sizeof(double));
    densities[1] = malloc((DIRECTIONS+1) * W * H * sizeof(double));
    v = malloc(H * W * sizeof(vec2));
    v_abs = malloc(H * W * sizeof(double));
    if (!densities[0] || !densities[1] || !v || !v_abs) {
        fprintf(stderr, "ERROR: Not enough memory...\n");
        exit(EXIT_FAILURE);
    }

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            for (int d = 0; d < DIRECTIONS+1; d++) {
                densities[1][d*H*W+y*W+x] = densities[0][d*H*W+y*W+x] = 1.0 / 6.0;
            }
        }
    }

    for( int_t d=0; d<DIRECTIONS; d++ )
    {
        e[d].x = cos ( M_PI * d / 3.0 );
        e[d].y = sin ( M_PI * d / 3.0 );
    }

    for (int_t i = 0; i < timesteps; i++) {
        collide();
        stream();
        if (i % 100 == 0) {
            printf("Iteration %lld/%lld\n", i, timesteps);
            save(i/100);
        }
    }

    free(sites);
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

    sites = malloc(H * W * sizeof(domain_t));
    if (!sites) {
        exit(EXIT_FAILURE);
    }

    for (int_t y = 0; y < H; y++) {
        for (int_t x = 0; x < W; x++) {
            int16_t value = geometry[3*(y*W+x)] + geometry[3*(y*W+x)+1]
                + geometry[3*(y*W+x)+2];

            sites[y*W+x] = value > 0 ? FLUID : SOLID;
        }
    }

    // All SOLID points that are next to FLUID points are categorized as WALL
    for (int y = 0; y < H; y++ ) {
        for (int x = 0; x < W; x++) {
            if (sites[y*W+x] != SOLID)
                continue;

            for (int i = 0; i < DIRECTIONS; i++) {
                int_t ny = (y + OFFSETS[y%2][i][0]+H)%H;
                int_t nx = (x + OFFSETS[y%2][i][1]+W)%W;

                if (sites[ny*W+nx] == FLUID)
                    sites[y*W+x] = WALL;
            }
        }
    }

    /* Bottom wall */
    for (int x = 0; x < W; x++) {
        sites[x] = WALL;
    }
    /* Top wall */
    for (int x = 0; x < W; x++) {
        sites[(H-1)*W+x] = WALL;
    }

    free(geometry);
    fclose(file);
}

void collide(void)
{
     int    i        = 0;    // Lattice site
     double rho      = 0.0;  // Density
     double ev       = 0.0;  // Dot product of e and v;
     double N_eq     = 0.0;  // Equilibrium at i
     double delta_N  = 0.0;  // Change

    for (int_t y = 0; y < H; y++) {
        for (int_t x = 0; x < W; x++) {
            assert(sites[y*W+x] == WALL || sites[y*W+x] == SOLID || sites[y*W+x] == FLUID);

            i = y*W+x;
            // Ignore solid sites
            if (sites[i] == SOLID)
                continue;

            /*rho = densities[0][6*H*W+y*W+x];*/

            rho = 0.0;
            v[i].x = v[i].y = 0.0;
            if (sites[i] == FLUID) {
                for (int d = 0; d < DIRECTIONS; d++) {
                    rho += densities[0][d*H*W+y*W+x];
                    v[i].x += e[d].x * densities[0][d*H*W+y*W+x];
                    v[i].y += e[d].y * densities[0][d*H*W+y*W+x];
                }
                assert(rho != 0.0);
                v[i].x /= rho;
                v[i].y /= rho;
            }

            for (int d = 0; d < DIRECTIONS; d++) {
                // Boundary condition: Reflect of walls
                if (sites[i] == WALL) {
                    densities[1][((d+3)%6)*H*W+y*W+x] = densities[0][d*H*W+y*W+x];
                    continue;
                }

                ev = e[d].x * v[i].x + e[d].y * v[i].y;
                N_eq =
                    // F_eq_i
                    rho*(1.0-ALPHA)/6.0
                    + rho/3.0*ev
                    + (2.0*rho/3.0)*ev*ev
                    - rho/6.0*(v[i].x*v[i].x+v[i].y*v[i].y)
                    // F_eq_0
                    + (ALPHA*rho/6.0 - rho/6.0*(v[i].x*v[i].x + v[i].y*v[i].y));

                delta_N = -(densities[0][d*H*W+y*W+x]-N_eq)/TAU;
                // Add external force
                if (x == 1)
                    delta_N += (1.0/3.0) * force.x * e[d].x;

                densities[1][d*H*W+y*W+x] = densities[0][d*H*W+y*W+x] + delta_N;
            }

            // Central point / Rest Particle
            /*if (sites[i] == FLUID) {*/
            /*    N_eq = ALPHA*rho - rho*(v[i].x*v[i].x + v[i].y*v[i].y);*/
            /*    densities[1][6*W*H+i] = densities[0][6*W*H+i] - (densities[0][6*W*H+i]-N_eq)/TAU;*/
            /*}*/
        }
    }
}

void stream(void)
{
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            for (int d = 0; d < DIRECTIONS; d++) {
                int ny = (y + OFFSETS[y%2][d][0]+H)%H;
                int nx = (x + OFFSETS[y%2][d][1]+W)%W;

                densities[0][d*H*W+ny*W+nx] = densities[1][d*H*W+y*W+x];
            }
        }
    }
}

void save(int iteration)
{
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int i = y * W + x;
            v_abs[i] = sqrt(v[i].x*v[i].x + v[i].y*v[i].y);
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
    for (int y = 0; y < H; y++) {
        f = y;
        fwrite(&f, sizeof(float), 1, file);
        for (int x = 0; x < W; x++) {
            f = v_abs[y*W+x];
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
