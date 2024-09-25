#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#define DIRECTIONS 6
#define ALPHA 0.0
#define TAU 1

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

int_t N, M;
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
    printf("  Width  (N)    %lld\n", N);
    printf("  Height (M)    %lld\n", M);

    densities[0] = malloc((DIRECTIONS+1) * N * M * sizeof(double));
    densities[1] = malloc((DIRECTIONS+1) * N * M * sizeof(double));
    v = malloc(M * N * sizeof(vec2));
    v_abs = malloc(M * N * sizeof(double));

    if (!densities[0] || !densities[1] || !v || !v_abs) {
        fprintf(stderr, "ERROR: Not enough memory...\n");
        exit(EXIT_FAILURE);
    }

    for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
            for (int d = 0; d < DIRECTIONS+1; d++) {
                densities[1][d*M*N+y*N+x] = densities[0][d*M*N+y*N+x] = 1.0 / (DIRECTIONS+1);
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
        if (i % 1 == 0) {
            printf("Iteration %lld/%lld\n", i, timesteps);
            save(i/1);
            printf("Iteration %lld/%lld stored\n", i, timesteps);
        }
    }

    free(sites);
    free(densities[0]);
    free(densities[1]);
    free(v);
    free(v_abs);

    return EXIT_SUCCESS;
}

void init_domain()
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

    N = width;
    M = height;

    sites = malloc(M * N * sizeof(uint8_t));
    if (!sites) {
        exit(EXIT_FAILURE);
    }

    for (int_t y = 0; y < M; y++) {
        for (int_t x = 0; x < N; x++) {
            int16_t value = geometry[3*(y*width+x)] + geometry[3*(y*width+x)+1]
                + geometry[3*(y*width+x)+2];

            sites[y*N+x] = value > 0 ? FLUID : SOLID;
        }
    }

    // All SOLID points that are next to FLUID points are categorized as WALL
    for (int_t y = 0; y < M; y++ ) {
        for (int_t x = 0; x < N; x++) {
            if (sites[y*N+x] != SOLID)
                continue;

            for (int_t i = 0; i < DIRECTIONS; i++) {
                int_t nx = (x + OFFSETS[y%2][i][1]+N)%N;
                int_t ny = (y + OFFSETS[y%2][i][0]+M)%M;

                if (sites[ny*N+nx] == FLUID)
                    sites[y*N+x] = WALL;
            }
        }
    }

    FILE *output = fopen("domain.ppm", "wb");
    fprintf(output,
        "P2\n"
        "# Created by GIMP version 2.10.38 PNM plug-in\n"
        "%lld %lld\n"
        "255\n",
        N,
        M
    );

    for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
            fprintf(output, "%d ", sites[y*N+x]*100);
        }
        fprintf(output, "\n");
    }

    fclose(output);

    free(geometry);
    fclose(file);
}

void collide(void)
{
    /* int    i        = 0;    // Lattice site */
    /* double rho      = 0.0;  // Density */
    /* double ev       = 0.0;  // Dot product of e and v; */
    /* double N_eq     = 0.0;  // Equilibrium */

    for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
            int i = y*N+x;

            // Ignore solid sites
            if (sites[i] == SOLID)
                continue;

            double N_eq, ev;
            double rho = 0;
            v[i].x = v[i].y = 0.0;
            if (sites[i] == FLUID) {
                for (int d = 0; d < DIRECTIONS; d++) {
                    rho += densities[0][d*M*N+y*N+x];
                    v[i].x += e[d].x * densities[0][d*M*N+y*N+x];
                    v[i].y += e[d].y * densities[0][d*M*N+y*N+x];
                }
                v[i].x /= rho;
                v[i].y /= rho;
            }

            for (int d = 0; d < DIRECTIONS; d++) {
                ev = e[d].x * v[i].x + e[d].y * v[i].y;

                // Boundary condition: Reflect of walls
                if (sites[i] == WALL) {
                    densities[1][((d+3)%6)*M*N+y*N+x] = densities[0][d*M*N+y*N+x];
                    continue;
                }

                N_eq =
                    rho*(1.0-ALPHA)/6.0
                    + rho/3.0*ev 
                    + 2.0/3.0*rho*ev*ev 
                    - rho/6.0*v[i].x*v[i].x+v[i].y*v[i].y;

                densities[1][d*M*N+y*N+x] = densities[0][d*M*N+y*N+x] - (densities[0][d*M*N+y*N+x]-N_eq)/TAU;

                // Add external force
                if (x == 0)
                    densities[1][d*M*N+y*N+x] += force.x * e[d].x;
                    /* densities[1][d*M*N+y*N+x] += force.x * (3.0/2.0)*(y*(M-y))/((M/2.0)*(M/2.0)); */
            }

            // Central point / Rest Particle
            /* N_eq = ALPHA*rho - rho*(v[i].x*v[i].x + v[i].y*v[i].y); */
            /* densities[1][6*N*M+i] = densities[0][6*N*M+i] - (N_eq - densities[0][6*N*M+i])/TAU; */
        }
    }
}


void stream(void)
{
    for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
            for (int d = 0; d < DIRECTIONS; d++) {
                int nx = (x + OFFSETS[y%2][d][1]+N)%N;
                int ny = (y + OFFSETS[y%2][d][0]+M)%M;

                densities[0][d*M*N+ny*N+nx] = densities[1][d*M*N+y*N+x];
            }
        }
    }
}

void save(int iteration)
{
    for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
            int i = y * N + x;
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
    f = N;
    fwrite(&f, sizeof(float), 1, file);
    for (int x = 0; x < N; x++) {
        f = x;
        fwrite(&f, sizeof(float), 1, file);
    }
    for (int y = 0; y < M; y++) {
        f = y;
        fwrite(&f, sizeof(float), 1, file);
        for (int x = 0; x < N; x++) {
            f = v_abs[y*N+x];
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
