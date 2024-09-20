#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#define DIRECTIONS 7
#define ALPHA 1

typedef int64_t int_t;

typedef enum {
    SOLID,
    WALL,
    FLUID
} domain_t;

typedef struct {
    double x, y;
} vec2;

double vec2_dot(vec2 a, vec2 b);

int OFFSETS[2][7][2] = {
    { {0, 0}, {0,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1}, {-1,0} }, /* Even rows */
    { {0, 0}, {0,1}, {1,1}, { 1,0}, {0,-1}, {-1, 0}, {-1,1} }  /* Odd rows */
};

void init_domain(void);
void collide(void);
void stream(void);

void options(int argc, char **argv);

int_t N, M;
int_t timesteps;
char *input = NULL;

domain_t *sites = NULL;
double *densities[2] = { NULL, NULL };
vec2 *v;
vec2 e[6];

int main(int argc, char **argv)
{
    options(argc, argv);

    N = 10;
    M = 10;
    timesteps = 1;

    init_domain();

    printf("=== Domain ===\n");
    printf("  Width  (N)    %lld\n", N);
    printf("  Height (M)    %lld\n", M);

    densities[0] = malloc(DIRECTIONS * N * M * sizeof(double));
    densities[1] = malloc(DIRECTIONS * N * M * sizeof(double));

    for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
            for (int d = 0; d < DIRECTIONS; d++) {
                densities[0][d * (y*N+x)] = 1.0 / 7.0;
            }
        }
    }

    for( int_t d=0; d<6; d++ )
    {
        e[d].x = sin ( M_PI * d / 3.0 );
        e[d].y = cos ( M_PI * d / 3.0 );
    }

    for (int_t i = 0; i < timesteps; i++) {
        collide();
        /*stream();*/
    }

    free(sites);
    free(densities[0]);
    free(densities[1]);

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
    for (int_t y = 0; y < M; y++) {
        for (int_t x = 0; x < N; x++) {
            int16_t value = geometry[3*(y*width+x)] + geometry[3*(y*width+x)+1]
                + geometry[3*(y*width+x)+2];

            sites[y*width+x] = value > 0 ? FLUID : SOLID;
        }
    }

    // All SOLID points that are next to FLUID points are categorized as WALL
    for (int_t y = 0; y < M; y++ ) {
        for (int_t x = 0; x < N; x++) {
            if (sites[y*width+x] != SOLID)
                continue;

            for (int_t i = 1; i < DIRECTIONS; i++) {
                int_t nx = x + OFFSETS[y%2][i][0];
                int_t ny = y + OFFSETS[y%2][i][1];

                bool in_bounds = (nx >= 0 || nx < N || ny >= 0 || ny < M);
                if (in_bounds && sites[ny*width+nx] == FLUID)
                    sites[y*width+x] = WALL;
            }
        }
    }

    free(geometry);
    fclose(file);
}

void collide(void)
{
    int i           = 0;    // Lattice site
    double rho      = 0.0;  // Density
    double ev       = 0.0;  // Dot product of e and v;
    double N_eq     = 0.0;  // Equilibrium

    for (int y = 0; y < M; y++) {
        for (int x = 0; x < N; x++) {
            i = y*N+x;

            // Ignore solid sites
            if (sites[i] == SOLID)
                continue;

            rho = 0;
            for (int d = 1; d < DIRECTIONS; d++) {
                rho += densities[0][d*i];
            }
            ev = vec2_dot(e[i], v[i]);
            for (int d = 0; d < DIRECTIONS; d++) {
                // Boundary condition: Reflect of walls
                if (sites[i] == WALL) {
                    densities[1][(d%3)*i] = densities[0][d*i];
                    continue;
                }

                N_eq = rho*(1.0-ALPHA)/6.0 + rho/3.0*ev + 2.0/3.0*rho*ev*ev - rho/6.0*v[i].x*v[i].y;
                switch (sites[d]) {
                    case FLUID:
                        densities[1][d*i] = N_eq - densities[0][d*i];
                        break;
                    case WALL:
                        densities[1][d*i] = N_eq - densities[0][d*i];
                        break;
                    case SOLID: break;
                    default:
                        assert(false && "Unknown domain type");
                        break;
                }
            }
        }
    }
}


/**
 *  * • Collisions. For each lattice site, do the following:
 * 1. From the nine microscopic densities ni
 * , compute the macroscopic density
 * ρ and velocity components ux and uy.
 * 2. From these three macroscopic variables, use equation (8) to compute the
 * equilibrium number densities n
 * eq
 * i
 * .
 * 3. Update each of the nine number densities according to equation (9).
 * • Streaming. Move all the moving molecules into adjacent or diagonal lattice
 * sites, by copying the appropriate ni values.
*/
void stream(void)
{
    FILE *out = fopen("out.ppm", "wr");
    fprintf(out,"Streaming\n");
    for (int_t y = 0; y < M; y++) {
        for (int_t x = 0; x < N; x++) {
            fprintf(out,"%d ", sites[y*N+x]);
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

void options(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "ih")) != -1 ) {
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
            case ':':
                printf("Hello %s\n", optarg);
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
