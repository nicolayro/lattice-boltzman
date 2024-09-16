#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

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

void options(int argc, char **argv);

int_t N, M;
int_t timesteps;
char *input = NULL;

domain_t *domain = NULL;

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

    for (int_t i = 0; i < timesteps; i++) {
        collide();
        stream();
    }

    free(domain);

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

    domain = malloc(M * N * sizeof(uint8_t));
    for (int_t y = 0; y < M; y++) {
        for (int_t x = 0; x < N; x++) {
            int16_t value = geometry[3*(y*width+x)] + geometry[3*(y*width+x)+1]
                + geometry[3*(y*width+x)+2];

            domain[y*width+x] = value > 0 ? FLUID : SOLID;
        }
    }

    // All SOLID points that are next to FLUID points are categorized as WALL
    for (int_t y = 0; y < M; y++ ) {
        for (int_t x = 0; x < N; x++) {
            if (domain[y*width+x] != SOLID)
                continue;

            for (int_t i = 0; i < 6; i++) {
                int_t nx = x + OFFSETS[y%2][i][0];
                int_t ny = x + OFFSETS[y%2][i][1];

                bool in_bounds = (nx >= 0 || nx < N || ny >= 0 || ny < M);
                if (in_bounds && domain[ny*width+nx] == FLUID)
                    domain[y*width+x] = WALL;
            }
        }
    }

    free(geometry);
    fclose(file);
}

void collide(void)
{
    FILE *out = fopen("out.ppm", "wr");
    fprintf(out,"Colliding\n");
    for (int_t y = 0; y < M; y++) {
        for (int_t x = 0; x < N; x++) {
            fprintf(out,"%d ", domain[y*N+x]);
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

void stream(void)
{
    FILE *out = fopen("out.ppm", "wr");
    fprintf(out,"Streaming\n");
    for (int_t y = 0; y < M; y++) {
        for (int_t x = 0; x < N; x++) {
            fprintf(out,"%d ", domain[y*N+x]);
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
