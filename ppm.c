#include <stdlib.h>

#include "ppm.h"

void image_from_ppm(const char *filename, struct pixel_t **image, int_t *o_height, int_t *o_width)
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "ERROR: Unable to open file %s.\n", filename);
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
    *o_height = height;
    *o_width = width;

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

    *image = malloc(height * width * sizeof(struct pixel_t));
    if (!*image) {
        fprintf(stderr, "ERROR: Not enough memory...\n");
        exit(EXIT_FAILURE);
    }

    if ((c = fgetc(file)) != '\n') {
        ungetc(c, file);
    }

    int_t read = fread(*image, sizeof(struct pixel_t), height * width, file);
    if (read != height*width) {
        fprintf(stderr, "ERROR: Unable to read file: only %lld/%d objects were read.\n", read, height*width);
        exit(EXIT_FAILURE);
    }

    fclose(file);
}

void image_print(FILE *out, struct pixel_t *image, int_t height, int_t width)
{
    fprintf(out, "(%lld, %lld)\n", width, height);
    for (int_t y = 0; y < height; y++) {
        for (int_t x = 0; x < width; x++) {
            struct pixel_t pixel = image[y*width+x];
            fprintf(out, "\033[48;2;%d;%d;%dm  ", pixel.r, pixel.g, pixel.b);
        }
        fprintf(out, "\033[0m\n");
    }
}

char value_to_ascii(int_t b)
{
    assert(b <= 255);

    if (b >= 230)
        return ' ';
    else if (b >= 200)
        return '.';
    else if (b>=180)
        return '\'';
    else if (b >= 160)
        return ':';
    else if (b >= 130)
        return 'o';
    else if (b >= 100)
        return '&';
    else if (b >= 70)
        return '8';
    else if (b >= 50)
        return '#';
    else if (b >= 0)
        return '@';
    else
        assert(false && "negative brightness");
}

void image_print_ascii(FILE *out, struct pixel_t *image, int_t height, int_t width)
{
    for (int_t y = 0; y < height; y++) {
        for (int_t x = 0; x < width; x++) {
            struct pixel_t pixel = image[y*width+x];
            int_t brightness = ((int_t)pixel.r + (int_t)pixel.g + (int_t)pixel.b) / 3;
            fprintf(out, "%c", value_to_ascii(brightness));
        }
        fprintf(out, "\033[0m\n");
    }
}
