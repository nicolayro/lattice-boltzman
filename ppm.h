#ifndef _PPM_H
#define _PPM_H
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

typedef int64_t int_t;

struct pixel_t {
    uint8_t r, g, b;
};

void image_from_ppm(const char *filename, struct pixel_t **image, int_t *height, int_t *width);
void image_print(FILE *out, struct pixel_t *image, int_t height, int_t width);
void image_print_ascii(FILE *out, struct pixel_t*image, int_t height, int_t width);

#endif
