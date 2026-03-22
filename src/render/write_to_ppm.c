/*
 * write_to_ppm.c  –  image output routines
 *
 * write_to_png_ex() is the main entry point.  It applies the full
 * render_config_t pipeline:
 *
 *   raw density  →  density_normalise()  →  colormap_eval()
 *                                        →  opacity_eval()
 *                →  background composite  →  RGBA PNG
 *
 * write_to_png() is kept as a thin wrapper that uses default settings,
 * preserving backward compatibility with callers that don't need the
 * extra control.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <png.h>

#include "header.h"
#include "colormap.h"

/* ------------------------------------------------------------------ */
/* pixel_at / save helpers (unchanged)                                  */
/* ------------------------------------------------------------------ */

pixel_t *pixel_at(bitmap_t *bitmap, int x, int y)
{
    return bitmap->pixels + bitmap->width * y + x;
}

/*
 * save_png_to_file_rgba: write an RGBA bitmap to a PNG file.
 * The pixel stride is 4 bytes (R, G, B, A).
 */
static int save_png_to_file_rgba(uint8_t *rgba, int width, int height,
                                  const char *path)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;

    png_structp png_ptr  = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                    NULL, NULL, NULL);
    if (!png_ptr) { fclose(fp); return -1; }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(fp); return -1;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp); return -1;
    }

    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr,
                 (png_uint_32)width, (png_uint_32)height,
                 8,                        /* bit depth          */
                 PNG_COLOR_TYPE_RGBA,      /* RGBA               */
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);

    png_byte **rows = (png_byte **)png_malloc(
                          png_ptr, height * sizeof(png_byte *));
    for (int y = 0; y < height; y++)
        rows[y] = rgba + y * width * 4;

    png_set_rows(png_ptr, info_ptr, rows);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    png_free(png_ptr, rows);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Legacy RGB save (kept for write_to_ppm / old callers)               */
/* ------------------------------------------------------------------ */

int save_png_to_file(bitmap_t *bitmap, const char *path)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                   NULL, NULL, NULL);
    if (!png_ptr) { fclose(fp); return -1; }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(fp); return -1;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp); return -1;
    }

    png_set_IHDR(png_ptr, info_ptr,
                 (png_uint_32)bitmap->width, (png_uint_32)bitmap->height,
                 8, PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);

    png_byte **rows = (png_byte **)png_malloc(
                          png_ptr, bitmap->height * sizeof(png_byte *));
    for (size_t y = 0; y < bitmap->height; y++) {
        png_byte *row = (png_byte *)png_malloc(
                            png_ptr, sizeof(uint8_t) * bitmap->width * 3);
        rows[y] = row;
        for (size_t x = 0; x < bitmap->width; x++) {
            pixel_t *p = pixel_at(bitmap, (int)x, (int)y);
            *row++ = p->red;
            *row++ = p->green;
            *row++ = p->blue;
        }
    }

    png_init_io(png_ptr, fp);
    png_set_rows(png_ptr, info_ptr, rows);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    for (size_t y = 0; y < bitmap->height; y++)
        png_free(png_ptr, rows[y]);
    png_free(png_ptr, rows);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    return 0;
}

/* ------------------------------------------------------------------ */
/* write_to_ppm: legacy stub (kept for compatibility)                   */
/* ------------------------------------------------------------------ */

void write_to_ppm(char *image_file, int width, int height,
                   int MaxColorComponentValue, float *data)
{
    (void)data;
    FILE *fp = fopen(image_file, "wb");
    if (!fp) return;
    fprintf(fp, "P5\t %d \t %d \t %d\n",
            width, height, MaxColorComponentValue);
    fclose(fp);
}

/* ------------------------------------------------------------------ */
/* write_to_png_ex: full-featured RGBA output                           */
/*                                                                      */
/* data[]  – normalised density values, assumed to be in the raw       */
/*            (not yet log/linear mapped) units that match cfg->vmin   */
/*            and cfg->vmax.  If render_image.c has already normalised  */
/*            to [0,1] set cfg->log_scale=0, cfg->vmin=0, cfg->vmax=1. */
/* ------------------------------------------------------------------ */

void write_to_png_ex(const char *image_file, int width, int height,
                      const float *data, const render_config_t *cfg)
{
    int i, j;

    uint8_t *rgba = (uint8_t *)malloc((size_t)width * height * 4);
    if (!rgba) {
        fprintf(stderr, "write_to_png_ex: out of memory\n");
        return;
    }

    /*
     * Auto-scale vmin/vmax from the data if requested.
     * We need a mutable copy of cfg for this; use a local copy.
     */
    render_config_t lcfg = *cfg;
    if (lcfg.auto_levels)
        auto_levels_from_data(&lcfg, data, width * height);
    cfg = &lcfg;

    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            float raw = *((data + j * width) + i);

            /* Map raw density to [0,1] */
            float t = density_normalise(cfg, raw);

            /* Colour */
            float r, g, b;
            colormap_eval(cfg, t, &r, &g, &b);

            /* Opacity */
            float a = opacity_eval(cfg, t);

            /*
             * Composite over background colour using standard
             * "over" operator:
             *   out = src * src_a + bg * bg_a * (1 - src_a)
             * If bg_a == 1 and global_alpha == 1 this simplifies to
             * the usual opaque blend, and the output alpha is always 1.
             */
            float out_a = a + cfg->bg_a * (1.0f - a);
            float out_r, out_g, out_b;
            if (out_a > 0.0f) {
                out_r = (r * a + cfg->bg_r * cfg->bg_a * (1.0f - a)) / out_a;
                out_g = (g * a + cfg->bg_g * cfg->bg_a * (1.0f - a)) / out_a;
                out_b = (b * a + cfg->bg_b * cfg->bg_a * (1.0f - a)) / out_a;
            } else {
                out_r = out_g = out_b = 0.0f;
            }

            /* Clamp and convert to 8-bit */
            uint8_t *px = rgba + (j * width + i) * 4;
            px[0] = (uint8_t)(out_r > 1.0f ? 255 : (out_r < 0.0f ? 0 : (int)(out_r * 255.0f + 0.5f)));
            px[1] = (uint8_t)(out_g > 1.0f ? 255 : (out_g < 0.0f ? 0 : (int)(out_g * 255.0f + 0.5f)));
            px[2] = (uint8_t)(out_b > 1.0f ? 255 : (out_b < 0.0f ? 0 : (int)(out_b * 255.0f + 0.5f)));
            px[3] = (uint8_t)(out_a > 1.0f ? 255 : (out_a < 0.0f ? 0 : (int)(out_a * 255.0f + 0.5f)));
        }
    }

    if (save_png_to_file_rgba(rgba, width, height, image_file) != 0)
        fprintf(stderr, "write_to_png_ex: failed to write '%s'\n", image_file);

    free(rgba);
}

/* ------------------------------------------------------------------ */
/* write_to_png: backward-compatible wrapper using default config       */
/* ------------------------------------------------------------------ */

void write_to_png(char *image_file, int width, int height, float *data)
{
    render_config_t cfg;
    render_config_default(&cfg);
    /*
     * render_image.c already maps density to [0,1] before calling here,
     * so tell density_normalise to pass through linearly in [0,1].
     */
    cfg.log_scale = 0;
    cfg.vmin      = 0.0f;
    cfg.vmax      = 1.0f;
    write_to_png_ex(image_file, width, height, data, &cfg);
}
