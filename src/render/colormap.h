#ifndef COLORMAP_H
#define COLORMAP_H

/*
 * colormap.h  –  colour palette and opacity control for render_image
 *
 * Everything that affects how the normalised density field [0,1] is
 * mapped to a final RGBA pixel lives in render_config_t.  Pass one of
 * these structs to write_to_png_ex() instead of the bare write_to_png()
 * call to get full control over:
 *
 *   • which colour palette is used
 *   • whether the palette is reversed
 *   • the transfer function used to turn density into opacity
 *   • a global opacity multiplier (useful for compositing)
 *   • the density normalisation range (vmin / vmax in log10 units)
 *   • whether to use linear or log10 density scaling
 *   • an optional background colour
 *
 * BUILT-IN PALETTES
 * -----------------
 * Pass the palette name string to colormap_by_name() to get the
 * corresponding enum value.
 *
 *   "coolwarm"    – blue → white → red  (diverging)
 *   "viridis"     – perceptually uniform purple→green→yellow
 *   "magma"       – perceptually uniform black→purple→white
 *   "inferno"     – perceptually uniform black→red→yellow→white
 *   "plasma"      – perceptually uniform blue→purple→yellow
 *   "hot"         – black → red → yellow → white
 *   "grayscale"   – black → white
 *   "fire"        – black → red → orange → yellow → white
 *   "ice"         – black → dark-blue → cyan → white
 *   "custom"      – loaded from a file (see colormap_load_file)
 *
 * OPACITY TRANSFER FUNCTIONS
 * --------------------------
 *   OPACITY_FLAT      – constant alpha = global_alpha
 *   OPACITY_LINEAR    – alpha ramps linearly from 0 at val=0 to
 *                       global_alpha at val=1
 *   OPACITY_SQRT      – alpha = sqrt(val) * global_alpha
 *   OPACITY_POWER     – alpha = val^opacity_gamma * global_alpha
 *   OPACITY_LOG       – alpha = log(1 + val*9)/log(10) * global_alpha
 *   OPACITY_THRESHOLD – alpha = 0 for val < opacity_threshold,
 *                       global_alpha otherwise
 */

/* ------------------------------------------------------------------ */
/* Enumerations                                                         */
/* ------------------------------------------------------------------ */

typedef enum {
    CMAP_COOLWARM = 0,
    CMAP_VIRIDIS,
    CMAP_MAGMA,
    CMAP_INFERNO,
    CMAP_PLASMA,
    CMAP_HOT,
    CMAP_GRAYSCALE,
    CMAP_FIRE,
    CMAP_ICE,
    CMAP_CUSTOM,
    CMAP_COUNT          /* sentinel */
} colormap_id_t;

typedef enum {
    OPACITY_FLAT = 0,
    OPACITY_LINEAR,
    OPACITY_SQRT,
    OPACITY_POWER,
    OPACITY_LOG,
    OPACITY_THRESHOLD
} opacity_func_t;

/* ------------------------------------------------------------------ */
/* Colour stop used for custom colourmap files                          */
/* ------------------------------------------------------------------ */
typedef struct {
    float pos;          /* position in [0,1]  */
    float r, g, b;      /* colour in [0,1]    */
} color_stop_t;

/* ------------------------------------------------------------------ */
/* Custom colourmap loaded from a file                                  */
/* ------------------------------------------------------------------ */
#define MAX_CUSTOM_STOPS 1024

typedef struct {
    int          n_stops;
    color_stop_t stops[MAX_CUSTOM_STOPS];
} custom_cmap_t;

/* ------------------------------------------------------------------ */
/* Main configuration struct                                            */
/* ------------------------------------------------------------------ */
typedef struct {
    /* --- Colour palette --- */
    colormap_id_t  colormap;          /* which built-in palette           */
    int            reverse;           /* 1 = reverse the palette          */

    /* --- Custom palette (used when colormap == CMAP_CUSTOM) --- */
    custom_cmap_t  custom;

    /* --- Opacity transfer function --- */
    opacity_func_t opacity_func;
    float          global_alpha;      /* overall opacity multiplier [0,1] */
    float          opacity_gamma;     /* exponent for OPACITY_POWER       */
    float          opacity_threshold; /* cutoff for OPACITY_THRESHOLD     */

    /* --- Density normalisation --- */
    float          vmin;              /* lower clip (log10 if log_scale)  */
    float          vmax;              /* upper clip                       */
    int            log_scale;         /* 1 = log10 density, 0 = linear    */

    /* --- Background --- */
    float          bg_r, bg_g, bg_b, bg_a; /* background RGBA [0,1]      */

    /* --- Auto-scaling --- */
    int            auto_levels;     /* 1 = derive vmin/vmax from data      */
    float          auto_pct_lo;     /* low  percentile to clip (default 0) */
    float          auto_pct_hi;     /* high percentile to clip (default 1) */
} render_config_t;

/* ------------------------------------------------------------------ */
/* API                                                                  */
/* ------------------------------------------------------------------ */

/* Fill *cfg with sensible defaults (viridis, flat alpha=1, log scale). */
void render_config_default(render_config_t *cfg);

/* Parse a single CLI-style key=value pair into *cfg.
   Returns 0 on success, -1 if the key is unrecognised.
   Usable in the argument-parsing loop of render_image.c. */
int render_config_parse_arg(render_config_t *cfg,
                             const char *key, const char *value);

/* Print a human-readable summary of the config to stdout. */
void render_config_print(const render_config_t *cfg);

/* Set cfg->vmin/vmax from the actual data range.
 * If auto_pct_lo/hi are set, clips the bottom/top percentile of
 * non-zero pixels before choosing the range.
 * Call this after smooth_to_mesh returns data[], before write_to_png_ex. */
void auto_levels_from_data(render_config_t *cfg,
                           const float *data, int npixels);

/* Look up a colormap_id_t by name string (case-insensitive).
   Returns CMAP_VIRIDIS on unknown names. */
colormap_id_t colormap_by_name(const char *name);

/*
 * Load a custom colourmap from a plain-text file.
 * Format: one entry per line, whitespace-separated:
 *     <pos>  <r>  <g>  <b>
 * All values in [0,1].  Lines starting with '#' are comments.
 * Returns 0 on success, -1 on error.
 */
int colormap_load_file(custom_cmap_t *cm, const char *path);

/* Evaluate the active colourmap at normalised value t in [0,1].
   Writes r,g,b in [0,1]. */
void colormap_eval(const render_config_t *cfg, float t,
                   float *r, float *g, float *b);

/* Evaluate the opacity transfer function at normalised value t in [0,1].
   Returns alpha in [0,1] (already multiplied by global_alpha). */
float opacity_eval(const render_config_t *cfg, float t);

/*
 * Normalise a raw density value to [0,1] according to cfg->vmin/vmax
 * and cfg->log_scale.  Values outside [vmin,vmax] are clamped.
 */
float density_normalise(const render_config_t *cfg, float raw);

#endif /* COLORMAP_H */
