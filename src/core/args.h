#ifndef ARGS_H
#define ARGS_H

/*
 * args.h — CLI argument struct and parser for render_image.
 *
 * Typical call sequence in main():
 *
 *   cli_args_t args;
 *   cli_args_default(&args);
 *   load_yaml_config(args.config_path, &args);  // YAML baseline
 *   parse_args(argc, argv, &args);              // CLI wins
 *
 * Interpolation mode selection (checked in this priority order):
 *
 *   hermite3_mode()  snap_prev + snap_a + snap_b set
 *                    → load_particles_hermite3()   Catmull-Rom, ID-matched
 *                    snap_next is optional; if set, gives a fully centred
 *                    difference at the right endpoint (4-snapshot mode).
 *
 *   hermite_mode()   snap_a + snap_b set (snap_prev empty)
 *                    → load_particles_hermite()    2-snapshot, index-matched
 *
 *   (neither)        → load_particles()            single snapshot
 */

#include "colormap.h"   /* render_config_t */

typedef struct {

    /* ----- Config file ----- */
    char config_path[512];      /* -config <file> */

    /* ----- I/O ----- */
    char file_root[256];        /* -input  <root>  (single snapshot or snap_a) */
    char image_file_root[256];  /* -output <root> */
    int  isHDF5;                /* -isHDF5 */

    /* ----- Hermite interpolation snapshots ----- */
    /*
     * 2-snapshot mode (hermite_mode):
     *   snap_a  — left  endpoint (t=0)
     *   snap_b  — right endpoint (t=1)
     *
     * 3-snapshot Catmull-Rom mode (hermite3_mode):
     *   snap_prev — snapshot before snap_a  (tangent at left  endpoint)
     *   snap_a    — left  endpoint (t=0)
     *   snap_b    — right endpoint (t=1)
     *   snap_next — snapshot after  snap_b  (tangent at right endpoint, optional)
     *               if omitted, a one-sided backward difference is used at snap_b
     *
     * All four paths are YAML-only (paths are too long for convenient CLI use).
     * snap_prev and snap_next can also be set via -snap_prev / -snap_next CLI flags.
     */
    char snap_prev[512];        /* snap_prev: path — snapshot before snap_a  */
    char snap_a[512];           /* snap_a:    path — first  snapshot (t=0)   */
    char snap_b[512];           /* snap_b:    path — second snapshot (t=1)   */
    char snap_next[512];        /* snap_next: path — snapshot after  snap_b  */

    /* ----- View / simulation ----- */
    int    boxunits;            /* -units 0|1     */
    double xcen, ycen, zcen;   /* -xc -yc -zc    */
    double lbox;                /* -lbox          */
    int    itmax;               /* -itmax N       */
    int    zoom;                /* -zoom  N       */
    float  zoom_factor;         /* -zoom_factor f */

    /* ----- Rotation ----- */
    float rot_axis[3];          /* -rot_axis x,y,z */
    float rot_dangle;           /* -rot_dangle deg  */

    /* ----- Particle selection ----- */
    int ptype_mask;             /* -ptype / -gas / -dark_matter / -stars / -all_types */

    /* ----- Colour / opacity ----- */
    render_config_t rcfg;       /* -colormap, -opacity_*, -vmin/vmax, -scene, … */
    int lock_levels;            /* -lock_levels */

    /* ----- Kernel smoothing ----- */
    int   num_ngb;              /* -num_ngb N      */
    float sph_eta;              /* -sph_eta val    */
    char  sph_cache[512];       /* -sph_cache file */
    int   fast_smooth;          /* -fast_smooth    */

    /* ----- Time interpolation ----- */
    float interp_frac;          /* -interp_frac f  — parameter t in [0,1]          */
    float snap_dt;              /* -snap_dt val    — legacy; ignored in Hermite3    */

    /* ----- 3-D grid ----- */
    int ngrid_z;                /* -ngrid_z N */

} cli_args_t;

/* ------------------------------------------------------------------ */
/* Function declarations                                                */
/* ------------------------------------------------------------------ */

void cli_args_default(cli_args_t *args);
int  parse_args(int argc, char *argv[], cli_args_t *args);
void apply_scene_preset(cli_args_t *args, const char *scene);

/* ------------------------------------------------------------------ */
/* Mode predicates                                                      */
/* ------------------------------------------------------------------ */

/*
 * hermite3_mode() — Catmull-Rom 3/4-snapshot interpolation.
 * Requires snap_prev + snap_a + snap_b.  snap_next is optional.
 * Takes priority over hermite_mode().
 */
static inline int hermite3_mode(const cli_args_t *args)
{
    return (args->snap_prev[0] != '\0' &&
            args->snap_a[0]   != '\0' &&
            args->snap_b[0]   != '\0');
}

/*
 * hermite_mode() — 2-snapshot finite-difference Hermite interpolation.
 * Requires snap_a + snap_b.  Only active when snap_prev is not set.
 */
static inline int hermite_mode(const cli_args_t *args)
{
    return (args->snap_prev[0] == '\0' &&
            args->snap_a[0]   != '\0' &&
            args->snap_b[0]   != '\0');
}

#endif /* ARGS_H */
