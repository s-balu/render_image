#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <png.h>
#include <hdf5.h>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "header.h"
#include "colormap.h"
#include "args.h"
#include "config.h"
#include "load_particles.h"
#include "compute_smoothing.h"
#include "postprocess.h"

#ifndef IMAGE_DIMENSIONX
#define IMAGE_DIMENSIONX 768
#endif
#ifndef IMAGE_DIMENSIONY
#define IMAGE_DIMENSIONY 768
#endif

/* ------------------------------------------------------------------ */
/* Rodrigues rotation                                                   */
/*                                                                      */
/* Rotates the particle array (x,y,z) around `axis` through the point  */
/* (cx,cy,cz) by `frame * dangle` degrees, writing results to (rx,ry,  */
/* rz).  For frame 0 the angle is 0 so rx/ry/rz == x/y/z exactly.     */
/*                                                                      */
/* Formula: v_rot = v*cos(a) + (k x v)*sin(a) + k*(k.v)*(1-cos(a))    */
/* where k is the unit rotation axis and a is the angle in radians.    */
/* ------------------------------------------------------------------ */
static void rotate_particles(const float *x, const float *y, const float *z,
                             float *rx, float *ry, float *rz,
                             long long n,
                             double cx, double cy, double cz,
                             const float axis[3], float dangle_deg, int frame)
{
    float angle_rad = (float)frame * dangle_deg * (float)M_PI / 180.0f;
    float ca = cosf(angle_rad), sa = sinf(angle_rad);
    float kx = axis[0], ky = axis[1], kz = axis[2];

    for (long long pi = 0; pi < n; pi++) {
        float vx = x[pi] - (float)cx;
        float vy = y[pi] - (float)cy;
        float vz = z[pi] - (float)cz;
        float kdotv = kx*vx + ky*vy + kz*vz;
        float cpx = ky*vz - kz*vy;   /* cross product k x v */
        float cpy = kz*vx - kx*vz;
        float cpz = kx*vy - ky*vx;
        rx[pi] = (float)cx + vx*ca + cpx*sa + kx*kdotv*(1.0f - ca);
        ry[pi] = (float)cy + vy*ca + cpy*sa + ky*kdotv*(1.0f - ca);
        rz[pi] = (float)cz + vz*ca + cpz*sa + kz*kdotv*(1.0f - ca);
    }
}

/* ------------------------------------------------------------------ */
/* Usage                                                                */
/* ------------------------------------------------------------------ */
static void print_usage(const char *prog)
{
    fprintf(stdout,
        "Usage: %s [-config <file>] -input <file> -output <file> [options]\n"
        "\n"
        "  -config <file>          Load parameters from YAML file (CLI overrides)\n"
        "  -input  <file>          Snapshot root (no extension)\n"
        "  -output <file>          Output image root; frames written as root.NNNN.png\n"
        "  [-isHDF5]               HDF5 format (default: Gadget binary)\n"
        "  [-units 0|1]            0=simulation units, 1=box-normalised coords\n"
        "  [-xc <val>]             X centre of view volume\n"
        "  [-yc <val>]             Y centre of view volume\n"
        "  [-zc <val>]             Z centre of view volume\n"
        "  [-lbox <val>]           Side length of view volume\n"
        "  [-itmax <N>]            N frames, fixed box (default 1)\n"
        "  [-zoom <N>]             N frames, box shrinks each iter\n"
        "  [-zoom_factor <f>]      Zoom factor per iter, 0<f<1 (default 0.5)\n"
        "  [-rot_dangle <deg>]     Degrees to rotate per frame\n"
        "  [-rot_axis x,y,z]       Rotation axis (default 0,0,1 = z)\n"
        "  [-ptype <N>]            Keep particle type N (repeatable)\n"
        "  [-all_types]            Keep all particle types\n"
        "  [-gas]                  Shorthand: keep type 0\n"
        "  [-dark_matter]          Shorthand: keep type 1\n"
        "  [-stars]                Shorthand: keep type 4\n"
        "\n"
        "Colour / opacity:\n"
        "  [-scene cluster|scattered|filament]  Apply a named preset\n"
        "  [-colormap <name>]      viridis magma inferno plasma hot fire\n"
        "                          ice grayscale coolwarm custom\n"
        "  [-reverse_colormap]     Reverse the palette\n"
        "  [-opacity <val>]        Global opacity [0,1] (default 1)\n"
        "  [-opacity_func <name>]  flat linear sqrt power log threshold\n"
        "  [-opacity_gamma <val>]  Exponent for 'power' function\n"
        "  [-opacity_threshold <v>] Cutoff for 'threshold' function\n"
        "  [-vmin <val>]           Lower density clip (log10 units)\n"
        "  [-vmax <val>]           Upper density clip (log10 units)\n"
        "  [-linear_scale]         Use linear (not log10) density\n"
        "  [-colormap_file <path>] Load custom palette from file\n"
        "  [-bg_color R,G,B[,A]]  Background colour in [0,1]\n"
        "\n"
        "Density scaling (auto-levels on by default):\n"
        "  [-no_auto_levels]       Use fixed vmin/vmax\n"
        "  [-lock_levels]          Auto-level on frame 0, lock for all frames\n"
        "  [-auto_pct_lo <val>]    Low  clip percentile 0-1 (default 0.001)\n"
        "  [-auto_pct_hi <val>]    High clip percentile 0-1 (default 0.999)\n"
        "\n"
        "Time interpolation:\n"
        "  [-interp_frac <f>]      Shift positions by f*snap_dt*velocity\n"
        "  [-snap_dt <val>]        Snapshot interval in simulation time units\n"
        "\n"
        "Kernel smoothing (requires -DKERNEL_SMOOTHING):\n"
        "  [-num_ngb <N>]          Neighbours for h estimate (default 32)\n"
        "  [-sph_eta <val>]        h = eta*dist_Nth_ngb (default 1.2)\n"
        "  [-sph_cache <file>]     Cache smoothing lengths to file\n"
        "  [-fast_smooth]          O(N) grid h estimator vs exact KNN\n"
        "\n"
        "3-D grid:\n"
        "  [-ngrid_z <N>]          Depth of 3-D density grid (default min(width,256))\n",
        prog);
    fflush(stdout);
}

/* ------------------------------------------------------------------ */
/* main                                                                 */
/* ------------------------------------------------------------------ */
int main(int argc, char *argv[])
{
    /* ---- Configuration: defaults → YAML → CLI ---- */
    cli_args_t cfg;
    cli_args_default(&cfg);

#ifdef ENABLE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
    ThisTask = 0;
    NTask    = 1;
#endif

    if (argc < 2) {
        if (ThisTask == 0) print_usage(argv[0]);
#ifdef ENABLE_MPI
        MPI_Finalize();
#endif
        return 0;
    }

    /* Check to see if there is a configuration file being passed in... */
    for (int k = 1; k < argc - 1; k++) {
        if (strcmp(argv[k], "-config") == 0) {
            snprintf(cfg.config_path, sizeof(cfg.config_path),
                     "%s", argv[k + 1]);
            break;
        }
    }
    /* If yes, load in the parameters */
    if (cfg.config_path[0])
        load_yaml_config(cfg.config_path, &cfg);   /* YAML baseline */
    /* Parse command line arguments, if any have been passed in */
    parse_args(argc, argv, &cfg);                  /* CLI wins      */

    /* Load snapshot header */
    struct sim_info header;
    memset(&header, 0, sizeof(header));
    long long NumPart      = 0;
    long long NumPartRead  = 0;
    char      filename[256] = "";
    int       isDistributed = 0;

    if (load_snapshot_header(&cfg, &header, &NumPart,
                             filename, &isDistributed) != 0) {
#ifdef ENABLE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    /* Derive simulation-unit centre */
    double xc = cfg.xcen, yc = cfg.ycen, zc = cfg.zcen;

    /* Load particles */
    float *x = NULL, *y = NULL, *z = NULL;
    int   *ptype = NULL;

    if (hermite_mode(&cfg)) {
        /*
         * Two-snapshot Hermite interpolation.
         * snap_a and snap_b are set in the YAML file.
         * interp_frac is the parameter t in [0,1]:
         *   t=0 → exact snapshot A positions
         *   t=1 → exact snapshot B positions
         *   0<t<1 → smooth cubic interpolation between them
         *
         * With -itmax N the render loop generates N evenly-spaced
         * frames: t = 0, 1/(N-1), 2/(N-1), ..., 1.
         * With -interp_frac f a single frame at t=f is produced.
         */
        if (load_particles_hermite(&cfg, &header,
                                   &x, &y, &z, &ptype,
                                   &NumPartRead) != 0) {
#ifdef ENABLE_MPI
            MPI_Finalize();
#endif
            return 1;
        }
    } else {
        /* Single-snapshot path — behaviour unchanged */
        if (load_particles(&cfg, &header, filename,
                           &x, &y, &z, &ptype, &NumPartRead) != 0) {
#ifdef ENABLE_MPI
            MPI_Finalize();
#endif
            return 1;
        }
    }
    
//    if (load_particles(&cfg, &header, filename,
//                       &x, &y, &z, &ptype, &NumPartRead) != 0) {
//#ifdef ENABLE_MPI
//        MPI_Finalize();
//#endif
//        return 1;
//    }

    /* Select volume and distribute across MPI tasks */
    select_particles(x, y, z, ptype, header.BoxSize, NumPartRead,
                     xc, yc, zc, cfg.lbox, cfg.ptype_mask, &NumPart);

#ifdef ENABLE_MPI
    split_across_tasks_as_slabs(x, y, z, &NumPart, xc, yc, zc,
                                (float)cfg.lbox, (float)header.BoxSize);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (NumPart == 0) {
        if (ThisTask == 0)
            fprintf(stdout, "Error: no particles selected (mask=0x%x)\n",
                    (unsigned)cfg.ptype_mask);
#ifdef ENABLE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    if (ThisTask == 0) {
        fprintf(stdout, "Plotting %lld particles (type mask 0x%x)...\n",
                NumPart, (unsigned)cfg.ptype_mask);
        fflush(stdout);
    }

    /* Smoothing lengths (SPH mode) or NULL (CIC mode) */
    float *smoothing_length =
        compute_smoothing_lengths(&cfg, x, y, z, NumPart, xc, yc, zc);

    /* ---- Allocate image buffers ---- */
    float *data        = (float *)calloc(IMAGE_DIMENSIONX * IMAGE_DIMENSIONY,
                                         sizeof(float));
    float *global_data = (float *)calloc(IMAGE_DIMENSIONX * IMAGE_DIMENSIONY,
                                         sizeof(float));

    /* Working coordinate arrays for per-frame rotation */
    float *rx = (float *)malloc(NumPart * sizeof(float));
    float *ry = (float *)malloc(NumPart * sizeof(float));
    float *rz = (float *)malloc(NumPart * sizeof(float));
    if (!data || !global_data || !rx || !ry || !rz) {
        fprintf(stderr, "malloc failed for image/rotation buffers\n");
        return 1;
    }

    /* ---- Render loop ---- */
    /*
     * Modes:
     *   -itmax N              N frames, fixed box
     *   -zoom  N              N frames, box shrinks by zoom_factor each iter
     *   -rot_dangle D         rotate D degrees per frame around rot_axis
     */
    float BoundingBox = (float)cfg.lbox;
    int   n_frames    = (cfg.zoom > 0) ? cfg.zoom : cfg.itmax;
    char  image_file[256] = "";
    clock_t tstart, tfinish;

    for (int iter = 0; iter < n_frames; iter++) {
        memset(data,        0, IMAGE_DIMENSIONX * IMAGE_DIMENSIONY * sizeof(float));
        memset(global_data, 0, IMAGE_DIMENSIONX * IMAGE_DIMENSIONY * sizeof(float));
        
        /*
         * For Hermite sequences driven by -itmax N:
         * override interp_frac per frame inside the render loop.
         *
         * Add this at the TOP of the for (int iter = 0; ...) loop body,
         * before the rotate_particles() call:
         */
        
        /* Inside render loop, before rotate_particles(): */
        if (hermite_mode(&cfg) && cfg.itmax > 1) {
            /*
             * Recompute interpolated positions for this frame.
             * t advances uniformly from 0 at iter=0 to 1 at iter=itmax-1.
             */
            cfg.interp_frac = (cfg.itmax > 1)
            ? (float)iter / (float)(cfg.itmax - 1)
            : 0.0f;
            
            if (load_particles_hermite(&cfg, &header,
                                       &x, &y, &z, &ptype,
                                       &NumPartRead) != 0) {
                fprintf(stderr, "Hermite load failed at frame %d\n", iter);
                break;
            }
            
            /* Re-select view volume after reloading */
            long long NP = NumPartRead;
            select_particles(x, y, z, ptype, header.BoxSize, NP,
                             xc, yc, zc, cfg.lbox, cfg.ptype_mask, &NumPart);
        }

        rotate_particles(x, y, z, rx, ry, rz, NumPart,
                         xc, yc, zc, cfg.rot_axis, cfg.rot_dangle, iter);

        /* ngrid_z: cap at 256 to avoid enormous allocations at high res */
        int effective_ngrid_z = cfg.ngrid_z > 0
                              ? cfg.ngrid_z
                              : (IMAGE_DIMENSIONX < 256 ? IMAGE_DIMENSIONX : 256);

        tstart = clock();
        smooth_to_mesh(NumPart, smoothing_length, rx, ry, rz,
                       xc, yc, zc, BoundingBox,
                       (float)effective_ngrid_z,
                       IMAGE_DIMENSIONX, IMAGE_DIMENSIONY, data);
        tfinish = clock();
        fprintf(stdout, "smooth_to_mesh - time: %.2f s\n",
                (double)(tfinish - tstart) / CLOCKS_PER_SEC);
        fflush(stdout);

#ifdef ENABLE_MPI
        MPI_Reduce(data, global_data,
                   IMAGE_DIMENSIONX * IMAGE_DIMENSIONY,
                   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        memcpy(global_data, data,
               IMAGE_DIMENSIONX * IMAGE_DIMENSIONY * sizeof(float));
#endif

        if (ThisTask == 0) {
            int npix_g = IMAGE_DIMENSIONX * IMAGE_DIMENSIONY;

            /*
             * Lock levels before inpainting: auto-levels must be computed
             * from pre-inpaint data so inpainted void pixels don't drag
             * vmin down and corrupt the colour scale.
             */
            if (cfg.rcfg.auto_levels && cfg.lock_levels && iter == 0) {
                auto_levels_from_data(&cfg.rcfg, global_data, npix_g);
                fprintf(stdout,
                        "Levels locked (pre-inpaint): vmin=%.3f vmax=%.3f\n",
                        cfg.rcfg.vmin, cfg.rcfg.vmax);
                fflush(stdout);
                cfg.rcfg.auto_levels = 0;   /* freeze for subsequent frames */
            }

            /* After locking levels, compute absolute noise floor in
             * linear density units.
             * vmin is log10 density, so the linear noise floor is
             * 10^vmin * NOISE_FLOOR_FRACTION. Using vmin directly means
             * the threshold tracks the locked colour scale exactly.
             */
            float noise_floor_abs = powf(10.0f, cfg.rcfg.vmin) * NOISE_FLOOR_FRACTION;

            postprocess_frame(global_data, IMAGE_DIMENSIONX, IMAGE_DIMENSIONY, noise_floor_abs);

            snprintf(image_file, sizeof(image_file),
                     "%s.%04d.png", cfg.image_file_root, iter);
            write_to_png_ex(image_file,
                            IMAGE_DIMENSIONX, IMAGE_DIMENSIONY,
                            global_data, &cfg.rcfg);
        }

        if (cfg.zoom > 0) BoundingBox *= cfg.zoom_factor;
    }

    /* Free up memory */
    free(rx); free(ry); free(rz);
    free(data); free(global_data);
    free(x); free(y); free(z); free(ptype);
    if (smoothing_length) free(smoothing_length);
    sim_info_free(&header);

    fprintf(stdout, "Finished.\n");
    fflush(stdout);

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
    return 0;
}
