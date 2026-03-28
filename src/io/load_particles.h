/*
 * load_particles.h — snapshot reading and position interpolation.
 *
 * Two interpolation modes are available:
 *
 *   load_particles_hermite()   — 2-snapshot Hermite, finite-difference tangents.
 *                                Tangent = (xB - xA), periodically wrapped.
 *                                Particles matched by array index (assumes
 *                                same ordering in both snapshots).
 *
 *   load_particles_hermite3()  — 3-snapshot Catmull-Rom Hermite.
 *                                Reads snap_prev, snap_a (middle), snap_b.
 *                                Particles matched by ParticleID (64-bit)
 *                                using radix sort — robust to reordering
 *                                between snapshots.
 *                                Tangents:
 *                                  m0 = wrap(xA - xPrev) / 2   (at left endpoint)
 *                                  m1 = wrap(xB - xA)    / 2   (at right endpoint,
 *                                                                one-sided backward diff)
 *                                With a 4th snapshot (snap_next) the right tangent
 *                                becomes the full centred difference:
 *                                  m1 = wrap(xNext - xA) / 2
 *
 * Both functions write periodically-wrapped interpolated positions into
 * the caller-supplied *x/*y/*z pointers (allocated internally).
 *
 * Mode selection in cli_args_t:
 *   snap_prev[0] != '\0'  →  hermite3_mode()  →  load_particles_hermite3()
 *   snap_a[0]   != '\0'   →  hermite_mode()   →  load_particles_hermite()
 *   otherwise             →  single-snapshot  →  load_particles()
 */

#ifndef LOAD_PARTICLES_H
#define LOAD_PARTICLES_H

#include "args.h"
#include "header.h"   /* struct sim_info */

/* ------------------------------------------------------------------ */
/* Header loading (shared by all paths)                                */
/* ------------------------------------------------------------------ */

int load_snapshot_header(cli_args_t      *cfg,
                         struct sim_info *header,
                         long long       *NumPart,
                         char            *filename_out,
                         int             *isDistributed_out);

/* ------------------------------------------------------------------ */
/* Single-snapshot path (legacy, linear extrapolation via velocity)    */
/* ------------------------------------------------------------------ */

int load_particles(const cli_args_t      *cfg,
                   const struct sim_info  *header,
                   const char             *filename,
                   float                **x,
                   float                **y,
                   float                **z,
                   int                  **ptype,
                   long long             *NumPartRead_out);

/* ------------------------------------------------------------------ */
/* Two-snapshot Hermite (index-matched, finite-difference tangents)    */
/* ------------------------------------------------------------------ */

int load_particles_hermite(const cli_args_t      *cfg,
                           const struct sim_info  *header,
                           float                **x,
                           float                **y,
                           float                **z,
                           int                  **ptype,
                           long long             *NumPartRead_out);

/* ------------------------------------------------------------------ */
/* Three-snapshot Catmull-Rom Hermite (ID-matched)                     */
/* ------------------------------------------------------------------ */

int load_particles_hermite3(const cli_args_t      *cfg,
                            const struct sim_info  *header,
                            float                **x,
                            float                **y,
                            float                **z,
                            int                  **ptype,
                            long long             *NumPartRead_out);

#endif /* LOAD_PARTICLES_H */
