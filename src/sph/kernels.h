#ifndef RENDER_IMAGE_KERNELS_H
#define RENDER_IMAGE_KERNELS_H

/*
 * kernels.h — SPH kernel and projection function prototypes.
 *
 * The cubic spline kernel W(r,h) with compact support at r=2h is used
 * throughout.  Two evaluation paths are provided:
 *
 *   get_kernel_value(u)      — 3D kernel at u = r/h
 *   get_kernel_value_u2(u2)  — 3D kernel at u^2 = (r/h)^2, no sqrtf
 *   get_projected_kernel_value(q) — analytically projected 2D kernel at
 *                                   projected radius q = R/h
 *
 * All hot evaluation paths use lookup tables initialised by the
 * tabulate_*() functions.  Call the appropriate tabulate function once
 * at startup before any kernel evaluations.
 */

/* ------------------------------------------------------------------ */
/* 3D kernel                                                            */
/* ------------------------------------------------------------------ */

/* Initialise the 3D kernel lookup table.  Call once before get_kernel_value(). */
void tabulate_kernel(void);

#ifdef KERNEL_SMOOTHING
/* u2-domain table: avoids sqrtf in the inner deposit loop. */
void tabulate_kernel_u2(void);

/* Evaluate W(u) at u = r/h via table lookup + linear interpolation. */
float get_kernel_value(float u);

/* Evaluate W(sqrt(u2)) at u2 = (r/h)^2 — no sqrtf needed. */
float get_kernel_value_u2(float u2);

/* Raw kernel function (slow — used only to populate tables). */
float cubic_spline_kernel(float u);
#endif /* KERNEL_SMOOTHING */

/* ------------------------------------------------------------------ */
/* 2D projected kernel                                                  */
/* ------------------------------------------------------------------ */

/* Initialise the projected kernel table.  Call once before
 * get_projected_kernel_value(). */
void tabulate_projected_kernel(void);

/* Initialise the q^2-domain projected kernel table.
 * Used by deposit_sph_2d to avoid sqrtf in the inner pixel loop. */
void tabulate_projected_kernel_u2(void);

/* Evaluate the 2D projected kernel at q = R/h. */
float get_projected_kernel_value(float q);

/* Evaluate the 2D projected kernel at q2 = (R/h)^2 — no sqrtf needed. */
float get_projected_kernel_value_u2(float q2);

/* Raw 2D projected kernel (slow — used only to populate tables). */
float cubic_spline_kernel_2d_proj(float q);

/* ------------------------------------------------------------------ */
/* Cumulative kernel integral (used by smooth_to_mesh CIC path)         */
/* ------------------------------------------------------------------ */
void  tabulate_integral(void);
float get_integral_value(float u);
float cumulative_cubic_spline_interpolant(float u);

#endif /* RENDER_IMAGE_KERNELS_H */
