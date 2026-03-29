#ifndef POSTPROCESS_H
#define POSTPROCESS_H

/*
 * postprocess.h — per-frame density map post-processing.
 *
 * Applied on the fully reduced global density map (after MPI_Reduce)
 * before colour mapping and PNG output.
 *
 * Steps (in order):
 *   1. Noise floor — zeros pixels below NOISE_FLOOR_FRACTION of peak,
 *      suppressing corner speckle from boundary particles with
 *      poorly-sampled kernels.
 *   2. CIC zero-hole inpainting — iterative 3×3 propagation fills
 *      empty pixels left by the sparse CIC stencil without touching
 *      or blurring any non-zero pixel.
 */

/* Fraction of peak density below which pixels are zeroed. */
#define NOISE_FLOOR_FRACTION 1.0e-4f

/*
 * Maximum number of inpainting passes.  Each pass expands the filled
 * frontier by one pixel, so this handles holes up to MAX_FILL_PASSES
 * pixels across.  64 comfortably covers worst-case CIC aliasing.
 */
#define MAX_FILL_PASSES 64

/*
 * postprocess_frame()
 *
 * Applies noise-floor zeroing and CIC hole inpainting in-place on
 * `data` (nx × ny floats, row-major).
 *
 * Must be called only on task 0, after MPI reduction.
 */
void postprocess_frame(float *data, int nx, int ny, float noise_floor_abs);

#endif /* POSTPROCESS_H */
