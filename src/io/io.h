#ifndef RENDER_IMAGE_IO_H
#define RENDER_IMAGE_IO_H

/*
 * io.h — Snapshot I/O function prototypes.
 *
 * Supports Gadget-2 binary and HDF5 (Gadget-4/SWIFT-compatible) formats.
 * Include this header in files that read snapshots or manage sim_info:
 *   io.c, render_image.c, header.c
 *
 * To add a new format (e.g. SWIFT-native, GIZMO):
 *   1. Add reader functions here following the read_particles_from_*
 *      naming convention.
 *   2. Implement them in a new io_<format>.c file.
 *   3. Dispatch from render_image.c based on the -format flag.
 */

#include "types.h"   /* struct sim_info, ngb_buf_t */

/* ------------------------------------------------------------------ */
/* sim_info lifecycle                                                   */
/* ------------------------------------------------------------------ */

/* Allocate the dynamic arrays inside h for num_types particle species. */
void sim_info_alloc(struct sim_info *h, int num_types);

/* Free the dynamic arrays inside h (does not free h itself). */
void sim_info_free(struct sim_info *h);

/* ------------------------------------------------------------------ */
/* Neighbour buffer lifecycle                                           */
/* ------------------------------------------------------------------ */

/* Allocate a neighbour buffer sized for num_ngb expected neighbours.
 * cap = 4*num_ngb gives comfortable headroom for the sphere filter. */
ngb_buf_t *ngb_buf_alloc(int num_ngb);
void       ngb_buf_free(ngb_buf_t *buf);

/* ------------------------------------------------------------------ */
/* File detection                                                       */
/* ------------------------------------------------------------------ */

/* Probe filename, set SnapFormat, detect multi-file distributions. */
void check_input_filenames(char *filename, char *buffer,
                            int isHDF5, int *isDistributed);

/* ------------------------------------------------------------------ */
/* HDF5 reader                                                          */
/* ------------------------------------------------------------------ */

void read_hdf5_header(char *filename, struct sim_info *header,
                      long long *NumPart);

void read_particles_from_hdf5(char *filename,
                               float *x, float *y, float *z,
                               float *vx, float *vy, float *vz,
                               int *ptype, int ptype_mask,
                               long long *NumPartRead,
                               int read_velocities_flag);

/* Read particle velocities from HDF5 into vx/vy/vz.
 * Must be called after read_particles_from_hdf5 with the same file
 * and ptype_mask so array indices correspond. */
void read_velocities_from_hdf5(char *filename,
                                float *vx, float *vy, float *vz,
                                int NumFiles, long long *NThisTask);

/* ------------------------------------------------------------------ */
/* Gadget binary reader                                                 */
/* ------------------------------------------------------------------ */

void read_gadget_binary_header(char *filename, struct sim_info *header,
                                long long *NumPart);

void read_particles_from_gadget_binary(char *filename,
                                        float *x, float *y, float *z,
                                        int *ptype, int ptype_mask,
                                        long long *NumPartRead);

#endif /* RENDER_IMAGE_IO_H */
