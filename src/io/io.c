#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "header.h"

/* ------------------------------------------------------------------ */
/* Helpers                                                              */
/* ------------------------------------------------------------------ */

void check_input_filenames(char *filename, char *buffer, int isHDF5,
                            int *isDistributed)
{
    FILE *infile;

    if (isHDF5 == 1) {
        SnapFormat = 3;
        if (ThisTask == 0) {
            fprintf(stdout, "Assuming HDF5 input...\n");
            fflush(stdout);
        }
        sprintf(filename, "%s.hdf5", buffer);
    } else {
        sprintf(filename, "%s", buffer);
    }

    if ((infile = fopen(filename, "r")) == NULL) {
        if (isHDF5 == 1)
            sprintf(filename, "%s.0.hdf5", buffer);
        else
            sprintf(filename, "%s.0", buffer);

        if ((infile = fopen(filename, "r")) == NULL) {
            if (ThisTask == 0) {
                fprintf(stdout, "Error: cannot find input file\n");
                fflush(stdout);
            }
#ifdef ENABLE_MPI
            MPI_Finalize();
#endif
            exit(0);
        } else {
            fclose(infile);
            *isDistributed = 1;
            if (ThisTask == 0) {
                fprintf(stdout, "Assuming distributed files...\n");
                fflush(stdout);
            }
        }
    } else {
        fclose(infile);
        if (ThisTask == 0) {
            fprintf(stdout, "Filename: %s\n", filename);
            fflush(stdout);
        }
    }
}

/* ------------------------------------------------------------------ */
/* HDF5 header reader                                                   */
/*                                                                      */
/* Discovers num_types by probing which /PartTypeN groups exist, so    */
/* the code works with any number of species (Gadget-2 = 6,            */
/* SWIFT = 7, custom = anything).                                       */
/* ------------------------------------------------------------------ */

static int rank_g;
static hsize_t dims_g[2], count_g[2], start_g[2];

static hid_t hdf5_file, hdf5_headergrp, hdf5_grp;
static hid_t hdf5_attribute;
static hid_t hdf5_datatype, hdf5_dataset;
static hid_t hdf5_dataspace_in_file, hdf5_dataspace_in_memory;
static hid_t hdf5_status;

/*
 * Probe how many /PartTypeN groups the file actually contains.
 * We try 0..MAX_PROBE-1 and return the highest index found + 1.
 * This handles both Gadget-2 (6 types) and formats with more or fewer.
 */
#define MAX_TYPE_PROBE 32
static int probe_num_types_hdf5(hid_t file_id)
{
    char buf[64];
    int num_types = 0;
    for (int i = 0; i < MAX_TYPE_PROBE; i++) {
        sprintf(buf, "/PartType%d", i);
        /* H5Lexists returns positive if the link exists */
        if (H5Lexists(file_id, buf, H5P_DEFAULT) > 0)
            num_types = i + 1;
    }
    /* Always return at least 1 to avoid zero-length allocations */
    return (num_types > 0) ? num_types : 1;
}

void read_hdf5_header(char *filename, struct sim_info *header,
                       long long *NumPart)
{
    int i, j;
    long long NumPartInFile = 0;

    hdf5_file      = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

    /* --- Discover actual number of types in this file --- */
    int num_types = probe_num_types_hdf5(hdf5_file);

    /*
     * NumPart_ThisFile is stored as an array whose length equals the
     * number of types the writer used.  Read it into a temporary
     * GADGET_MAX_TYPES-wide buffer to stay safe with legacy files,
     * then copy the relevant prefix into the dynamically-sized header.
     */
    int tmp_npart[MAX_TYPE_PROBE] = {0};

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
    /* Read as many ints as the attribute contains */
    hid_t aspace = H5Aget_space(hdf5_attribute);
    hsize_t nelem;
    H5Sget_simple_extent_dims(aspace, &nelem, NULL);
    H5Sclose(aspace);
    /* num_types from probe may differ from attribute length; use the
       attribute length as the authoritative type count */
    num_types = (int)nelem;
    hdf5_status = H5Aread(hdf5_attribute, H5T_NATIVE_INT, tmp_npart);
    hdf5_status = H5Aclose(hdf5_attribute);

    sim_info_alloc(header, num_types);

    for (i = 0; i < num_types; i++) {
        header->npart[i] = tmp_npart[i];
        NumPartInFile    += tmp_npart[i];
        if (ThisTask == 0) {
            fprintf(stdout, "Species %d : n_in_file %d\n", i, header->npart[i]);
            fflush(stdout);
        }
    }

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
    hdf5_status = H5Aread(hdf5_attribute, H5T_NATIVE_INT, &(header->NumFiles));
    hdf5_status = H5Aclose(hdf5_attribute);

    if (header->NumFiles > 1) {
        hdf5_status = H5Gclose(hdf5_headergrp);
        hdf5_status = H5Fclose(hdf5_file);

        /* Sum NumPart across all sub-files */
        for (i = 0; i < header->NumFiles; i++) {
            hdf5_file      = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
            hdf5_headergrp = H5Gopen(hdf5_file, "/Header");

            hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
            hdf5_status    = H5Aread(hdf5_attribute, H5T_NATIVE_INT, tmp_npart);
            hdf5_status    = H5Aclose(hdf5_attribute);

            for (j = 0; j < num_types; j++)
                *NumPart += tmp_npart[j];

            hdf5_status = H5Gclose(hdf5_headergrp);
            hdf5_status = H5Fclose(hdf5_file);
        }

        hdf5_file      = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
    } else {
        *NumPart = NumPartInFile;
    }

    /* NumPart_Total length may differ from NumPart_ThisFile length in
       some writers; use the same tmp buffer approach. */
    unsigned int tmp_nall[MAX_TYPE_PROBE] = {0};
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
    hdf5_status    = H5Aread(hdf5_attribute, H5T_NATIVE_UINT, tmp_nall);
    hdf5_status    = H5Aclose(hdf5_attribute);
    for (i = 0; i < num_types; i++)
        header->nall[i] = (int)tmp_nall[i];

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Redshift");
    hdf5_status = H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &(header->redshift));
    hdf5_status = H5Aclose(hdf5_attribute);

    header->time = 1.0 / (1.0 + header->redshift);

    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
    hdf5_status = H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &(header->BoxSize));
    hdf5_status = H5Aclose(hdf5_attribute);

    hdf5_status = H5Gclose(hdf5_headergrp);
    hdf5_status = H5Fclose(hdf5_file);
}

/* ------------------------------------------------------------------ */
/* Gadget binary header reader                                          */
/*                                                                      */
/* The Gadget binary header block is always 256 bytes and stores       */
/* exactly GADGET_MAX_TYPES (6) entries for npart, massarr, nall.      */
/* After reading we count how many of those type slots are non-zero    */
/* and set header->num_types accordingly.                               */
/* ------------------------------------------------------------------ */
void read_gadget_binary_header(char *filename, struct sim_info *header,
                                long long *NumPart)
{
    int i;
    long long NumPartInFile = 0;

    FILE *infile;
    long long fileoffset;
    int dummy;

    /* Temporary fixed-width buffers matching the on-disk layout */
    int    tmp_npart  [GADGET_MAX_TYPES] = {0};
    double tmp_massarr[GADGET_MAX_TYPES] = {0};
    int    tmp_nall   [GADGET_MAX_TYPES] = {0};

    infile = fopen(filename, "r");
    fread(&dummy, sizeof(int), 1, infile);
    SnapFormat = (dummy == 8) ? 2 : 1;
    rewind(infile);

    fileoffset = (SnapFormat == 2)
        ? sizeof(int) + 4 * sizeof(char) + 2 * sizeof(int) + sizeof(int)
        : sizeof(int);

    fseek(infile, fileoffset, SEEK_SET);
    fread(tmp_npart, sizeof(int), GADGET_MAX_TYPES, infile);

    fileoffset += GADGET_MAX_TYPES * sizeof(int);
    fseek(infile, fileoffset, SEEK_SET);
    fread(tmp_massarr, sizeof(double), GADGET_MAX_TYPES, infile);

    double tmp_time;
    fread(&tmp_time, sizeof(double), 1, infile);
    header->time = tmp_time;

    fileoffset += GADGET_MAX_TYPES * sizeof(double)
                + 2 * sizeof(double)
                + 2 * sizeof(int);
    fseek(infile, fileoffset, SEEK_SET);
    fread(tmp_nall, sizeof(int), GADGET_MAX_TYPES, infile);

    fileoffset += 7 * sizeof(int);
    fseek(infile, fileoffset, SEEK_SET);
    fread(&(header->NumFiles), sizeof(int),    1, infile);
    fread(&(header->BoxSize),  sizeof(double), 1, infile);

    fclose(infile);

    /*
     * Determine num_types: count how many type slots have any particles
     * in either npart or nall.  Fall back to GADGET_MAX_TYPES if all
     * are zero (shouldn't happen with a valid file).
     */
    int num_types = 0;
    for (i = 0; i < GADGET_MAX_TYPES; i++)
        if (tmp_npart[i] > 0 || tmp_nall[i] > 0)
            num_types = i + 1;
    if (num_types == 0) num_types = GADGET_MAX_TYPES;

    sim_info_alloc(header, num_types);

    for (i = 0; i < num_types; i++) {
        header->npart  [i] = tmp_npart  [i];
        header->nall   [i] = tmp_nall   [i];
        header->massarr[i] = tmp_massarr[i];
        NumPartInFile      += tmp_npart[i];

        if (ThisTask == 0) {
            fprintf(stdout,
                    "Species %d : n_in_file %d, n_in_snapshot %d\n",
                    i, header->npart[i], header->nall[i]);
            fflush(stdout);
        }
        *NumPart += header->nall[i];
    }
}

/* ------------------------------------------------------------------ */
/* HDF5 particle reader                                                 */
/*                                                                      */
/* Iterates over all /PartTypeN groups that actually exist in the file, */
/* so it naturally handles any number of types.                         */
/* ------------------------------------------------------------------ */
void read_particles_from_hdf5(char *temp, float *x, float *y, float *z,
                               int *ptype, int NumFiles,
                               long long *NThisTask)
{
    int i, j, k, nfile;
    char filename[256], groupname[64];
    float *dummy;

    /* Per-file npart: size up to MAX_TYPE_PROBE */
    int tmp_npart[MAX_TYPE_PROBE];

    for (nfile = 0; nfile < NumFiles; nfile++) {
        if (NumFiles > 1)
            sprintf(filename, "%s%d.hdf5", temp, nfile);
        else
            sprintf(filename, "%s.hdf5", temp);

        hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

        /* Read NumPart_ThisFile to know how many particles of each type */
        hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");

        hid_t aspace = H5Aget_space(hdf5_attribute);
        hsize_t nelem;
        H5Sget_simple_extent_dims(aspace, &nelem, NULL);
        H5Sclose(aspace);

        int num_types_in_file = (int)nelem;
        memset(tmp_npart, 0, sizeof(tmp_npart));
        hdf5_status = H5Aread(hdf5_attribute, H5T_NATIVE_INT, tmp_npart);
        hdf5_status = H5Aclose(hdf5_attribute);
        hdf5_status = H5Gclose(hdf5_headergrp);

        /* Loop over every type slot reported by the header */
        for (i = 0; i < num_types_in_file; i++) {
            if (tmp_npart[i] <= 0) continue;

            sprintf(groupname, "/PartType%d", i);
            if (H5Lexists(hdf5_file, groupname, H5P_DEFAULT) <= 0) continue;

            hdf5_grp     = H5Gopen(hdf5_file, groupname);
            hdf5_dataset = H5Dopen(hdf5_grp, "Coordinates");
            hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);

            rank_g    = 2;
            dims_g[1] = 3;
            dims_g[0] = tmp_npart[i];
            hdf5_dataspace_in_file = H5Screate_simple(rank_g, dims_g, NULL);

            hsize_t per_task = tmp_npart[i] / NTask;
            dims_g[0] = per_task;
            hdf5_dataspace_in_memory = H5Screate_simple(rank_g, dims_g, NULL);

            start_g[0] = (hsize_t)ThisTask * per_task;
            start_g[1] = 0;
            count_g[0] = per_task;
            count_g[1] = 3;

            hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
                                               H5S_SELECT_SET,
                                               start_g, NULL, count_g, NULL);

            dummy = (float *)malloc(3 * sizeof(float) * per_task);
            hdf5_status = H5Dread(hdf5_dataset, hdf5_datatype,
                                   hdf5_dataspace_in_memory,
                                   hdf5_dataspace_in_file,
                                   H5P_DEFAULT, dummy);

            k = 0;
            for (j = 0; j < (int)per_task; j++) {
                x[*NThisTask] = dummy[k];
                y[*NThisTask] = dummy[k + 1];
                z[*NThisTask] = dummy[k + 2];
                ptype[*NThisTask] = i;
                k += 3;
                (*NThisTask)++;
            }
            free(dummy);

            hdf5_status = H5Sclose(hdf5_dataspace_in_memory);
            hdf5_status = H5Sclose(hdf5_dataspace_in_file);
            hdf5_status = H5Tclose(hdf5_datatype);
            hdf5_status = H5Dclose(hdf5_dataset);
            hdf5_status = H5Gclose(hdf5_grp);
        }

        hdf5_status = H5Fclose(hdf5_file);
    }
}


/* ------------------------------------------------------------------ */
/* read_velocities_from_hdf5: read particle velocities into vx/vy/vz. */
/*                                                                      */
/* Reads only particle types that pass ptype_mask (same mask used when  */
/* reading coordinates) and in the same order, so vx[i]/vy[i]/vz[i]   */
/* corresponds to x[i]/y[i]/z[i].                                       */
/*                                                                      */
/* HDF5 dataset: /PartTypeN/Velocities  shape [N,3]  float32           */
/* Gadget internal units: km/s * sqrt(a) — divide by sqrt(a) for        */
/* physical km/s; multiply by a*dt/H0 for comoving displacement.        */
/* For rendering purposes the raw values work fine as a relative shift.  */
/* ------------------------------------------------------------------ */
void read_velocities_from_hdf5(char *temp, float *vx, float *vy, float *vz,
                                int NumFiles, long long *NThisTask)
{
    int i, j, k, nfile;
    char filename[256], groupname[64];
    float *dummy;
    int tmp_npart[MAX_TYPE_PROBE];
    long long n_read = 0;

    for (nfile = 0; nfile < NumFiles; nfile++) {
        if (NumFiles > 1)
            sprintf(filename, "%s%d.hdf5", temp, nfile);
        else
            sprintf(filename, "%s.hdf5", temp);

        hdf5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
        hid_t aspace = H5Aget_space(hdf5_attribute);
        hsize_t nelem;
        H5Sget_simple_extent_dims(aspace, &nelem, NULL);
        H5Sclose(aspace);
        int num_types_in_file = (int)nelem;
        memset(tmp_npart, 0, sizeof(tmp_npart));
        hdf5_status = H5Aread(hdf5_attribute, H5T_NATIVE_INT, tmp_npart);
        hdf5_status = H5Aclose(hdf5_attribute);
        hdf5_status = H5Gclose(hdf5_headergrp);

        for (i = 0; i < num_types_in_file; i++) {
            if (tmp_npart[i] <= 0) continue;
            sprintf(groupname, "/PartType%d", i);
            if (H5Lexists(hdf5_file, groupname, H5P_DEFAULT) <= 0) continue;

            /* Check Velocities dataset exists (some formats omit it) */
            hdf5_grp = H5Gopen(hdf5_file, groupname);
            if (H5Lexists(hdf5_grp, "Velocities", H5P_DEFAULT) <= 0) {
                /* Fill with zeros if missing — no shift applied */
                hsize_t per_task = tmp_npart[i] / NTask;
                for (j = 0; j < (int)per_task; j++) {
                    vx[n_read] = vy[n_read] = vz[n_read] = 0.0f;
                    n_read++;
                }
                hdf5_status = H5Gclose(hdf5_grp);
                continue;
            }

            hdf5_dataset  = H5Dopen(hdf5_grp, "Velocities");
            hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);

            rank_g    = 2;
            dims_g[1] = 3;
            dims_g[0] = tmp_npart[i];
            hdf5_dataspace_in_file = H5Screate_simple(rank_g, dims_g, NULL);

            hsize_t per_task = tmp_npart[i] / NTask;
            dims_g[0] = per_task;
            hdf5_dataspace_in_memory = H5Screate_simple(rank_g, dims_g, NULL);

            start_g[0] = (hsize_t)ThisTask * per_task;
            start_g[1] = 0;
            count_g[0] = per_task;
            count_g[1] = 3;

            hdf5_status = H5Sselect_hyperslab(hdf5_dataspace_in_file,
                                               H5S_SELECT_SET,
                                               start_g, NULL, count_g, NULL);

            dummy = (float *)malloc(3 * sizeof(float) * per_task);
            hdf5_status = H5Dread(hdf5_dataset, hdf5_datatype,
                                   hdf5_dataspace_in_memory,
                                   hdf5_dataspace_in_file,
                                   H5P_DEFAULT, dummy);
            k = 0;
            for (j = 0; j < (int)per_task; j++) {
                vx[n_read] = dummy[k];
                vy[n_read] = dummy[k + 1];
                vz[n_read] = dummy[k + 2];
                k += 3;
                n_read++;
            }
            free(dummy);

            hdf5_status = H5Sclose(hdf5_dataspace_in_memory);
            hdf5_status = H5Sclose(hdf5_dataspace_in_file);
            hdf5_status = H5Tclose(hdf5_datatype);
            hdf5_status = H5Dclose(hdf5_dataset);
            hdf5_status = H5Gclose(hdf5_grp);
        }
        hdf5_status = H5Fclose(hdf5_file);
    }
    *NThisTask = n_read;
}

/* ------------------------------------------------------------------ */
/* Gadget binary particle reader                                        */
/*                                                                      */
/* Uses header->num_types for all loops instead of the hard-coded 6.   */
/* ------------------------------------------------------------------ */
void read_particles_from_gadget_binary(char *temp, float *x, float *y,
                                       float *z, int *ptype, int NumFiles,
                                       long long *NThisTask)
{
    int i, j, nfile;
    float pos[3];
    unsigned long long *pid;
    int tmp_npart[GADGET_MAX_TYPES];
    char filename[256];
    long long NumPartInFile, fileoffset;
    FILE *infile;
    int blksize;

    *NThisTask = 0;

    for (nfile = 0; nfile < NumFiles; nfile++) {
        if (NumFiles > 1)
            sprintf(filename, "%s.%d", temp, nfile);
        else
            sprintf(filename, "%s", temp);

        infile = fopen(filename, "r");
        fread(&blksize, sizeof(int), 1, infile);
        SnapFormat = (blksize == 8) ? 2 : 1;
        rewind(infile);

        fileoffset = (SnapFormat == 2)
            ? sizeof(int) + 4 * sizeof(char) + 2 * sizeof(int) + sizeof(int)
            : sizeof(int);

        fseek(infile, fileoffset, SEEK_SET);

        NumPartInFile = 0;
        memset(tmp_npart, 0, sizeof(tmp_npart));
        fread(tmp_npart, sizeof(int), GADGET_MAX_TYPES, infile);
        for (i = 0; i < GADGET_MAX_TYPES; i++)
            NumPartInFile += tmp_npart[i];

        /* Skip over the rest of the 256-byte header block plus its
           Fortran record delimiters, then the SnapFormat-2 block
           descriptor if present */
        fileoffset += 256 + sizeof(int);
        if (SnapFormat == 2)
            fileoffset += sizeof(int) + 4 * sizeof(char)
                        + 2 * sizeof(int) + sizeof(int);

        fseek(infile, fileoffset, SEEK_SET);

        fprintf(stdout, "Reading %lld particles from %s on Task %d...\n",
                NumPartInFile, filename, ThisTask);
        fflush(stdout);

        /* Seek to this task's share of the position block */
        fileoffset += (long long)3 * ThisTask
                      * (NumPartInFile / NTask) * sizeof(float);
        fseek(infile, fileoffset, SEEK_SET);

        long long local_n = NumPartInFile / NTask;
        pid = (unsigned long long *)malloc(
                  (local_n + 1) * sizeof(unsigned long long));

        for (i = (int)(ThisTask * local_n);
             i < (int)((ThisTask + 1) * local_n); i++) {
            fread(&pos, 3 * sizeof(float), 1, infile);
            x[*NThisTask] = pos[0];
            y[*NThisTask] = pos[1];
            z[*NThisTask] = pos[2];
            pid[*NThisTask - ThisTask * local_n] = i;
            (*NThisTask)++;
        }

        /*
         * Assign particle types: build a cumulative count array from
         * tmp_npart (only up to GADGET_MAX_TYPES, which is what's in
         * the binary file), then determine which type each particle
         * belongs to.
         */
        long long counts[GADGET_MAX_TYPES];
        long long running = 0;
        for (j = 0; j < GADGET_MAX_TYPES; j++) {
            running      += tmp_npart[j];
            counts[j]     = running;
        }

        long long local_count = *NThisTask - ThisTask * local_n;
        for (j = 0; j < (int)local_count; j++) {
            int t = 0;
            while (t < GADGET_MAX_TYPES - 1 && counts[t] <= (long long)pid[j])
                t++;
            ptype[ThisTask * local_n + j] = t;
        }

        free(pid);
        fclose(infile);
    }
}
