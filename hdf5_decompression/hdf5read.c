/* Copyright (C) 2017-2018 Intel Corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted for any purpose (including commercial purposes)
 * provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions, and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions, and the following disclaimer in the
 *    documentation and/or materials provided with the distribution.
 *
 * 3. In addition, redistributions of modified forms of the source or binary
 *    code must carry prominent notices stating that the original code was
 *    changed and the date of the change.
 *
 *  4. All publications or advertising materials mentioning features or use of
 *     this software are asked, but not required, to acknowledge that it was
 *     developed by Intel Corporation and credit the contributors.
 *
 * 5. Neither the name of Intel Corporation, nor the name of any Contributor
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <hdf5.h>

#define BATCH_SIZE 128

void timespec_diff(struct timespec *start,
                   struct timespec *stop,
                   struct timespec *result)
{
       if ((stop->tv_nsec - start->tv_nsec) < 0) {
           result->tv_sec = stop->tv_sec - start->tv_sec - 1;
           result->tv_nsec = stop->tv_nsec - start->tv_nsec + 1000000000;
       } else {
           result->tv_sec = stop->tv_sec - start->tv_sec;
           result->tv_nsec = stop->tv_nsec - start->tv_nsec;
       }
}

int main(int argc, char **argv)
{
//    char *path = "/global/cscratch1/sd/afarbin/"
//                 "h5_files_2D_3D/2D_h5/antielectron_207.2d.h5";
    char *path = "antielectron_207.2d_raw.h5";
    int i, eindex = 0, rank;
    hsize_t dims[4], offsets[] = {0, 0, 0, 0};
    hsize_t counts[] = {BATCH_SIZE, 2, 240, 4096};
    hid_t fid, dsid, dspid, mspid, tid;
    herr_t status;
    void *dsmem = NULL;
    struct timespec tstart, tend, tres;

    if(MPI_SUCCESS != MPI_Init(&argc, &argv))
        return -1;
    dsmem = malloc(2 * 2 * 240 * 4096L * BATCH_SIZE);
     if (dsmem == NULL)
        goto cleanup;
    fid = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    dsid = H5Dopen(fid, "/images", H5P_DEFAULT);
    tid = H5Dget_type(dsid);
    dspid = H5Dget_space(dsid);
    mspid = H5Screate_simple(4, counts, NULL);
    H5Sget_simple_extent_dims(dspid, dims, NULL);
    fprintf(stdout, "# of Samples: %d\n", dims[0]);
    for (i=0; i<(dims[0]/BATCH_SIZE); i++) {
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        offsets[0] = i * (hsize_t) BATCH_SIZE;
        H5Sselect_hyperslab(dspid, H5S_SELECT_SET, offsets, NULL, counts, NULL);
        status = H5Dread(dsid, tid, mspid, dspid, H5P_DEFAULT, dsmem);
        clock_gettime(CLOCK_MONOTONIC, &tend);
        timespec_diff(&tstart, &tend, &tres);
        fprintf(stdout, "Batch %d, Time: %lld.%.9ld\n", i,
                        (long long)tres.tv_sec, tres.tv_nsec);
        fflush(stdout);
    }
    H5Sclose(mspid);
    H5Tclose(tid);
    H5Sclose(dspid);
    H5Dclose(dsid);
    H5Fclose(fid);

cleanup:
    if (dsmem)
        free (dsmem);
    MPI_Finalize();
    return 0;
}
