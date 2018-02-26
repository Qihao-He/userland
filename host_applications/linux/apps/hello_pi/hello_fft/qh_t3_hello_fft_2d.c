/*
BCM2835 "GPU_FFT" release 2.0
Copyright (c) 2014, Andrew Holme.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/*
Created: 2/24/2018
desciption:
hello_fft with function of auto size and output to .csv
Author:Qihao He
*/

#include <string.h>
#include <stdio.h>
#include <time.h>

#include "gpu_fft_trans.h"
#include "hello_fft_2d_bitmap.h"

#define GPU_FFT_ROW(fft, io, y) ((fft)->io+(fft)->step*(y))

char Usage[] =
    "Usage: hello_fft.bin log2_N [jobs [loops [RMS_C]]]\n"
    "log2_N = log2(FFT_length),       log2_N = 8...11\n"
    "log2_M = log2(FFT_length),       log2_M > log2_N\n"
    "loops  = number of test repeats, loops>0,       default 1\n"
    "RMS_C  = RMS_controller, T(1),F(0),     default 0\n"
    "BMP_C  = BMP_controller, T(1),F(0),       default 0\n";

unsigned Microseconds(void);
void REL_RMS_ERR_init(int span_log2_N, int loops, double **REL_RMS_ERR);
void output_RMS(struct GPU_FFT *fft, struct GPU_FFT_COMPLEX *base,
  int span_log2_N, double **REL_RMS_ERR, int N, int j, int k);
void print_RMS(int span_log2_N, int loops, int log2_N, double **REL_RMS_ERR);
// void time_elapsed_init(int span_log2_N, int loops);

int main(int argc, char *argv[]) {
    int x, y, l, k, ret, mb = mbox_open(), log2_N, log2_M, log2_P, span_log2_N,
     loops, N, RMS_C, BMP_C;
    unsigned t[4];
    double tsq[2];

    struct GPU_FFT_COMPLEX *row;
    struct GPU_FFT_TRANS *trans;
    struct GPU_FFT *fft_pass[2];

    log2_N = argc>1? atoi(argv[1]) : 12; // 8 <= log2_N <= 11
    log2_M = argc>2? atoi(argv[2]) : log2_N + 1; // 8 <= log2_N <= 11
    loops  = argc>3? atoi(argv[3]) : 1;  // test repetitions
    RMS_C  = argc>4? atoi(argv[4]) : 0;  // RMS_controller
    BMP_C  = argc>5? atoi(argv[5]) : 0;  // BMP_controller
    //
    if (!(argc >=2 && argc <= 6) || loops < 1 || log2_N >= log2_M ||
    !(log2_N >= 8 && log2_N <= 11 && log2_M <= 12) ||
    !(RMS_C >= 0 && RMS_C <= 1) || !(BMP_C >= 0 && BMP_C <= 1)){
        printf(Usage);
        return -1;
    }

    span_log2_N = log2_M - log2_N;
    REL_RMS_ERR = (double **)malloc(span_log2_N * sizeof(double *));
    if(REL_RMS_ERR == NULL){
      printf("Malloc failed\n");
      exit(-1);
    }
    for (i = 0; i < span_log2_N; i++){
          REL_RMS_ERR[i] = (double *)malloc(loops * sizeof(double));
          if(REL_RMS_ERR[i] == NULL){
             printf("Malloc failed on loop %d",i);
             exit(-1);
          }
    }
    // initializing 2D, 3D array to 0
    REL_RMS_ERR_init(span_log2_N, loops, (double **)REL_RMS_ERR);
    // time_elapsed_init(span_log2_N, loops);
    // print out lables for .csv file
    printf("log2_N,N,1st FFT_T,Transpose_T,2nd FFT_T\n");

    for(l = 0; l < span_log2_N; l++){
        log2_P = log2_N + l;
        N = 1 << log2_P; // FFT length

        if (BMP_C == 1){
            BITMAPFILEHEADER bfh;
            BITMAPINFOHEADER bih;

            // Create Windows bitmap file
            FILE *fp = fopen("hello_fft_2d.bmp", "wb");
            if (!fp) return -666;

            // Write bitmap header
            memset(&bfh, 0, sizeof(bfh));
            bfh.bfType = 0x4D42; //"BM"
            bfh.bfSize = N*N*3;
            bfh.bfOffBits = sizeof(bfh) + sizeof(bih);
            fwrite(&bfh, sizeof(bfh), 1, fp);

            // Write bitmap info
            memset(&bih, 0, sizeof(bih));
            bih.biSize = sizeof(bih);
            bih.biWidth = N;
            bih.biHeight = N;
            bih.biPlanes = 1;
            bih.biBitCount = 24;
            bih.biCompression = BI_RGB;
            fwrite(&bih, sizeof(bih), 1, fp);
        }

        // Prepare 1st FFT pass
        ret = gpu_fft_prepare(mb, log2_P, GPU_FFT_REV, N, fft_pass+0);
        if (ret) {
            return ret;
        }
        // Prepare 2nd FFT pass
        ret = gpu_fft_prepare(mb, log2_P, GPU_FFT_REV, N, fft_pass+1);
        if (ret) {
            gpu_fft_release(fft_pass[0]);
            return ret;
        }
        // Transpose from 1st pass output to 2nd pass input
        ret = gpu_fft_trans_prepare(mb, fft_pass[0], fft_pass[1], &trans);
        if (ret) {
            gpu_fft_release(fft_pass[0]);
            gpu_fft_release(fft_pass[1]);
            return ret;
        }

        for (k = 0; k < loops; k++) {
            // Clear input array
            for (y=0; y<N; y++) {
                row = GPU_FFT_ROW(fft_pass[0], in, y);
                for (x=0; x<N; x++) row[x].re = row[x].im = 0;
            }

            // Setup input data
            GPU_FFT_ROW(fft_pass[0], in, 0)[0].re = 1;

            // ==> FFT() ==> T() ==> FFT() ==>
            usleep(1); /* yield to OS */   t[0] = Microseconds();
            gpu_fft_execute(fft_pass[0]);  t[1] = Microseconds();
            gpu_fft_trans_execute(trans);  t[2] = Microseconds();
            gpu_fft_execute(fft_pass[1]);  t[3] = Microseconds();


            if (RMS_C == 1){
              output_RMS(fft, base, span_log2_N, REL_RMS_ERR, N, l, k);
            }

            printf( "%i,%i,%d,%d,%d\n", log2_P, N, t[3] - t[2], t[2] - t[1],
            t[1] - t[0]);
        }
        // print output
        for (y = 0; y < N; y ++) {
            row = GPU_FFT_ROW(fft_pass[1], out, y);
            for (x = 0; x < N; x ++) {
                printf("value is %lf + j%lf\n",row[x].re,row[x].im);
            }
        }
        // Write output to bmp file
        if (BMP_C == 1){
            for (y = 0; y < N; y ++) {
                row = GPU_FFT_ROW(fft_pass[1], out, y);
                for (x = 0; x < N; x ++) {
                    fputc(128 + row[x].re, fp); // blue
                    fputc(128 + row[x].re, fp); // green
                    fputc(128 + row[x].re, fp); // red
                }
            }
        }
        // Clean-up properly.  Videocore memory lost if not freed !
        gpu_fft_release(fft_pass[0]);
        gpu_fft_release(fft_pass[1]);
        gpu_fft_trans_release(trans);
    }
    if (RMS_C == 1) print_RMS(span_log2_N, loops, log2_N, REL_RMS_ERR);
    return 0;
}

unsigned Microseconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec*1000000 + ts.tv_nsec/1000;
}

void REL_RMS_ERR_init(int span_log2_N, int loops, double **REL_RMS_ERR){
    int i, j;
    for(i = 0; i < span_log2_N; i++){
        for(j = 0; j < loops; j++){
            REL_RMS_ERR[i][j] = 0;
        }
    }
}

// output REL_RMS_ERR
void output_RMS(struct GPU_FFT *fft, struct GPU_FFT_COMPLEX *base,
  int span_log2_N, double **REL_RMS_ERR, int N, int l, int k){
    int i, j, freq;
    double tsq[2], a;
    tsq[0] = tsq[1] = 0;
    a = 2 * GPU_FFT_PI / N;
    base = fft->out + fft->step;
    for (i = 0; i < N; i++) {
        double re = cos(a * i);
        tsq[0] += pow(re, 2);
        tsq[1] += pow(re - base[i].re, 2) + pow(base[i].im, 2);
    }
    REL_RMS_ERR[l][k] = sqrt(tsq[1] / tsq[0]);
}
// print out REL_RMS_ERR
void print_RMS(int span_log2_N, int loops, int log2_N, double **REL_RMS_ERR){
    int i,j;
    for (i = 0; i < span_log2_N; i++) {
        printf("REL_RMS_ERR for log2_N:%d\n", log2_N + i);
            for (j = 0; j < loops; j++) {
                printf("%.10e,",REL_RMS_ERR[i][j]);
            }
            printf("\n");
    }
    printf("\n");
}
//
// void time_elapsed_init(int span_log2_N, int loops){
//     int i,j,k;
//     double time_elapsed[span_log2_N][loops][4]; //3D array
//     for(i = 0; i < span_log2_N; i++){
//         for(j = 0; j < loops; j++){
//             for(k = 0; k < 4; k++){
//                 time_elapsed[i][j][k] = 0;
//             }
//         }
//     }
// }
