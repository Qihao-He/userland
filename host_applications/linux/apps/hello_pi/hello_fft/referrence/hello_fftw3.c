/*
Created: 2/17/2018
desciption:
hello_fftw3.c with function of auto size and output to .csv
Author:Qihao He
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "fftw3.h"
// #include <unistd.h>
#define REAL 0
#define IMAG 1

char Usage[] =
    "Usage: hello_fftw3.bin log2_N [log2_M [loops [RMS_C]]]\n"
    "log2_N = log2(FFT_length),       log2_N = 8...22\n"
    "log2_M = log2(FFT_length),       log2_M > log2_N\n"
    "loops  = number of test repeats, loops>0,       default 1\n"
    "RMS_C  = number of test repeats, T(1),F(0),     default 1\n";

// not sure about the time profiling for the RPI function
unsigned Microseconds(void);
void REL_RMS_ERR_init(int span_log2_N, int loops, double **REL_RMS_ERR);
void time_elapsed_init(int span_log2_N, int loops);
void input_buffer(fftw_complex* in, int N); // input buffer
// output REL_RMS_ERR
void output_RMS(fftw_complex *out, int span_log2_N, double **REL_RMS_ERR, int N,
   int j, int k);
   // print out REL_RMS_ERR
void print_RMS(int span_log2_N, int loops, int log2_N, double **REL_RMS_ERR);

/* global array for holding the RMS error */
// double REL_RMS_ERR[span_log2_N][loops]; //2D array

int main(int argc, char *argv[]){
    int i, j, k, l, loops, freq, log2_N, log2_M, log2_P, N, RMS_C, span_log2_N;
    double **REL_RMS_ERR;
    unsigned t[4];

    fftw_complex *in, *out; //in, out buffer
    fftw_plan p; //fftw_plan prepare

    log2_N = argc>1? atoi(argv[1]) : 12; // 8 <= log2_N <= 22
    log2_M = argc>2? atoi(argv[2]) : log2_N + 1; // 8 <= log2_N <= 22
    loops  = argc>3? atoi(argv[3]) : 1;  // test repetitions
    RMS_C  = argc>4? atoi(argv[4]) : 1;  // RMS_controller

    if (!(argc>=2 && argc<=5) || loops<1 || !(RMS_C>=0 && RMS_C<=1) ||
    log2_N>=log2_M) {
        printf(Usage);
        return -1;
    }

    span_log2_N = log2_M - log2_N;
    REL_RMS_ERR=(double**)malloc(span_log2_N*sizeof(double *));
    if(REL_RMS_ERR==NULL)
    {
      printf("Malloc failed\n");
      exit(-1);
    }
    for (i=0;i<span_log2_N;i++)
    {
          REL_RMS_ERR[i]=(double *)malloc(loops*sizeof(double));
          if(REL_RMS_ERR[i]==NULL)
          {
             printf("Malloc failed on loop %d",i);
             exit(-1);
          }
    }

// initializing 2D, 3D array to 0
    REL_RMS_ERR_init(span_log2_N, loops, (double **)REL_RMS_ERR);
    time_elapsed_init(span_log2_N, loops);
// print out lables for .csv file
    printf("log2_N,N,Init_T,FFT_T,RMS_T,Total_T\n");

    for(l=0; l<span_log2_N; l++){
        log2_P = log2_N + l;
        N = 1<<log2_P; // initializing FFT length: N
        in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

        for (k=0; k<loops; k++) {
            t[0] = Microseconds();
            input_buffer(in, N);

            t[1] = Microseconds();
            fftw_execute(p); /* repeat as needed */
            t[2] = Microseconds();

            if(RMS_C == 1) output_RMS(out, span_log2_N, REL_RMS_ERR, N, l, k);
            t[3] = Microseconds();
            printf("%i,%i,%d,%d,%d,%d\n",log2_P,N,t[1] - t[0],t[2] - t[1],
            t[3] - t[2],t[3] - t[0]); //print for .csv file
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
    }
    // print out REL_RMS_ERR
    if(RMS_C == 1) print_RMS(span_log2_N, loops, log2_N, REL_RMS_ERR);
    return 0;
}

unsigned Microseconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec*1000000 + ts.tv_nsec/1000;
}

void REL_RMS_ERR_init(int span_log2_N, int loops, double **REL_RMS_ERR){
    int i, j;
    for(i=0; i<span_log2_N; i++){
        for(j=0; j<loops; j++){
            REL_RMS_ERR[i][j] = 0;
        }
    }
}

void time_elapsed_init(int span_log2_N, int loops){
    int i,j,k;
    double time_elapsed[span_log2_N][loops][4]; //3D array
    for(i=0; i<span_log2_N; i++){
        for(j=0; j<loops; j++){
            for(k=0; k<4; k++){
                time_elapsed[i][j][k] = 0;
            }
        }
    }
}
// input buffer
void input_buffer(fftw_complex *in, int N){
    int i;
    for (i = 0; i < N; i++) in[i][REAL] = in[i][IMAG] = 0;
    in[1][REAL] = in[N-1][REAL] = 0.5;
}

// output REL_RMS_ERR
void output_RMS(fftw_complex *out, int span_log2_N, double **REL_RMS_ERR, int N,
   int j, int k){
    int i;
    double tsq[2], a;
    tsq[0]=tsq[1]=0;
    a = 2 * M_PI / N;
    for (i=0; i<N; i++) {
        double re = cos(a * i);
        tsq[0] += pow(re, 2);
        tsq[1] += pow(re - out[i][REAL], 2) + pow(out[i][IMAG], 2);
    }
    REL_RMS_ERR[j][k] = sqrt(tsq[1] / tsq[0]);
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
