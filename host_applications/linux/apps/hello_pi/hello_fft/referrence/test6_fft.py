#!/usr/bin/python
"""
Created: 2/14/2018
desciption:
scipy doing FFT and calcualte Relative RMS
Author:Qihao He
"""
# import libraries
import sys
import numpy as np
import scipy as sp
import math
from scipy.fftpack import ifft, fft
import time

# Usage
Usage = """Usage: hello_scipy_fft.py log2_N [log2_M [loops [RMS_C]]]
        log2_N = log2(FFT_length),       log2_N = 1...28
        log2_M = log2(FFT_length),       log2_M > log2_N
        loops  = number of test repeats, loops>0,       default 1
        RMS_C = True(1), False(0),    default 0"""

# Default values for optional arguments
log2_N = 8 #default value to be 8
if len(sys.argv) > 1:
    log2_N = int(sys.argv[1])

log2_M = log2_N + 1 #default value to be 9
if len(sys.argv) > 2:
    log2_M = int(sys.argv[2])

loops = 1 #default value to be 1
if len(sys.argv) > 3:
    loops = int(sys.argv[3])

RMS_C = 0
if len(sys.argv) > 4:
    RMS_C = int(sys.argv[4])

if (not 2 <= len(sys.argv) <= 5 or log2_M <= log2_N  or loops < 1 or not
0 <= RMS_C <= 1):
    print(Usage)
    sys.exit()

print "The jobsize for the scipy-FFT is 2^", log2_N
print "The jobsize_end for the scipy-FFT is 2^", log2_M
print "Repeat times:", loops
print "REL_RMS_ERR", RMS_C

# array of the log2_FFT_length
span_N = log2_M - log2_N
log2_FFT_length = np.zeros(span_N, dtype = np.int) # 1D span_Nrray
if RMS_C == 1:
    REL_RMS_ERR = np.zeros((loops, span_N), dtype = np.float64) # 2D array
time_elapsed = np.zeros((loops, span_N, 4), dtype = np.float64) # 3D array

for i in range(0, span_N): log2_FFT_length[i] = i + log2_N
print "log2_N,","N,""Init_T:,""FFT_T:,","RMS_T:,""Total_T"
for k in range(loops):
    for j in range(span_N):
        t0 = time.time()# Time counter

        N = 1 << int(log2_FFT_length[j]) #fft length
        # input buffer
        x = np.zeros((N, ), dtype = np.complex64)
        x.real[1] = x.real[N - 1] = np.float32(0.5)
        tsq0 = 0
        tsq1 = 0

        t1 = time.time()
        # fft execute
        y = fft(x)
        t2 = time.time()
        # output buffer and rel_rms_err
        if RMS_C == 1:
            l = 2 * math.pi / N
            for i in range(N):
                re = np.cos(l * i) # True solution
                tsq0 += re * re
                a = re - y.real[i]
                b = y.imag[i]
                tsq1 += a * a + b * b
            REL_RMS_ERR[k][j] = math.sqrt(tsq1 / tsq0)

        t3 = time.time()
        print  j + log2_N,",",N,",",t1 - t0,",",t2 - t1,",",t3 - t2,",",t3 - t0
        time_elapsed[k][j][0] = t1 - t0
        time_elapsed[k][j][1] = t2 - t1
        time_elapsed[k][j][2] = t3 - t2
        time_elapsed[k][j][3] = t3 - t0
if RMS_C == 1:
    print"rel_rms_err = ", REL_RMS_ERR
print"time_elapsed = ", time_elapsed
