#ifndef _FFT_H_
#define _FFT_H_


#include <complex>
#include <iostream>
#include <valarray>

using namespace std;

// the fft and the ifft are from https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B

const double PI = 3.141592653589793238460;

typedef complex<double> Complex;
typedef valarray<Complex> CArray;

void fft(CArray &x);
void ifft(CArray& x);
void test();

#endif //_FFT_H_