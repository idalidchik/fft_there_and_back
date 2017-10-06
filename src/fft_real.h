// fft_real.h
//
#ifndef __FFT_REAL_H__
#define __FFT_REAL_H__

#pragma once

namespace numeric {

/*
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; 
or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as -1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn. 
nn MUST be an integer power of 2 (this is not checked for!).
*/
void four1(float data[], unsigned long nn, int isign);

/*
Calculates the Fourier transform of a set of n real-valued data points. Replaces this data (which
is stored in array data[1..n]) by the positive frequency half of its complex Fourier transform.
The real-valued first and last components of the complex transform are returned as elements
data[1] and data[2], respectively. n must be a power of 2. This routine also calculates the
inverse transform of a complex data array if it is the transform of real data. (Result in this case
must be multiplied by 2/n.)
*/
void realft(float data[], unsigned long n, int isign);

/*
Given two real input arrays data1[1..n] and data2[1..n], this routine calls four1 and
returns two complex output arrays, fft1[1..2n] and fft2[1..2n], each of complex length
n (i.e., real length 2*n), which contain the discrete Fourier transforms of the respective data
arrays. n MUST be an integer power of 2.
*/
void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned long n);

//http://en.wikipedia.org/wiki/Complex_conjugate
// data is a complex array of length nn or, equivalently, a real array of length 2*nn. 
inline void conjugate(float data[], unsigned long nn)
{
#if 0
    for (unsigned long i = 0; i < nn; ++i)
    {
        data[2*i+1] *= -1;
    }
#else
	const unsigned long imax = (nn << 1) + 1;
	for (unsigned long i = 1; i < imax; i += 2)
	{
		data[i] *= -1;
	}
#endif
}

template<typename T, size_t nn_2>
inline void conjugate_t(T(&data)[nn_2])
{
	T * p = data + 1;
	T * end = p + nn_2;
	for (; p != end; p += 2) { //assert(p < data + nn_2);
		(*p) *= -1;
	}
}

} // namespace numeric

#endif //__FFT_REAL_H__
