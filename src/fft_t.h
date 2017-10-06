// fft_t.h
//
#ifndef __FFT_T_H__
#define __FFT_T_H__

#pragma once

#include <math.h>
    
namespace numeric {

template<unsigned long x> struct is_pow_2 {
    enum { value = x && !(x & (x - 1)) };
};

/*
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; 
or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as -1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn. 
nn MUST be an integer power of 2 (this is not checked for!).
*/
template<unsigned long nn, int isign>
void four1_t(float data[])
{
    static_assert(is_pow_2<nn>::value, "four1_t");

	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta; // Double precision for the trigonometric recurrences.
	float tempr,tempi; 

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { // This is the bit-reversal section of the routine.
		if (j > i) {
            tempr=data[j];data[j]=data[i];data[i]=tempr;
            tempr=data[j+1];data[j+1]=data[i+1];data[i+1]=tempr;
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	// Here begins the Danielson-Lanczos section of the routine.
	mmax=2;
	while (n > mmax) { // Outer loop executed log2(nn) times.
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); // Initialize the trigonometric recurrence.
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { // Here are the two nested inner loops.
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; // This is the Danielson-Lanczos formula:
#pragma warning(disable:4244)
                tempr=wr*data[j]-wi*data[j+1]; 
				tempi=wr*data[j+1]+wi*data[j];
#pragma warning(default:4244)
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; // Trigonometric recurrence.
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
} // four1_t

template<unsigned long n, int isign>
void realft_t(float data[])
{
    static_assert(is_pow_2<n>::value, "realft_t");

	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta; // Double precision for the trigonometric recurrences.
	theta=3.141592653589793/(double) (n>>1); // Initialize the recurrence.
	if (isign == 1)	{
		c2 = -0.5;
		four1_t<n/2,1>(data); // The forward transform is here.
	} else {
		c2=0.5; // Otherwise set up for an inverse transform.
		theta = -theta; 
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) { // Case i=1 done separately below.
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]); // The two separate transforms are separated out of data.
		h1i=c1*(data[i2]-data[i4]); 
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
#pragma warning(disable:4244)
		data[i1]=h1r+wr*h2r-wi*h2i; // Here they are recombined to form the true transform of the original real data.
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
#pragma warning(default:4244)
		wr=(wtemp=wr)*wpr-wi*wpi+wr; // The recurrence.
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2]; // Squeeze the first and last data together to get them all within the	original array.
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1_t<n/2,-1>(data); // This is the inverse transform for the case isign=-1.
	} 
} // realft_t

} // namespace numeric

#endif //__FFT_T_H__