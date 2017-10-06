// fft_real.cpp : 
//
#include "fft_real.h"
#include <math.h>
#include <assert.h>

// warning C4244: '=' : conversion from 'double' to 'float'
#pragma warning(disable:4244)

//http://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf
//http://en.wikipedia.org/wiki/Numerical_Recipes

namespace numeric {

#define NUMBER_IS_2_POW_K(x) ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...

#define FAST_SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/*
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; 
or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as -1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn. 
nn MUST be an integer power of 2 (this is not checked for!).
*/
void four1(float data[], unsigned long nn, int isign)
{
    assert(NUMBER_IS_2_POW_K(nn));

	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta; // Double precision for the trigonometric recurrences.
	float tempr,tempi; 

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { // This is the bit-reversal section of the routine.
		if (j > i) {
			FAST_SWAP(data[j],data[i]); // Exchange the two complex numbers.
			FAST_SWAP(data[j+1],data[i+1]);
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
				tempr=wr*data[j]-wi*data[j+1]; 
				tempi=wr*data[j+1]+wi*data[j];
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
} // four1

/*
Calculates the Fourier transform of a set of n real-valued data points. Replaces this data (which
is stored in array data[1..n]) by the positive frequency half of its complex Fourier transform.
The real-valued first and last components of the complex transform are returned as elements
data[1] and data[2], respectively. n must be a power of 2. This routine also calculates the
inverse transform of a complex data array if it is the transform of real data. (Result in this case
must be multiplied by 2/n.)
*/
void realft(float data[], unsigned long n, int isign)
{
    assert(NUMBER_IS_2_POW_K(n));

	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta; // Double precision for the trigonometric recurrences.
	theta=3.141592653589793/(double) (n>>1); // Initialize the recurrence.
	if (isign == 1)	{
		c2 = -0.5;
		four1(data,n>>1,1); // The forward transform is here.
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
		data[i1]=h1r+wr*h2r-wi*h2i; // Here they are recombined to form the true transform of the original real data.
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr; // The recurrence.
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2]; // Squeeze the first and last data together to get them all within the	original array.
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1); // This is the inverse transform for the case isign=-1.
	} 
} // realft

/*
Given two real input arrays data1[1..n] and data2[1..n], this routine calls four1 and
returns two complex output arrays, fft1[1..2n] and fft2[1..2n], each of complex length
n (i.e., real length 2*n), which contain the discrete Fourier transforms of the respective data
arrays. n MUST be an integer power of 2.
*/
void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned long n)
{
    assert(NUMBER_IS_2_POW_K(n));

    unsigned long nn3,nn2,jj,j;
    float rep,rem,aip,aim;

    nn3=1+(nn2=2+n+n);
    for (j=1,jj=2;j<=n;j++,jj+=2) { //Pack the two real arrays into one complex array.
        fft1[jj-1]=data1[j]; 
        fft1[jj]=data2[j];
    }
    four1(fft1,n,1); // Transform the complex array.
    fft2[1]=fft1[2];
    fft1[2]=fft2[2]=0.0;
    for (j=3;j<=n+1;j+=2) {
        rep=0.5*(fft1[j]+fft1[nn2-j]); // Use symmetries to separate the two transforms.
        rem=0.5*(fft1[j]-fft1[nn2-j]); 
        aip=0.5*(fft1[j+1]+fft1[nn3-j]);
        aim=0.5*(fft1[j+1]-fft1[nn3-j]);
        fft1[j]=rep; // Ship them out in two complex arrays.
        fft1[j+1]=aim;
        fft1[nn2-j]=rep;
        fft1[nn3-j] = -aim;
        fft2[j]=aip;
        fft2[j+1] = -rem;
        fft2[nn2-j]=aip;
        fft2[nn3-j]=rem;
    }
} // twofft

} // namespace numeric
