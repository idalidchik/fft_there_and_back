// main.cpp : Defines the entry point for the console application.
//
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <assert.h>

#include "fft_real.h"
#include "fft_t.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

namespace test {

////////////////////////////////////////////////////////////////////////////////

struct signal
{
    template<typename T>
    static void example(T& data)
    {
        const double PI_2 = M_PI * 2;
        const int N = 8;
        const double rate = 8000;
        const double time_step = 1.0 / rate;
        
		data.clear();
        data.resize(N);

        for (int i = 0; i < N; ++i)
        {
            data[i] = std::sin(PI_2 * 1000.0 * i * time_step) +
                0.5 * std::sin(PI_2 * 2000.0 * i * time_step + (3.0 * M_PI / 4.0)); 
        }
    }

    /*template<typename T>
	static void print(T& data)
	{
		for (size_t i = 0; i < data.size(); ++i)
		{
			std::cout << "[" << i << "] = ";
			if (std::abs(data[i]) < 1e-6)
				std::cout << 0;
			else
				std::cout << data[i];
			std::cout << "\n";
		}
	}*/

    template<typename T>
	static void print(T& data)
	{
        size_t i = 0;
		for (const auto & d : data) {
			std::cout << "[" << i++ << "] = ";
			if (std::abs(d) < 1e-6)
				std::cout << 0;
			else
				std::cout << d;
			std::cout << "\n";
		}
	}

	template<typename T>
	static void realft(T& data, const int isign, const bool conjugate)
	{
        static_assert(std::is_same<float, typename T::value_type>::value, "");
		assert((isign == 1) || (isign == -1));

		float * const p = &(data[0]);

		if (conjugate && (isign == -1)) {
			numeric::conjugate(p, data.size());
		}
		numeric::realft(p-1, data.size(), isign);

		if (conjugate && (isign == 1)) {
			numeric::conjugate(p, data.size());
		}
		if (isign == -1) { // inverse transform
			const float factor = 2.0 / data.size();
			for (size_t i = 0; i < data.size(); ++i) {
				data[i] *= factor;
			}
		}
	}

	template<int isign, typename T>
	static void realft(T& data)
	{
		static_assert((isign == 1) || (isign == -1), "");
		realft(data, isign, true);
	}

};

} // namespace test

int main(int argc, char* argv[])
{
    if (1) {
      	std::vector<float> data;
        if (1) {
	        for (size_t i = 0; i < 2; ++i) 
	        {
		        const bool conjugate = (0 == i);

		        std::cout << "\nconjugate = " << conjugate << "\n";

		        test::signal::example(data);

		        std::cout << "\ntime domain:\n"; 
		        test::signal::print(data);

		        test::signal::realft(data, 1, conjugate);

		        std::cout << "\nfrequency domain:\n";
		        test::signal::print(data);

		        test::signal::realft(data, -1, conjugate);

		        std::cout << "\ntime domain:\n"; 
		        test::signal::print(data);

		        std::cout << "\n----------------------------------------------\n";
	        }
        }
        if (1) {
	        test::signal::example(data);

	        std::cout << "\ntime domain:\n"; 
	        test::signal::print(data);

	        test::signal::realft<1>(data);

	        std::cout << "\nfrequency domain:\n";
	        test::signal::print(data);

	        test::signal::realft<-1>(data);

	        std::cout << "\ntime domain:\n"; 
	        test::signal::print(data);
        }
    }
    if (1) {
	    enum { nn = 1 << 10 }; 
	    enum { isign = 1 };
	    float data[nn] = {}; // filled by 0

        std::cout << "\nnumeric::realft_t:\n"; 
        numeric::realft_t<nn, isign>(data - 1);
        //test::signal::print(data);
    }
    return 0;
}

