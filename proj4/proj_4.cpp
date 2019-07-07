#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>
#include "simd.p4.h"

// how many tries to discover the maximum performance:
#ifndef NUMTRIES
#define NUMTRIES	10
#endif

// setting the number of trials in the monte carlo simulation:
#ifndef ARR_SIZE
#define ARR_SIZE	1000
#endif

// function prototypes:
float		Ranf(float, float);
int			Ranf(int, int);
void		TimeOfDaySeed();
void		cpp_mul(float *, float *, float *, int len);

int main(int argc, char *argv[]) {
	TimeOfDaySeed();		// seed the random number generator

	// cpp mul
	{
		float *a = new float[ARR_SIZE];
		float *b = new float[ARR_SIZE];
		float *c = new float[ARR_SIZE];

		// fill the random-value arrays:
		for (int n = 0; n < ARR_SIZE; n++)
		{
			a[n] = Ranf(float(-1.0), float(1.0));
			b[n] = Ranf(float(-1.0), float(1.0));
		}

		float maxPerformance = 0.;

		for (int t = 0; t < NUMTRIES; t++) {
			double time0 = omp_get_wtime();

			//
			cpp_mul(a, b, c, ARR_SIZE);

			double time1 = omp_get_wtime();
			double megaTrialsPerSecond = (double)ARR_SIZE / (time1 - time0) / 1000000.;
			if (megaTrialsPerSecond > maxPerformance)
				maxPerformance = megaTrialsPerSecond;
		}

		delete[] a;
		delete[] b;
		delete[] c;

		printf("%f\t", maxPerformance);
	}

	// SIMD mul
	{
		float *a = new float[ARR_SIZE];
		float *b = new float[ARR_SIZE];
		float *c = new float[ARR_SIZE];

		// fill the random-value arrays:
		for (int n = 0; n < ARR_SIZE; n++)
		{
			a[n] = Ranf(float(-1.0), float(1.0));
			b[n] = Ranf(float(-1.0), float(1.0));
		}

		float maxPerformance = 0.;

		for (int t = 0; t < NUMTRIES; t++) {
			double time0 = omp_get_wtime();

			//
			SimdMul(a, b, c, ARR_SIZE);

			double time1 = omp_get_wtime();
			double megaTrialsPerSecond = (double)ARR_SIZE / (time1 - time0) / 1000000.;
			if (megaTrialsPerSecond > maxPerformance)
				maxPerformance = megaTrialsPerSecond;
		}

		delete[] a;
		delete[] b;
		delete[] c;

		printf("%f\t", maxPerformance);
	}

	// cpp mul + reduction
	{
		float *a = new float[ARR_SIZE];
		float *b = new float[ARR_SIZE];
		float *c = new float[ARR_SIZE];

		// fill the random-value arrays:
		for (int n = 0; n < ARR_SIZE; n++)
		{
			a[n] = Ranf(float(-1.0), float(1.0));
			b[n] = Ranf(float(-1.0), float(1.0));
		}

		float maxPerformance = 0.;

		for (int t = 0; t < NUMTRIES; t++) {
			double time0 = omp_get_wtime();

			//
			NonSimdMulSum(a, b, ARR_SIZE);

			double time1 = omp_get_wtime();
			double megaTrialsPerSecond = (double)ARR_SIZE / (time1 - time0) / 1000000.;
			if (megaTrialsPerSecond > maxPerformance)
				maxPerformance = megaTrialsPerSecond;
		}

		delete[] a;
		delete[] b;
		delete[] c;

		printf("%f\t", maxPerformance);
	}

	// SIMD mul + reduction
	{
		float *a = new float[ARR_SIZE];
		float *b = new float[ARR_SIZE];
		float *c = new float[ARR_SIZE];

		// fill the random-value arrays:
		for (int n = 0; n < ARR_SIZE; n++)
		{
			a[n] = Ranf(float(-1.0), float(1.0));
			b[n] = Ranf(float(-1.0), float(1.0));
		}

		float maxPerformance = 0.;

		for (int t = 0; t < NUMTRIES; t++) {
			double time0 = omp_get_wtime();

			//
			SimdMulSum(a, b, ARR_SIZE);

			double time1 = omp_get_wtime();
			double megaTrialsPerSecond = (double)ARR_SIZE / (time1 - time0) / 1000000.;
			if (megaTrialsPerSecond > maxPerformance)
				maxPerformance = megaTrialsPerSecond;
		}

		delete[] a;
		delete[] b;
		delete[] c;

		printf("%f\n", maxPerformance);
	}
}

float
Ranf(float low, float high)
{
	float r = (float)rand();				// 0 - RAND_MAX
	float t = r / (float)RAND_MAX;			// 0. - 1.

	return   low + t * (high - low);
}

int
Ranf(int ilow, int ihigh)
{
	float low = (float)ilow;
	float high = ceil((float)ihigh);

	return (int)Ranf(low, high);
}

void
TimeOfDaySeed()
{
	struct tm y2k = { 0 };
	y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

	time_t  timer;
	time(&timer);
	double seconds = difftime(timer, mktime(&y2k));
	unsigned int seed = (unsigned int)(1000.*seconds);    // milliseconds
	srand(seed);
}

void cpp_mul(float *a, float *b, float *c, int len){
	for (int i = 0; i < len; i++) {
		c[i] = a[i] * b[i];
	}
}