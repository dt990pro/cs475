#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>

#define XMIN	 0.
#define XMAX	 3.
#define YMIN	 0.
#define YMAX	 3.

#define TOPZ00  0.
#define TOPZ10  1.
#define TOPZ20  0.
#define TOPZ30  0.

#define TOPZ01   1.
#define TOPZ11  12.
#define TOPZ21   1.
#define TOPZ31   0.

#define TOPZ02  0.
#define TOPZ12  1.
#define TOPZ22  0.
#define TOPZ32  4.

#define TOPZ03  3.
#define TOPZ13  2.
#define TOPZ23  3.
#define TOPZ33  3.

#define BOTZ00  0.
#define BOTZ10  -3.
#define BOTZ20  0.
#define BOTZ30  0.

#define BOTZ01  -2.
#define BOTZ11  10.
#define BOTZ21  -2.
#define BOTZ31  0.

#define BOTZ02  0.
#define BOTZ12  -5.
#define BOTZ22  0.
#define BOTZ32  -6.

#define BOTZ03  -3.
#define BOTZ13   2.
#define BOTZ23  -8.
#define BOTZ33  -3.

// setting the number of threads:
#ifndef NUMNODES
#define NUMNODES	4
#endif

// how many tries to discover the maximum performance:
#ifndef NUMTRIES
#define NUMTRIES	10
#endif

// setting the number of threads:
#ifndef NUMT
#define NUMT		4
#endif

float Height(int, int);

int main(int argc, char *argv[])
{
#ifndef _OPENMP
	fprintf(stderr, "No OpenMP support!\n");
	return 1;
#endif

	omp_set_num_threads(NUMT);	// set the number of threads to use in the for-loop:`

	// the area of a single full-sized tile:
	float fullTileArea = (((XMAX - XMIN) / (float)(NUMNODES - 1))  *
	((YMAX - YMIN) / (float)(NUMNODES - 1)));

	// get ready to record the maximum performance and the volume:
	float maxPerformance = 0.;
	float v = 0.;

	// sum up the weighted heights into the variable "volume"
	// using an OpenMP for loop and a reduction:
	// looking for the maximum performance:
	for (int t = 0; t < NUMTRIES; t++)
	{
		double time0 = omp_get_wtime();

		#pragma omp parallel for default(none) shared(fullTileArea) reduction(+:v)
		for (int i = 0; i < NUMNODES*NUMNODES; i++)
		{
			int iu = i % NUMNODES;
			int iv = i / NUMNODES;

			bool iu_on_edge = iu == 0 || iu == NUMNODES - 1;  // iu on the edge
			bool iv_on_edge = iv == 0 || iv == NUMNODES - 1;  // iv on the edge

			// on corner
			if (iu_on_edge && iv_on_edge)
				v += Height(iu, iv) * (fullTileArea / 4.);
			// on edge
			else if (iu_on_edge || iv_on_edge)
				v += Height(iu, iv) * (fullTileArea / 2.);
			// inside
			else
				v += Height(iu, iv) * fullTileArea;
		}

		double time1 = omp_get_wtime();
		double megaTrialsPerSecond = (double)(NUMNODES*NUMNODES) / (time1 - time0) / 1000000.;
		if (megaTrialsPerSecond > maxPerformance)
			maxPerformance = megaTrialsPerSecond;
	}

	// Print out: (1) the number of threads, (2) the number of trials, (3) the volume, (4) the MegaTrialsPerSecond.
	printf("%d\n%d\n%f\n%f\n", NUMT, NUMNODES, v / NUMTRIES, maxPerformance);
}

float
Height(int iu, int iv)	// iu,iv = 0 .. NUMNODES-1
{
	float u = (float)iu / (float)(NUMNODES - 1);
	float v = (float)iv / (float)(NUMNODES - 1);

	// the basis functions:

	float bu0 = (1. - u) * (1. - u) * (1. - u);
	float bu1 = 3. * u * (1. - u) * (1. - u);
	float bu2 = 3. * u * u * (1. - u);
	float bu3 = u * u * u;

	float bv0 = (1. - v) * (1. - v) * (1. - v);
	float bv1 = 3. * v * (1. - v) * (1. - v);
	float bv2 = 3. * v * v * (1. - v);
	float bv3 = v * v * v;

	// finally, we get to compute something:


	float top = bu0 * (bv0*TOPZ00 + bv1*TOPZ01 + bv2*TOPZ02 + bv3*TOPZ03)
		+ bu1 * (bv0*TOPZ10 + bv1*TOPZ11 + bv2*TOPZ12 + bv3*TOPZ13)
		+ bu2 * (bv0*TOPZ20 + bv1*TOPZ21 + bv2*TOPZ22 + bv3*TOPZ23)
		+ bu3 * (bv0*TOPZ30 + bv1*TOPZ31 + bv2*TOPZ32 + bv3*TOPZ33);

	float bot = bu0 * (bv0*BOTZ00 + bv1*BOTZ01 + bv2*BOTZ02 + bv3*BOTZ03)
		+ bu1 * (bv0*BOTZ10 + bv1*BOTZ11 + bv2*BOTZ12 + bv3*BOTZ13)
		+ bu2 * (bv0*BOTZ20 + bv1*BOTZ21 + bv2*BOTZ22 + bv3*BOTZ23)
		+ bu3 * (bv0*BOTZ30 + bv1*BOTZ31 + bv2*BOTZ32 + bv3*BOTZ33);

	return top - bot;	// if the bottom surface sticks out above the top surface
						// then that contribution to the overall volume is negative
}