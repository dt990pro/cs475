// Array multiplication: C = A * B:

// System includes
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

// CUDA runtime
#include <cuda_runtime.h>

// Helper functions and utilities to work with CUDA
#include "helper_functions.h"
#include "helper_cuda.h"


#ifndef BLOCKSIZE
#define BLOCKSIZE		32		// number of threads per block (16, 32, 64)
#endif

#ifndef SIZE
#define SIZE			1*1024*1024	// array size(16K, 32K, 64K, 128K, 256K, and 512K)
#endif

#ifndef NUMTRIALS
#define NUMTRIALS		100
#endif

#ifndef TOLERANCE
#define TOLERANCE		0.00001f	// tolerance to relative error
#endif

// ranges for the random numbers:
const float XCMIN = 0.;
const float XCMAX = 2.0;
const float YCMIN = 0.0;
const float YCMAX = 2.0;
const float RMIN = 0.5;
const float RMAX = 2.0;

__host__ void
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

__host__ float
Ranf(float low, float high)
{
	float r = (float)rand();				// 0 - RAND_MAX
	float t = r / (float)RAND_MAX;			// 0. - 1.

	return   low + t * (high - low);
}

// array multiplication (CUDA Kernel) on the device: C = A * B

__global__  void Monte_Carlo( float *A, float *B, float *R, float *C )
{
	__shared__ float prods[BLOCKSIZE];

	unsigned int numItems = blockDim.x;
	unsigned int tnum = threadIdx.x;
	unsigned int wgNum = blockIdx.x;
	unsigned int gid = blockIdx.x*blockDim.x + threadIdx.x;

	/*prods[tnum] = A[gid] * B[gid] * R[gid];*/

	// solve for the intersection using the quadratic formula:
	float a = 2.;
	float b = -2.*(A[gid] + B[gid]);
	float c = A[gid] * A[gid] + B[gid] * B[gid] - R[gid] * R[gid];
	float d = b*b - 4.*a*c;

	// If d is less than 0., then the circle was completely missed. (Case A) no hit.
	if (d < 0) prods[tnum] = 0;
	else {
		// hits the circle:
		// get the first intersection:
		d = sqrt(d);
		float t1 = (-b + d) / (2.*a);		// time to intersect the circle
		float t2 = (-b - d) / (2.*a);		// time to intersect the circle
		float tmin = t1 < t2 ? t1 : t2;		// only care about the first intersection

		// If tmin is less than 0., then the circle completely engulfs the laser pointer. (Case B) no hit.
		if (tmin < 0) prods[tnum] = 0;
		else {
			// where does it intersect the circle?
			float xcir = tmin;
			float ycir = tmin;

			// get the unitized normal vector at the point of intersection:
			float nx = xcir - A[gid];
			float ny = ycir - B[gid];
			float n = sqrt(nx*nx + ny*ny);
			nx /= n;	// unit vector
			ny /= n;	// unit vector

			// get the unitized incoming vector:
			float inx = xcir - 0.;
			float iny = ycir - 0.;
			float in = sqrt(inx*inx + iny*iny);
			inx /= in;	// unit vector
			iny /= in;	// unit vector

			// get the outgoing (bounced) vector:
			float dot = inx*nx + iny*ny;
			//float outx = inx - 2.*nx*dot;	// angle of reflection = angle of incidence`
			float outy = iny - 2.*ny*dot;	// angle of reflection = angle of incidence`

											// find out if it hits the infinite plate:
			float t = (0. - ycir) / outy;

			// If t is less than 0., then the reflected beam went up instead of down. no hit.
			if (t < 0) prods[tnum] = 0;

			// Otherwise, this beam hit the infinite plate. (Case D) Increment the number of hits.
			else prods[tnum] = 1;
		}
	}
	

	for (int offset = 1; offset < numItems; offset *= 2)
	{
		int mask = 2 * offset - 1;
		__syncthreads();
		if ((tnum & mask) == 0)
		{
			prods[tnum] += prods[tnum + offset];
		}
	}

	__syncthreads();
	if (tnum == 0)
		C[wgNum] = prods[0];
}


// main program:

int
main( int argc, char* argv[ ] )
{
	int dev = findCudaDevice(argc, (const char **)argv);

	/**************************Setting Up the Memory for the Arrays*************************/
	// allocate host memory:

	float * hA = new float [ SIZE ];				// xcs
	float * hB = new float [ SIZE ];				// ycs
	float *hR = new float[SIZE];					// rs
	float * hC = new float [ SIZE/BLOCKSIZE ];

	TimeOfDaySeed();		// seed the random number generator

	for( int i = 0; i < SIZE; i++ )
	{
		/*hA[i] = hB[i] = (float) sqrt(  (float)(i+1)  );*/
		hA[i] = Ranf(XCMIN, XCMAX);
		hB[i] = Ranf(YCMIN, YCMAX);
		hR[i] = Ranf(RMIN, RMAX);
	}

	// allocate device memory:

	float *dA, *dB, *dR, *dC;

	dim3 dimsA( SIZE, 1, 1 );
	dim3 dimsB( SIZE, 1, 1 );
	dim3 dimsR(SIZE, 1, 1);
	dim3 dimsC( SIZE/BLOCKSIZE, 1, 1 );

	//__shared__ float prods[SIZE/BLOCKSIZE];


	cudaError_t status;
	status = cudaMalloc( reinterpret_cast<void **>(&dA), SIZE*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc( reinterpret_cast<void **>(&dB), SIZE*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc(reinterpret_cast<void **>(&dR), SIZE * sizeof(float));
		checkCudaErrors(status);
	status = cudaMalloc( reinterpret_cast<void **>(&dC), (SIZE/BLOCKSIZE)*sizeof(float) );
		checkCudaErrors( status );


	/*************Copying the Arrays from the Host to the Device***********/
	// copy host memory to the device:

	status = cudaMemcpy( dA, hA, SIZE*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );
	status = cudaMemcpy( dB, hB, SIZE*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );
	status = cudaMemcpy(dR, hR, SIZE * sizeof(float), cudaMemcpyHostToDevice);
		checkCudaErrors(status);


	/**********************Getting Ready to Execute***********************/
	// setup the execution parameters:

	dim3 threads(BLOCKSIZE, 1, 1 );
	dim3 grid( SIZE / threads.x, 1, 1 );

	// Create and start timer

	cudaDeviceSynchronize( );

	// allocate CUDA events that we'll use for timing:

	cudaEvent_t start, stop;
	status = cudaEventCreate( &start );
		checkCudaErrors( status );
	status = cudaEventCreate( &stop );
		checkCudaErrors( status );

	// record the start event:

	status = cudaEventRecord( start, NULL );
		checkCudaErrors( status );


	/*****************************Executing the Kernel******************************/
	// execute the kernel:

	for( int t = 0; t < NUMTRIALS; t++)
	{
		Monte_Carlo <<< grid, threads >>>( dA, dB, dR, dC );
	}


	/****************************Getting the Stop Time*******************************/
	// record the stop event:

	status = cudaEventRecord( stop, NULL );
		checkCudaErrors( status );

	// wait for the stop event to complete:

	status = cudaEventSynchronize( stop );
		checkCudaErrors( status );


	/******************************Printing the Performance*****************************/
	float msecTotal = 0.0f;
	status = cudaEventElapsedTime( &msecTotal, start, stop );
		checkCudaErrors( status );

	// compute and print the performance

	double secondsTotal = 0.001 * (double)msecTotal;
	double multsPerSecond = (float)SIZE * (float)NUMTRIALS / secondsTotal;
	double gigaMultsPerSecond = multsPerSecond / 1000000000.;
	fprintf( stderr, "\t\t\t\t\tgigaTrialsPerSecond/Second = %10.2lf\t\tArray Size = %10d\tBLOCKSIZE = %d\n", SIZE, BLOCKSIZE, gigaMultsPerSecond );


	/***************************Copying the Array from the Device to the Host*****************************/
	// copy result from the device to the host:

	status = cudaMemcpy( hC, dC, (SIZE/BLOCKSIZE)*sizeof(float), cudaMemcpyDeviceToHost );
		checkCudaErrors( status );

	// check the sum :

	double sum = 0.;
	for(int i = 0; i < SIZE/BLOCKSIZE; i++ )
	{
		//fprintf(stderr, "hC[%6d] = %10.2f\n", i, hC[i]);
		sum += (double)hC[i];
	}
	fprintf( stderr, "sum = %10.2lf, probability = %10.2lf\n", sum, sum / double(SIZE) );


	/************************Cleaning Up******************************/
	// clean up memory:
	delete [ ] hA;
	delete [ ] hB;
	delete[] hR;
	delete [ ] hC;

	status = cudaFree( dA );
		checkCudaErrors( status );
	status = cudaFree( dB );
		checkCudaErrors( status );
	status = cudaFree(dR);
		checkCudaErrors(status);
	status = cudaFree( dC );
		checkCudaErrors( status );


	return 0;
}