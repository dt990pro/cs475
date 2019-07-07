#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>

//The "state" of the system consists of the following global variables :
int		NowYear;		// 2019 - 2024
int		NowMonth;		// 0 - 11

float	NowPrecip;		// inches of rain per month
float	NowTemp;		// temperature this month
float	NowHeight;		// grain height in inches
int		NowNumDeer;		// number of deer in the current population
int		num_behemoth;

int		month_count;

//Your basic time step will be one month.Interesting parameters that you need are :
const float		GRAIN_GROWS_PER_MONTH = 8.0;
const float		ONE_DEER_EATS_PER_MONTH = 0.5;

const float		AVG_PRECIP_PER_MONTH = 6.0;		// average
const float		AMP_PRECIP_PER_MONTH = 6.0;		// plus or minus
const float		RANDOM_PRECIP = 2.0;			// plus or minus noise

const float		AVG_TEMP = 50.0;				// average
const float		AMP_TEMP = 20.0;				// plus or minus
const float		RANDOM_TEMP = 10.0;				// plus or minus noise

const float		MIDTEMP = 40.0;
const float		MIDPRECIP = 10.0;
												//Units of grain growth are inches.
												//Units of temperature are degrees Fahrenheit(°„F).
												//Units of precipitation are inches.

// On visual studio
//omp_lock_t	Lock;
//int			NumInThreadTeam;
//int			NumAtBarrier;
//int			NumGone;
//// specify how many threads will be in the barrier:
////	(also init's the Lock)
//void InitBarrier(int n)
//{
//	NumInThreadTeam = n;
//	NumAtBarrier = 0;
//	omp_init_lock(&Lock);
//}
//
//// have the calling thread wait here until all the other threads catch up:
//void WaitBarrier()
//{
//	omp_set_lock(&Lock);
//	{
//		NumAtBarrier++;
//		if (NumAtBarrier == NumInThreadTeam)
//		{
//			NumGone = 0;
//			NumAtBarrier = 0;
//			// let all other threads get back to what they were doing
//			// before this one unlocks, knowing that they might immediately
//			// call WaitBarrier( ) again:
//			while (NumGone != NumInThreadTeam - 1);
//			omp_unset_lock(&Lock);
//			return;
//		}
//	}
//	omp_unset_lock(&Lock);
//
//	while (NumAtBarrier != 0);	// this waits for the nth thread to arrive
//
//#pragma omp atomic
//	NumGone++;			// this flags how many threads have returned
//}

// useful func
float SQR(float x) {
	return x*x;
}
float Ranf(unsigned int *seedp, float low, float high) {
	float r = (float)rand_r(seedp);              // 0 - RAND_MAX

	return(low + r * (high - low) / (float)RAND_MAX);
}
int Ranf(unsigned int *seedp, int ilow, int ihigh) {
	float low = (float)ilow;
	float high = (float)ihigh + 0.9999f;

	return (int)(Ranf(seedp, low, high));
}
unsigned int seed = 0;  // a thread-private variable


// setting the number of threads:
#ifdef behemoth
#define NUMT		4
#else
#define NUMT		3
#endif

void GrainDeer();
void Grain();
void Watcher();
void Behemoth();

int main(int argc, char *argv[]) {
#ifndef _OPENMP
	fprintf(stderr, "No OpenMP support!\n");
	return 1;
#endif

	// starting date and time:
	NowMonth = 0;
	NowYear = 2019;
	month_count = 1;

	// starting state (feel free to change this if you want):
	NowNumDeer = 1;
	num_behemoth = 0;
	NowHeight = 1.;
	NowTemp = 21.;
	NowPrecip = 7.;
	
	printf("thread = %d\n", NUMT);

	omp_set_num_threads(NUMT);	// same as # of sections

	#pragma omp parallel sections
	{
		#pragma omp section
		{
			GrainDeer();
		}

		#pragma omp section
		{
			Grain();
		}

		#pragma omp section
		{
			Watcher();
		}

		#ifdef behemoth
		#pragma omp section
		{
			Behemoth();	// your own
		}
		#endif
	}       // implied barrier -- all functions must return in order
			// to allow any of them to get past here
}

void GrainDeer() {
	while (NowYear < 2025)
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:
		float	temp_Height = NowHeight;
		int		temp_NumDeer = NowNumDeer;
		int		temp_num_behemoth = num_behemoth;

		/*If the number of graindeer exceeds this value at the end of a month, decrease the number of graindeer by one.
		If the number of graindeer is less than this value at the end of a month, increase the number of graindeer by one.*/
		if (temp_NumDeer > temp_Height) temp_NumDeer--;
		else if (temp_NumDeer < temp_Height) temp_NumDeer++;

		#ifdef behemoth
		if (temp_num_behemoth > 0) temp_NumDeer -= temp_num_behemoth * 3;
		if (temp_NumDeer < 0) temp_NumDeer = 0;
		#endif
		// DoneComputing barrier:
		#pragma omp barrier
		

		NowNumDeer = temp_NumDeer;
		// DoneAssigning barrier:
		#pragma omp barrier
		

		// DonePrinting barrier:
		#pragma omp barrier
			
	}
}

void Grain() {
	while (NowYear < 2025)
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:
		float	temp_Height = NowHeight;
		float	tempFactor = exp(-SQR((NowTemp - MIDTEMP) / 10.));
		float	precipFactor = exp(-SQR((NowPrecip - MIDPRECIP) / 10.));

		temp_Height += tempFactor * precipFactor * GRAIN_GROWS_PER_MONTH;
		temp_Height -= (float)NowNumDeer * ONE_DEER_EATS_PER_MONTH;
		if (temp_Height < 0) temp_Height = 0.;
		// DoneComputing barrier:
		#pragma omp barrier


		NowHeight = temp_Height;
		// DoneAssigning barrier:
		#pragma omp barrier


		// DonePrinting barrier:
		#pragma omp barrier

	}
}

void Watcher() {
	while (NowYear < 2025)
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:


		// DoneComputing barrier:
		#pragma omp barrier


		// DoneAssigning barrier:
		#pragma omp barrier


		// print
		#ifdef behemoth
		printf("%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n", NowYear, NowMonth, month_count, NowPrecip, NowTemp, NowHeight, NowNumDeer, num_behemoth);
		#else
		printf("%d\t%d\t%d\t%f\t%f\t%f\t%d\n", NowYear, NowMonth, month_count, NowPrecip, NowTemp, NowHeight, NowNumDeer);
		#endif

		if (NowMonth == 11) {
			NowMonth = 0;
			NowYear++;
		}
		else NowMonth++;

		month_count++;

		// The temperature and precipitation are a function of the particular month :
		float ang = (30.*(float)NowMonth + 15.) * (M_PI / 180.);

		float temp = AVG_TEMP - AMP_TEMP * cos(ang);
		NowTemp = temp + Ranf(&seed, -RANDOM_TEMP, RANDOM_TEMP);

		float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin(ang);
		NowPrecip = precip + Ranf(&seed, -RANDOM_PRECIP, RANDOM_PRECIP);
		if (NowPrecip < 0.)
			NowPrecip = 0.;
		// DonePrinting barrier:
		#pragma omp barrier

	}
}

void Behemoth() {
	while (NowYear < 2025)
	{
		// compute a temporary next-value for this quantity
		// based on the current state of the simulation:
		int		temp_NumDeer = NowNumDeer;
		int		temp_num_behemoth = num_behemoth;
		int		x = Ranf(&seed, 1, 2);

		if (temp_NumDeer / (1 + temp_num_behemoth) > 3) temp_num_behemoth += x;
		else temp_num_behemoth--;
		if (temp_num_behemoth < 0) temp_num_behemoth = 0;
		
		// DoneComputing barrier:
		#pragma omp barrier


		num_behemoth = temp_num_behemoth;
		// DoneAssigning barrier:
		#pragma omp barrier


		// DonePrinting barrier:
		#pragma omp barrier

	}
}