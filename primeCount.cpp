#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <iostream>
#include <mpi.h>
//using namespace std;

int CountPrimes(unsigned long long int n0, unsigned long long int n1, unsigned long int *basePrimes, int nBasePrimes )
{
	// Use the sieve of Eratosthenes to count the number of primes between n0 and n1-1 inclusive
	// Assumptions:
	// * basePrimes is a pointer to a list of known primes in ascending order
	// n0 > max(basePrimes)
	// * n1 <= max(basePrimes)^2 (a composite number must have a prime factor <= its square root)
	// * n1-n0 is small enough so we can allocate an array of (n1-n0)*64 bits.
	int nPrimes = 0;
	int i;
	unsigned long long int *working; 
	working = (unsigned long long int*)malloc( (n1-n0)*sizeof(unsigned long long int));
	memset( working, 0, (n1-n0)*sizeof(unsigned long long int) );

	for ( i = 0; i < nBasePrimes; i++ )
	{
		// find starting point
		long long int i0;
		if ( (n0 % basePrimes[i]) == 0 )
		{
			i0 = n0;
			//printf("i0 = %d\n", i0);
		}
		else
		{
			i0 = (n0 / basePrimes[i] + 1)*basePrimes[i];
		}
		while ( i0 < n1 )
		{
			working[i0 - n0] = 1;
			i0 += basePrimes[i];
		}
	}
	for ( i = 0; i < (n1-n0); i++ )
		if ( !working[i] )
			nPrimes++;
	free( working );
	return nPrimes;
}

int main(int argc, char **argv)
{
	#define MAX_PRIME 100000 // we are counting primes up to 1e10, so we need a base primes list up to sqrt(1e10) = 1e5
	#define SQRT_MAX_PRIME 317 // the largest number we will need to cross out multiples of for our base prime list
	unsigned long long int nPrimes = 0;
	unsigned long long int *working = (unsigned long long int*)malloc( (MAX_PRIME + 1)*sizeof(unsigned long  int) );
	unsigned long int *primes = (unsigned long int *)malloc( MAX_PRIME*sizeof(unsigned long int) );
	
	//partial buffer used by the ranks for computation
	unsigned long long int *subIp = (unsigned long long int*)malloc(MAX_PRIME*sizeof(unsigned long long int)); 
	
	// construct the base prime list (primes up to 1e5) using the sieve of Eratosthenes
	int i, j;
	int mpiSize, mpiRank;
	unsigned long long int ip;
	unsigned long long int delta = 100000000/10;//00; delta ca be divided up to 1000
	unsigned long long int pmax = 10000000000;
	unsigned long long int countIp = 1;
	double startTime, endtime;
	
	//open the file
        FILE *myFile;
        myFile = fopen("serialMeasurements.txt", "a");
	
	//Initialise the mpi
	MPI_Init(&argc, &argv);
	
	//start timing...
	startTime = MPI_Wtime();	

	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	
	if(mpiRank == 0)
	{
		for ( i = 0; i < MAX_PRIME; i++ )
		{
			working[i] = 1;
		}
		working[0] = working[1] = 0;

		// sieve
		for ( i = 2; i < SQRT_MAX_PRIME; i++ )
		{
			if ( working[i] )
			{
				for ( j = 2*i; j < MAX_PRIME; j += i )
				{
					working[j] = 0;
				}
			}
		}
		nPrimes = 0;

		for ( i = 0; i < MAX_PRIME; i++ )
			if ( working[i] )
				primes[nPrimes++] = i;
		free( working );
	}

	//Broadcast the elements that needs to be passed to the rank for the computation
	MPI_Bcast(&nPrimes,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(primes,nPrimes,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);	
	
	ip = primes[nPrimes-1] + 1;
	
	//setting the value in all ranks
       	subIp[0] = ip + delta * mpiRank;        
                   
	//adding the delta to the values already in the ranks delta * mpiSize
  	int x = 1;
		
        while(subIp[x-1] + delta * mpiSize < pmax)
        {	
           	subIp[x] = subIp[x-1] + delta * mpiSize;
                countIp++;
                x++;
        }
	
	unsigned long long int nBase = nPrimes;
	unsigned long long int sumPrimes = 0;
	
	for(int i = 0; i < countIp; i++)
	{
		unsigned long long int next = subIp[i] + delta;

		if ( next > pmax )
			next = pmax;
		printf("counting from %lld to %lld\n",subIp[i],next);
		sumPrimes += CountPrimes( subIp[i], next, primes, nBase );
        }
	     
	//create a variable to hold the ranks results
	unsigned long long int totPrimes = 0;

	//reducing process that will sum all the local sums
	MPI_Reduce(&sumPrimes,&totPrimes,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	
	free(primes);
	free(subIp);
	
	if(mpiRank == 0)
	{
		// add the first part of the sieve to the total
		printf("prime counting function(%lld) = %lld\n",pmax,totPrimes + nPrimes);
		
		//end time....
		endtime = MPI_Wtime();
		
		printf("The time to find the prime numbers took: %f seconds\n",endtime - startTime);
		
		//print to .txt file for assignment report
		fprintf(myFile,"The time to complete the computation of the primes is: %f", endtime - startTime);	
		
		//close the file
		fclose(myFile);	
	}
		
	//Finalise the mpi
	MPI_Finalize();
	return 0;
}

