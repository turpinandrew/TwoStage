#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dist.h"

extern double myRandom(void);

/*
** Return the dog value at dist * distanceScale
*/
double 
dog(DogParams dp, double dist, double distanceScale) {

    double distToCell = dist / distanceScale ;
    distToCell *= dist / distanceScale ;

    double c2 = distToCell/dp.cWidth/dp.cWidth;
    double s2 = distToCell/dp.sWidth/dp.sWidth;
    return dp.cHeight * exp(-c2) - dp.sHeight * exp(-s2);
}//dog()

// ***************************************************************
// Poisson stuff
// ***************************************************************
double
fact(int x) {
    double r = 1.0;
    for(double i = 2 ; i <= x ; i++)
        r *= i;
    return r;
}//fact()

// return poisson prob for x 
double 
poisson(int x, double lambda) {
    return exp(-lambda) * pow(lambda, x) / fact(x);
}//poisson()

/*
** Return a lookup table which contains a mapping from uniform random
** on the range 0..POISSON_TABLE_LENGTH to a number of events
** drawn from the poisson distribution.
*/
PoissonTable * 
poissonInitialise(double lambda) {

    PoissonTable *pTable = (PoissonTable *)malloc(sizeof(PoissonTable));
    pTable->table = (int *)malloc(sizeof(int)*POISSON_TABLE_LENGTH);
    if (pTable->table == NULL) {
        fprintf(stderr, "Out of memory for poisson table lambda=%10.6f\n",lambda);
        exit(-1);
    }

    int bin = 0;
    for(int i = 0; i < POISSON_TABLE_LENGTH ; ) {
        int extent = (int)round(poisson(bin, lambda) * (double)POISSON_TABLE_LENGTH);
        if (extent < 1)
            extent = 1;
        int end = i + extent;
        if (end > POISSON_TABLE_LENGTH)
            end = POISSON_TABLE_LENGTH;
        for( ; i < end ; i++ ) {
            pTable->table[i] = bin;
        }
        bin++;
    }

/*
int A[10];
for(int i = 0 ; i < 10 ; i++)
    A[i] = 0;
for(int i = 0 ; i < POISSON_TABLE_LENGTH ; i++)
    A[table[i]]++;
for(int i = 0 ; i < 10 ; i++)
    printf("%8d %8d\n",i,A[i]);
*/

    for(int i = 0 ; i < POISSON_TABLE_LENGTH ; i++) {
        if (pTable->table[i] == 0)
            pTable->lastZeroIndex = i;
        if (pTable->table[i] == 1)
            pTable->lastOneIndex = i;
    }
        
    return pTable;
}//poissonInitialise()

/*
** Return a random number from poisson distribution "table"
** as created by poissonInitialise()
*/
int
poissonRand(PoissonTable *pTable) { 
    int r = (int)round(myRandom() * (double)POISSON_TABLE_LENGTH / (double)RAND_MAX); 
    if (r == POISSON_TABLE_LENGTH) 
        r = POISSON_TABLE_LENGTH - 1; 
    if (r <= pTable->lastZeroIndex) return 0;
    if (r <= pTable->lastOneIndex) return 1;
    return pTable->table[r];
}//poissonRand()

/*
** Return the max possible value that can be returned by a random
** sample from table[].
*/
int
poissonMax(PoissonTable *pTable) {
    return pTable->table[POISSON_TABLE_LENGTH - 1];
}
