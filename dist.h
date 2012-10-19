#include <math.h>

typedef struct dogParams {
    double cHeight, cWidth, sHeight, sWidth;
} DogParams;

double dog(DogParams dp, double dist, double distanceScale);

typedef struct poissonTable {
    int *table;
    int lastZeroIndex;
    int lastOneIndex;
} PoissonTable;

#define POISSON_TABLE_LENGTH (1 << 16)
PoissonTable *poissonInitialise(double lambda);
int poissonRand(PoissonTable *poissonTable);
int poissonMax(PoissonTable *table);

/*
#define POISSON_RAND(_table, _v)                                                    \
do {                                                                                \
    int r = (int)round(rand() * (double)POISSON_TABLE_LENGTH / (double)RAND_MAX);   \
    if (r == POISSON_TABLE_LENGTH)                                                  \
        r = POISSON_TABLE_LENGTH - 1;                                               \
    _v = ((_table)->table)[r];                                                               \
} while (0);
*/
#define POISSON_RAND(_table, _v) _v = poissonRand(_table);
