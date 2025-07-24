#include <time.h>
#include <stdlib.h>

float rand0(long *idum);
int randint(long *idum, int range);
float rand2(long *idum);
unsigned long randqd1(unsigned long *idum);
int random1(long *seed);

/* randint:                                             */
/* POST: FCTVAL == random integer from 0 to (range - 1) */
int randint(long *idum, int range)
{
    return (int)((float)range * rand0(idum));
}

/* Thanks to "Numerical Recipes in C", 2nd Ed., Cambridge Univ. Press */
/* for random number generators (rand0, rand2, and randqd1).          */

#define IA 16807
#define IM0 2147483647
#define AM0 (1.0 / IM0)
#define IQ 127773
#define IR 2836
#define RAND_MASK 123459876

/* rand0:                               */
/* FCTVAL == random float number        */
float rand0(long *idum)
{
    long k;
    float ans;

    *idum ^= RAND_MASK;
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0)
        *idum += IM0;
    ans = AM0 * (*idum);
    *idum ^= RAND_MASK;
    return ans;
}

#define IM1 2147483563
#define IM2 2147483399
#define AM2 (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

/* rand2:                               */
/* FCTVAL == random float number        */
/* Execution time: 2.0 * rand0 (approx) */
float rand2(long *idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--)
        {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0)
        *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp = AM2 * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

/* randqd1:                             */
/* FCTVAL == random integer number      */
/* Execution time: 0.1 * rand0 (approx) */
unsigned long randqd1(unsigned long *idum)
{
    *idum ^= RAND_MASK;
    *idum = 1664525L * (*idum) + 1013904223L;
    *idum ^= RAND_MASK;
    *idum >>= 1;
    return *idum;
}

/* random:                                     */
/* using rand function from standard C library */
int random1(long *seed)
{
    if (*seed == 0)
    {
        time_t t;
        t = (long)time((long *)seed);
        *seed = (int)t;
    }
    srand(*seed);
    *seed = rand();
    return (int)*seed;
}
