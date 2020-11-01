//
// Created by Chris Eberle on 29.10.20.
//

#include "random.h"
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#define LOW 0.02425
#define HIGH 0.97575

static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00};
static const double b[] = {-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01};
static const double c[] = {-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00};
static const double d[] = {7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00};


long getRandomSeed(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    long seed = 0;
    // ensure non-zero seed to avoid xorshift fixed point
    while (seed == 0) {
        long sec = (long) ts.tv_sec;
        long nsec = (long) ts.tv_nsec;
        int shift = nsec % 9;
        seed = (sec ^ (nsec << shift));
    }
    return seed;
}

void fillArrayXor64Int(int size, int array[size]) {
    xorRngState64* rng = (xorRngState64*) malloc(sizeof(xorRngState64));
    rng->state = getRandomSeed();
    for (int i = 0; i < size; ++i) {
        array[i] = xorshift64(rng);
    }
    free(rng);
}

void fillArrayXor64Unit(int size, double array[size]) {
    xorRngState64* rng = (xorRngState64*) malloc(sizeof(xorRngState64));
    rng->state = getRandomSeed();
    for (int i = 0; i < size; ++i) {
        array[i] = toUnit64(xorshift64(rng));
    }
    free(rng);
}

void seedXorRngState128(xorRngState128* rng) {
    rng->state_a = getRandomSeed();
    rng->state_b = getRandomSeed();
    rng->state_c = getRandomSeed();
    rng->state_d = getRandomSeed();
}

uint32_t xorshift32(xorRngState32* rng) {
    uint32_t x = rng->state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng->state = x;
    return x;
}

uint64_t xorshift64(xorRngState64* rng) {
    uint64_t x = rng->state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng->state = x;
    return x;
}

uint32_t xorshift128(xorRngState128* rng) {
    uint32_t t = rng->state_d;
    uint32_t const s = rng->state_a;
    rng->state_d = rng->state_c;
    rng->state_c = rng->state_b;
    rng->state_b = s; // s = rng->state_a;
    t ^= t << 11;
    t ^= t >> 8;
    rng->state_a = t ^ s ^ (s >> 19);
    return rng->state_a;
}

double toUnit64(uint64_t x) {return (double)x/UINT64_MAX;};

double toUnit32(uint32_t x) {return (double)x/UINT32_MAX;};

double normal(double mu, double sigma) {
    xorRngState64* rng = (xorRngState64*)malloc(sizeof(xorRngState64));
    rng->state = getRandomSeed();
    double x = mu + sigma*normInvCDF(toUnit64(xorshift64(rng)));
    free(rng);
    return x;
}

void fillArrayNormal(double mu, double sigma, int size, double array[size]) {
    xorRngState64* rng = (xorRngState64*) malloc(sizeof(xorRngState64));
    rng->state = getRandomSeed();
    for (int i = 0; i < size; ++i) {
        array[i] = mu + sigma*normInvCDF(toUnit64(xorshift64(rng)));
    }
    free(rng);
}

double normInvCDF(double p) {
    double q, r;

    errno = 0;

    if (p < 0 || p > 1)
    {
        errno = EDOM;
        return 0.0;
    }
    else if (p == 0)
    {
        errno = ERANGE;
        return -HUGE_VAL /* minus "infinity" */;
    }
    else if (p == 1)
    {
        errno = ERANGE;
        return HUGE_VAL /* "infinity" */;
    }
    else if (p < LOW)
    {
        /* Rational approximation for lower region */
        q = sqrt(-2*log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    else if (p > HIGH)
    {
        /* Rational approximation for upper region */
        q  = sqrt(-2*log(1-p));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    else
    {
        /* Rational approximation for central region */
        q = p - 0.5;
        r = q*q;
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }
}

transitionProbability* init_pvector(int n_transitions, char transitionNames[n_transitions][100], double probabilities[n_transitions]) {
    transitionProbability* pvector = malloc(n_transitions*sizeof(transitionProbability));
    for (int i = 0; i < n_transitions; ++i) {
        strcpy(pvector[i].transition, transitionNames[i]);
        pvector[i].p = probabilities[i];
    }
    return pvector;
}

double get_transition_probability(int pvector_size, transitionProbability* pvector, char transitionName[100]) {
    for (int i = 0; i < pvector_size; ++i) {
        if (strcmp(pvector[i].transition, transitionName) == 0) {
            return pvector[i].p;
        }
    }
    return 0;
}