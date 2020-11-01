//
// Created by Chris Eberle on 29.10.20.
//

#ifndef MONTECARLO_RANDOM_H
#include <math.h>
#include <assert.h>
#include <stdint.h>
#define MONTECARLO_RANDOM_H
#define LOW 0.02425
#define HIGH 0.97575

typedef struct XORRNGState32 { uint32_t state; } xorRngState32;
typedef struct XORRNGState64 { uint64_t state; } xorRngState64;
typedef struct XORRNGState128 { uint32_t state_a, state_b, state_c, state_d; } xorRngState128;
typedef struct TransitionProbability { char transition[100]; double p; } transitionProbability;

transitionProbability* init_pvector(int n_transitions, char transitionNames[n_transitions][100], double probabilities[n_transitions]);
double get_transition_probability(int pvector_size, transitionProbability* pvector, char transitionName[100]);

uint32_t xorshift32(xorRngState32* rng);
uint64_t xorshift64(xorRngState64* rng);
uint32_t xorshift128(xorRngState128* rng);
void seedXorRngState128(xorRngState128* rng);

void fillArrayXor64Int(int size, int array[size]);
void fillArrayXor64Unit(int size, double array[size]);
void fillArrayNormal(double mu, double sigma, int size, double array[size]);

extern long getRandomSeed(void);
double toUnit32(uint32_t x);
double toUnit64(uint64_t x);

double normInvCDF(double p);
double normal(double mu, double sigma);

#endif //MONTECARLO_RANDOM_H
