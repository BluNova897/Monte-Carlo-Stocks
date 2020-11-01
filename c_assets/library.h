#ifndef MONTECARLO_LIBRARY_H
#define MONTECARLO_LIBRARY_H

void cumsum_d(int size, double arr[size]);
void GBM(int steps, double max_t, double W[steps], double S0, double mu, double sigma, double out_array[steps]);
void _gmbThreadedInit(int n_sims, int n_steps, int n_threads, double S0, double W0, double mu, double sigma, double ret_array[n_sims*n_steps]);
void* _gbmThreaded(void *fargs);
#endif //MONTECARLO_LIBRARY_H
