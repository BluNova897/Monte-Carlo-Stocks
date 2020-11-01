#include "library.h"
#include <math.h>
#include "random.h"
#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

double* times;
double* gmbTrajectories;

double square(double x) {
    return x*x;
}

double* linspace(double start, double end, int steps) {
    assert(end > start);
    assert(start >= 0 && end > 0 && steps > 0);
    double* out = (double*)malloc(steps*sizeof(int));
    double delta = (end-start)/(steps-1);
    for (int i = 0; i < steps; ++i) {
        out[i] = (start + i*delta);
    }
    return out;
}

void cumsum_d2(int size, double arr[size]) {
    for (int i = 1; i < size+1; ++i) {
        arr[i] += arr[i-1];
    }
}

void cumsum_d(int size, double arr[size]) {
    double* out = (double*)malloc(size*sizeof(double));
    for (int i = 0; i < size; ++i) {
        double temp = 0;
        for (int j = 0; j < i+1; ++j) {
            temp += arr[j];
        }
        out[i] = temp;
    }
    for (int i = 0; i < size; ++i) {
        arr[i] = out[i];
    }
    free(out);}

void calcBrownianIncrements(int steps, double W0, double * W, double * b) {
    double dt = 1/((double) steps);
    fillArrayNormal(0, sqrt(dt), steps, b);
    W[0] = W0;
    for (int i = 1; i < steps+1; ++i) {
        W[i] = b[i-1];
    }
    cumsum_d2(steps+1, W);
}

void GBM(int steps, double max_t, double W[steps], double S0, double mu, double sigma, double out_array[steps]) {
    double drift, diffusion;
    out_array[0] = S0;
    double* tau = malloc(steps * sizeof(double)); // allocate array containing the time steps t ∈ [0, max_t]
    double delta_tau = max_t/(steps-1);
    for (int i = 0; i < steps; ++i) {
        tau[i] = (double) i*delta_tau;
    }
    for (int i = 1; i < steps; ++i) {
        drift = (mu-square(sigma)/2)*tau[i];
        diffusion = sigma*W[i-1];
        out_array[i] = S0*exp(drift+diffusion);
    }
    free(tau);
}

struct gmbArgType {
    double W0;
    double S0;
    double mu;
    double sigma;
    int steps;
    int ptr_start;
    int ptr_end;
    int n_ptr;
    int thread_id;
};

void printArr(int size, double* arr) {
    printf("[%f", arr[0]);
    for (int i = 1; i < size; ++i) {
        printf(",%f", arr[i]);
    }
    printf("]\n");
}

void* _gbmThreaded(void *fargs) {
    errno = 0;
    struct gmbArgType* args = (struct gmbArgType*) fargs;
    int ptr_length = (args->ptr_end - args->ptr_start) / args->n_ptr;
    for (int i = 0; i < args->n_ptr; ++i) {
        // simulation of geometric brownian process start here
        double drift, diffusion;
        // calculate brownian increments first
        double dt = 1/((double) args->steps);
        double* b = malloc(args->steps*sizeof(double));
        double* W = malloc((args->steps+1)*sizeof(double));
        W[0] = args->W0;
        fillArrayNormal(0, sqrt(dt), args->steps, b);
        W[0] = args->W0;
        for (int j = 1; j < args->steps+1; ++j) {
            W[j] = b[j-1];
        }
        cumsum_d2(args->steps+1, W);
        for (int j = i*ptr_length; j < (i+1)*ptr_length; ++j) {
            drift = (args->mu-square(args->sigma)/2)*times[j%ptr_length];
            diffusion = args->sigma*W[(j%ptr_length)];
            double S_t = args->S0*exp(drift+diffusion);
            if (isinf(S_t)) {
                printf("Error at %d from thread %d\n", args->ptr_start+j, args->thread_id);
                printf("b = ");
                printArr(args->steps, b);
                printf("W = ");
                printArr(args->steps, W);
                printf("accessing W @ index %d", (j%ptr_length));
                printf("diffusion = %f", diffusion);
                printf("drift = %f", drift);
                printf("%s\n", strerror(errno));
                assert(0);
            }
            assert(S_t != 0);
            gmbTrajectories[args->ptr_start + j] = S_t;
        }
        free(b);
        free(W);
    }
    free(fargs);
    pthread_exit(NULL);
}

void _gmbThreadedInit(int n_sims, int n_steps, int n_threads, double S0, double W0, double mu, double sigma, double ret_array[n_sims*n_steps]) {
    times = malloc(n_steps * sizeof(double)); // allocate array of time steps
    for (int i = 0; i < n_steps; ++i) {
        times[i] = i/((double) n_steps-1);
    }
    gmbTrajectories = calloc(n_sims*n_steps, sizeof(double)); // allocate array containing n_sims simulations of length n_steps
    pthread_t threads[n_threads];
    int sims_per_thread = n_sims/n_threads; // # of simulations each thread will carry out
    int sims_remainder = n_sims%n_threads; // remainder of simulations will be carried out by last thread
    for (int i = 0; i < n_threads; ++i) {
        struct gmbArgType* args = malloc(sizeof(struct gmbArgType));
        args->mu = mu;
        args->sigma = sigma;
        args->S0 = S0;
        args->thread_id=i;
        if (i == n_threads-1) {
            args->ptr_start = i*sims_per_thread*n_steps;
            args->ptr_end = n_sims*n_steps;
            args->n_ptr = sims_per_thread+sims_remainder;
            args->steps = n_steps;
            args->W0 = W0;
        }
        else {
            args->ptr_start = i*sims_per_thread*n_steps;
            args->ptr_end = (i+1)*sims_per_thread*n_steps;
            args->n_ptr = sims_per_thread;
            args->steps = n_steps;
            args->W0 = W0;
        }
        pthread_create(&threads[i], NULL, _gbmThreaded, (void*) args);
    }

    for (int i = 0; i < n_threads; ++i) {
        pthread_join(threads[i], NULL);
    }
    for (int i = 0; i < n_sims*n_steps; ++i) {
        ret_array[i] = gmbTrajectories[i];
    }
    free(gmbTrajectories);
    free(times);
}