/**
 * Runs a simulation of the n-body problem in 3D.
 * 
 * To compile the program:
 *   gcc -Wall -fopenmp -O3 -march=native nbody-p3.c matrix.c util.c -o build/nbody-p3 -lm
 * 
 * To run the program:
 *   ./nbody-p3 time-step total-time outputs-per-body input.npy output.npy [opt: num-threads]
 * where:
 *   - time-step is the amount of time between steps (Δt, in seconds)
 *   - total-time is the total amount of time to simulate (in seconds)
 *   - outputs-per-body is the number of positions to output per body
 *   - input.npy is the file describing the initial state of the system (below)
 *   - output.npy is the output of the program (see below)
 *   - last argument is an optional number of threads (a reasonable default is
 *     chosen if not provided)
 * 
 * input.npy has a n-by-7 matrix with one row per body and the columns:
 *   - mass (in kg)
 *   - initial x, y, z position (in m)
 *   - initial x, y, z velocity (in m/s)
 * 
 * output.npy is generated and has a (outputs-per-body)-by-(3n) matrix with each
 * row containing the x, y, and z positions of each of the n bodies after a
 * given timestep.
 * 
 * See the PDF for implementation details and other requirements.
 * 
 * AUTHORS: David Olsakowski,
 */

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <omp.h>

#include "matrix.h"
#include "util.h"

// Gravitational Constant in N m^2 / kg^2 or m^3 / kg / s^2
#define G 6.6743015e-11

// Softening factor to reduce divide-by-near-zero effects
#define SOFTENING 1e-9

void universal_gravitation(double* x, double* y, double* z, double* vx, double* vy, double* vz, double* masses, size_t n, double time_step,size_t num_threads){
  //printf("%d",omp_get_num_threads());

 // omp_set_num_threads(num_threads);
  //#pragma omp parallel for default(none) shared(x,y,z,vx,vy,vz,n,masses,time_step) 

    for(size_t i=0;i<n;i++){
        // go to i in s3
   // #pragma omp parallel for default(none) shared(i,x,y,z,vx,vy,vz,n,masses,time_step) private(j)

        for(size_t j=0;j<n;j++){
            if (i == j) { continue; }
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            double dist = sqrt(dx*dx + dy*dy + dz*dz + SOFTENING);
            double F = G * time_step / (dist*dist*dist);
            double Fi = F * masses[j];
            vx[i] -= Fi * dx;
            vy[i] -= Fi * dy;
            vz[i] -= Fi * dz;
        }
    }
}

void save_outputs(Matrix* output, double* x, double* y, double* z, size_t n, size_t row){
    for(size_t i=0;i<n;i++){
        output->data[row*3*n + 0 + 3*i] = x[i];
        output->data[row*3*n + 1 + 3*i] = y[i];
        output->data[row*3*n + 2 + 3*i] = z[i];
    }
}
int main(int argc, const char* argv[]) {
    // parse arguments
    if (argc != 6 && argc != 7) { fprintf(stderr, "usage: %s time-step total-time outputs-per-body input.npy output.npy [num-threads]\n", argv[0]); return 1; }
    double time_step = atof(argv[1]), total_time = atof(argv[2]);
    if (time_step <= 0 || total_time <= 0 || time_step > total_time) { fprintf(stderr, "time-step and total-time must be positive with total-time > time-step\n"); return 1; }
    size_t num_outputs = atoi(argv[3]);
    if (num_outputs <= 0) { fprintf(stderr, "outputs-per-body must be positive\n"); return 1; }
    size_t num_threads = argc == 7 ? atoi(argv[6]) : get_num_cores_affinity()/2; // TODO: you may choose to adjust the default value
    if (num_threads <= 0) { fprintf(stderr, "num-threads must be positive\n"); return 1; }
    Matrix* input = matrix_from_npy_path(argv[4]);
    if (input == NULL) { perror("error reading input"); return 1; }
    if (input->cols != 7) { fprintf(stderr, "input.npy must have 7 columns\n"); return 1; }
    size_t n = input->rows;
    if (n == 0) { fprintf(stderr, "input.npy must have at least 1 row\n"); return 1; }
    if (num_threads > n) { num_threads = n; }
    size_t num_steps = (size_t)(total_time / time_step + 0.5);
    if (num_steps < num_outputs) { num_outputs = 1; }
    size_t output_steps = num_steps/num_outputs;
    num_outputs = (num_steps+output_steps-1)/output_steps;

    // variables available now:
    //   time_step    number of seconds between each time point
    //   total_time   total number of seconds in the simulation
    //   num_steps    number of time steps to simulate (more useful than total_time)
    //   num_outputs  number of times the position will be output for all bodies
    //   output_steps number of steps between each output of the position
    //   input        n-by-7 Matrix of input data
    //   n            number of bodies to simulate

    // start the clock
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    double* masses = malloc(n * sizeof(double));
    double* x = malloc((n*3) * sizeof(double));
    double* y = malloc((n*3) * sizeof(double));
    double* z = malloc((n*3) * sizeof(double));
    double* vx = malloc((n*3) * sizeof(double));
    double* vy = malloc((n*3) * sizeof(double));
    double* vz = malloc((n*3) * sizeof(double));


    Matrix* next_matrix = matrix_create_raw(n*3, 7);

    // malloc the masses and velocities
    for (size_t i = 0; i < n; i++) {
        masses[i] = input->data[7*i];
        x[i] = input->data[7*i+1];
        y[i] = input->data[7*i+2];
        z[i] = input->data[7*i+3];
        vx[i] = input->data[7*i+4];
        vy[i] = input->data[7*i+5];
        vz[i] = input->data[7*i+6];
    }

    // allocate output matrix as num_outputs x
    Matrix* output = matrix_create_raw(num_outputs, 3*n);

    // save positions to row `0` of output
    save_outputs(output, x, y, z, n, 0);
    
    //printf("%li\n",num_threads);
    // Run simulation for each time step 
     omp_set_num_threads(num_threads);
  #pragma omp parallel for default(none) schedule(static) shared(x,y,z,vx,vy,vz,n,masses,time_step,num_steps,num_threads,output_steps,output) 

    for (size_t t = 1; t < num_steps; t++) { 
        // update the velocities
        universal_gravitation(x, y, z, vx, vy, vz, masses, n, time_step,num_threads);


        // update the positions
        for (size_t i = 0; i < n; i++) {
            x[i] += vx[i] * time_step;
            y[i] += vy[i] * time_step;
            z[i] += vz[i] * time_step;
        }

        // Periodically copy the positions to the output data 
        if (t % output_steps == 0) { 
            // save positions to row `t/output_steps` of output 
            save_outputs(output, x, y, z, n, t/output_steps);
        } 
    } 
    
    // Save the final set of data if necessary 
    if (num_steps % output_steps != 0) { 
        //  save positions to row `num_outputs-1` of output 
        save_outputs(output, x, y, z, n, num_outputs-1);
    }

    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    // save results
    matrix_to_npy_path(argv[5], output);

    // cleanup
    matrix_free(input);
    matrix_free(next_matrix);


    return 0;
}
