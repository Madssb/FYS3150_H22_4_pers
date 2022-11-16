/*
For now, a test program to implement parallel coding using OpenMP
*/
/*
Main code for simulating the Ising model
*/

#include "omp.h"        // OpenMP header
#include "lattice.hpp"
#include "mcmc.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>

using namespace std;
using namespace arma;

int main(int argc, const char* argv[])
{
  const int threads = atoi(argv[2]);          // No. of threads
  const int L = atoi(argv[3]);                // Lattice length
  const int nCyclesPerThread = atoi(argv[4]); // No. of MCMC cycles per thread
  const double T = atof(argv[5]);             // Temperature [J/k]
  const string ordered_str = argv[6];         // Weather the lattice is ordered or not
  const string outputFilename = argv[7];      // Filename for output
  const string simulationType = argv[8];      // Phase transition or not
  const int N = L * L;                        // No. of spin particles

  // If ordered, all spins in the lattice are set to +1 (spin up)
  bool ordered;

  if (ordered_str == "ordered")
  {
    ordered = true;
  }

  else
  {
    ordered = false;
  }

  double E = 0;
  double M = 0;

  mat lattice = initialize_lattice(L, T, E, M, ordered);

  // Storing values to matrix
  cube results;

  // Base seed to be used either as base for independent
  // seeds per thread or as seed for serial computing
  unsigned int baseSeed = std::chrono::system_clock::now().time_since_epoch().count();

  if (simulationType == "phasetransition")
  {
    vec tempList = linspace(2.1, 2.4, 6);
    vec L_list = vec{40, 60, 80, 100};

    results = cube(6, 2, 4);
    double time;

    #ifdef _OPENMP
    {
      #pragma omp parallel
      {
        const int nThreads = omp_get_num_threads();
        const int threadID = omp_get_thread_num();
        const int weight = 1. / nThreads;

        mt19937 generator;
        unsigned int threadSeed = baseSeed + threadID;
        generator.seed(threadSeed);

        mat threadLattice;
        int L_phase;
        int N_phase;
        double T_phase;
        double threadE, threadM;
        double sumE, sumEE, sumM, sumMM;
        double avgEps, avgEps_sqrd, avgM, avgM_sqrd;
        double heatCap, X;

        for (int l = 0; l < 4; l++)
        {
          L_phase = L_list(l);
          N_phase = L_phase * L_phase;

          auto t1 = chrono::high_resolution_clock::now();

          for (int t = 0; t < 6; t++)
          {
            T_phase = tempList(t);
            threadE = E;
            threadM = M;
            sumE = 0;
            sumM = 0;
            sumMM = 0;
            threadLattice = initialize_lattice(L_phase, T_phase, threadE, threadM, false);

            for (int n = 0; n < nCyclesPerThread; n++)
            {
              mcmc(threadLattice, generator, L_phase, threadE, threadM, T_phase);

              #pragma omp critical
              {
                sumE += weight * threadE;
                sumEE += weight * threadE * threadE;
                sumM += weight * fabs(threadM);
                sumMM += weight * threadM * threadM;
              }

              avgEps = sumE / N_phase;
              avgEps_sqrd = sumEE / (N_phase * N_phase);
              avgM = sumM / N_phase;
              avgM_sqrd = sumMM / (N_phase * N_phase);

              heatCap = N_phase * (avgEps_sqrd - avgEps * avgEps) / (T_phase * T_phase);
              X = N_phase * (avgM_sqrd - avgM * avgM) / T_phase;

              results(t, 0, l) = heatCap / nCyclesPerThread;
              results(t, 1, l) = X / nCyclesPerThread;
            }
          }
          }

          auto t2 = chrono::high_resolution_clock::now();
          double duration = chrono::duration<double>(t2 - t1).count();

          #pragma omp critical
          {
            time = weight * duration;

            cout << "Completed simulation for L=" << L_phase << "x" << L_phase
                 << " after " << int(time) / 60
                 << " min, "  << int(time) % 60
                 << " sec"
                 << endl;
          }
        }
      }

      results.save(outputFilename, raw_ascii);
    }
   }


  #else
  {
    results = cube(nCyclesPerThread, 4, threads);

    // If the -fopenmp flag is passed to the compiler, the code
    // will run in parallel
    #ifdef _OPENMP
    {
      // Start parallell region
      #pragma omp parallel
      {
        const int threadID = omp_get_thread_num();

        mt19937 generator;
        unsigned int threadSeed = baseSeed + threadID;
        generator.seed(threadSeed);

        // Each thread gets their own copy of the lattice to work with
        // as well as variables to update
        mat threadLattice = lattice;

        double threadE = E;
        double threadM = M;
        double sumE = 0;
        double sumEE = 0;
        double sumM = 0;
        double sumMM = 0;
        double avgEps, avgEps_sqrd, avgM, avgM_sqrd, heatCap, X;

        // Parallelized for loop. Note that only the first loop is
        // parallelized.
        #pragma omp for
        for (int i = 0; i < threads; i++)
        {
          for (int n = 1; n <= nCyclesPerThread; n++)
          {
            mcmc(threadLattice, generator, L, threadE, threadM, T);

            sumE += threadE;
            sumEE += (threadE * threadE);
            sumM += fabs(threadM);
            sumMM += (threadM * threadM);

            avgEps = sumE / (n * N);
            avgEps_sqrd = sumEE / (n * N * N);
            avgM = sumM / (n * N);
            avgM_sqrd = sumMM / (n * N * N);
            heatCap = N * (avgEps_sqrd - avgEps * avgEps) / (T * T);
            X = N * (avgM_sqrd - avgM * avgM) / T;

            // Saving results to individual slices belonging to each thread
            results.slice(i).row(n - 1) = rowvec{avgEps, avgM, heatCap, X};
          }
        }
      }
    }
  }


    // If the flag -fopenmp is not passed during compiling, the
    // program will run in series
    #else
    {
      mt19937 generator;
      generator.seed(baseSeed);

      double serialE = E;
      double serialM = M;
      double sumE = 0;
      double sumEE = 0;
      double sumM = 0;
      double sumMM = 0;
      double avgEps, avgEps_sqrd, avgM, avgM_sqrd, heatCap, X;

      for (int n = 1; n <= nCyclesPerThread; n++)
      {
        mcmc(lattice, generator, L, serialE, serialM, T);

        sumE += serialE;
        sumEE += (serialE * serialE);
        sumM += fabs(serialM);
        sumMM += (serialM * serialM);

        avgEps = sumE / (n * N);
        avgEps_sqrd = sumEE / (n * N * N);
        avgM = sumM / (n * N);
        avgM_sqrd = sumMM / (n * N * N);
        heatCap = N * (avgEps_sqrd - avgEps * avgEps) / (T * T);
        X = N * (avgM_sqrd - avgM * avgM) / T;

        results.slice(0).row(n - 1) = rowvec{avgEps, avgM, heatCap, X};
      }
    }

    #endif
    }
    #endif
    }

    // Saving the results to matrix in raw ascii (mostly to avoid the
    // header of the matrix)
    results.save(outputFilename, raw_ascii);



  return 0;
}
