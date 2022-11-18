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

  // Storing values to matrix
  cube results;

  // Base seed to be used either as base for independent
  // seeds per thread or as seed for serial computing
  unsigned int baseSeed = std::chrono::system_clock::now().time_since_epoch().count();

  mat lattice = initialize_lattice(L, T, E, M, ordered);

  if (simulationType == "phasetransition")
  {
    int numTempElements = 6;
    vec tempList = linspace(2.1, 2.4, numTempElements);

    mat phaseRes = mat(numTempElements, 4);
    double burnIn = nCyclesPerThread / 4.;
    double numSamples = nCyclesPerThread - burnIn;

    #ifdef _OPENMP
    {
      #pragma omp parallel
      {
        const int threadID = omp_get_thread_num();

        mt19937 generator;
        unsigned int threadSeed = baseSeed + threadID;

        generator.seed(threadSeed);

        // Timing algorithm
        auto t1 = chrono::high_resolution_clock::now();

        #pragma omp for
        for (int t = 0; t < numTempElements; t++)
        {

          double T_phase = tempList(t);
          double threadE = 0, threadM = 0;
          double sumE = 0, sumEE = 0, sumM = 0, sumMM = 0;
          double avgE = 0, avgEE = 0, avgM = 0, avgMM = 0;
          double avg_e = 0, avg_m = 0;
          double heatCap = 0, X = 0;

          mat threadLattice = initialize_lattice(L, T_phase, threadE, threadM, false);

          for (int burn = 0; burn < burnIn; burn++)
          {
            mcmc(threadLattice, generator, L, threadE, threadM, T_phase);
          }

          for (int n = burnIn; n < nCyclesPerThread; n++)
          {
            mcmc(threadLattice, generator, L, threadE, threadM, T_phase);

            sumE += threadE;
            sumEE += threadE * threadE;
            sumM += fabs(threadM);
            sumMM += threadM * threadM;
          }

          avg_e = sumE / (numSamples * N);
          avg_m = sumM / (numSamples * N);

          avgE = sumE / numSamples;
          avgEE = sumEE / numSamples;
          avgM = sumM / numSamples;
          avgMM = sumMM / numSamples;

          heatCap = (avgEE - avgE * avgE) / (N * T_phase * T_phase);
          X = (avgMM - avgM * avgM) / (N * T_phase);

          phaseRes(t, 0) = avg_e;
          phaseRes(t, 1) = avg_m;
          phaseRes(t, 2) = heatCap;
          phaseRes(t, 3) = X;

          auto t2 = chrono::high_resolution_clock::now();
          double duration = chrono::duration<double>(t2 - t1).count();

          cout << "Part " << t + 1
               << " finished in "
               << int(duration) / 60 << " min, "
               << int(duration) % 60 << " sec"
               << endl;
        }
      }
    }

    #endif

    phaseRes.save(outputFilename, raw_ascii);
  }

  else
  {
    results = cube(nCyclesPerThread, 6, threads);

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
        double avgE, avgEE, avgM, avgMM;
        double avg_e, avg_m;
        double heatCap, X;

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

            avg_e = sumE / (n * N);
            avg_m = sumM / (n * N);

            avgE = sumE / n;
            avgEE = sumEE / n;
            avgM = sumM / n;
            avgMM = sumMM / n;

            heatCap = (avgEE - avgE * avgE) / (T * T);
            X = (avgMM - avgM * avgM) / T;

            // Saving results to individual slices belonging to each thread
            results.slice(i).row(n - 1) = rowvec{avg_e, avg_m, heatCap, X, threadE / N, threadM / N};
          }
        }
      }
    }

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
      double avgE, avgEE, avgM, avgMM;
      double avg_e, avg_ee, avg_m, avg_mm;
      double heatCap, X;

      for (int n = 1; n <= nCyclesPerThread; n++)
      {
        mcmc(lattice, generator, L, serialE, serialM, T);

        sumE += serialE;
        sumEE += (serialE * serialE);
        sumM += fabs(serialM);
        sumMM += (serialM * serialM);

        avg_e = sumE / (n * N);
        avg_ee = sumEE / (n * N * N);
        avg_m = sumM / (n * N);
        avg_mm = sumMM / (n * N * N);

        avgE = sumE / n;
        avgEE = sumEE / n;
        avgM = sumM / n;
        avgMM = sumMM / n;

        heatCap = (avgEE - avgE * avgE) / (N * T * T);
        X = (avgMM - avgM * avgM) / (N * T);

        results.slice(0).row(n - 1) = rowvec{avg_e, avg_m, heatCap, X, serialE / N, serialM / N};
      }
    }

    #endif

    // Saving the results to matrix in raw ascii (mostly to avoid the
    // header of the matrix)
    results.save(outputFilename, raw_ascii);
  }

  return 0;
}
