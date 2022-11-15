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
  const int L = atoi(argv[2]);                // Lattice length
  const int nCyclesPerThread = atoi(argv[3]); // No. of MCMC cycles per thread
  const double T = atof(argv[4]);             // Temperature [J/k]
  const string ordered_str = argv[5];         // Weather the lattice is ordered or not
  const string fileOutName = argv[6];         // Filename for output file
  const int N = L * L;                        // No. of spin particles

  bool ordered;

  if (ordered_str == "ordered")
  {
    ordered = true;
  }

  else
  {
    ordered = false;
  }

  int width = 15;

  ofstream outfile;

  outfile.open(fileOutName, ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "<e>"
          << setw(width) << "<|m|>"
          << setw(width) << "C_V"
          << setw(width) << "X"
          << endl;

  double E = 0;
  double M = 0;
  double sumE = 0;
  double sumEE = 0;
  double sumM = 0;
  double sumMM = 0;
  double avgEps, avgEps_sqrd, avgM, avgM_sqrd, heatCap, X;

  mat lattice = initialize_lattice(L, T, E, M, ordered);

  mt19937 generator;
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator.seed(seed);

  #ifdef _OPENMP
  {
    // Start parallell region
    #pragma omp parallel
    {
      const int n_threads = omp_get_num_threads();
      const int threadID = omp_get_thread_num();
      const double weight = 1. / n_threads;

      mt19937 threadGen;
      unsigned int threadSeed = seed + threadID;

      threadGen.seed(threadSeed);

      mat threadLattice = lattice;

      // Variables to be calculated before averaged across threads
      double threadE = E;
      double threadM = M;
      double threadSumE = 0;
      double threadSumEE = 0;
      double threadSumM = 0;
      double threadSumMM = 0;

      for (int n = 1; n <= nCyclesPerThread; n++)
      {
        mcmc(threadLattice, threadGen, L, threadE, threadM, T);

        // Average across threads
        #pragma omp critical
        {
          sumE += weight * threadE;
          sumEE += weight * (threadE * threadE);
          sumM += weight * fabs(threadM);
          sumMM += weight * (threadM * threadM);

          avgEps = sumE / (n * N);
          avgEps_sqrd = sumEE / (n * N * N);
          avgM = sumM / (n * N);
          avgM_sqrd = sumMM / (n * N * N);
          heatCap = N * (avgEps_sqrd - avgEps * avgEps) / (T * T);
          X = N * (avgM_sqrd - avgM * avgM) / T;

          outfile << setw(width) << n
                  << setw(width) << avgEps
                  << setw(width) << avgM
                  << setw(width) << heatCap
                  << setw(width) << X
                  << endl;
        }
      }
    }
  }

  #else
  {
    for (int n = 1; n <= nCyclesPerThread; n++)
    {
      mcmc(lattice, generator, L, E, M, T);

      sumE += E;
      sumEE += (E * E);
      sumM += fabs(M);
      sumMM += (M * M);

      avgEps = sumE / (n * N);
      avgEps_sqrd = sumEE / (n * N * N);
      avgM = sumM / (n * N);
      avgM_sqrd = sumMM / (n * N * N);
      heatCap = N * (avgEps_sqrd - avgEps * avgEps) / (T * T);
      X = N * (avgM_sqrd - avgM * avgM) / T;

      outfile << setw(width) << n
              << setw(width) << avgEps
              << setw(width) << avgM
              << setw(width) << heatCap
              << setw(width) << X
              << endl;
    }
  }
  #endif

  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "<e>"
          << setw(width) << "<|m|>"
          << setw(width) << "C_V"
          << setw(width) << "X"
          << endl;

  outfile.close();

  return 0;
}
