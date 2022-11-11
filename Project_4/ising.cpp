/*
Main code for simulating the Ising model
*/

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
  int L = atoi(argv[2]);        // Lattice length
  int nCycles = atoi(argv[3]);  // No. of MCMC cycles
  double T = atof(argv[4]);     // Temperature [J/k]
  string ordered_str = argv[5]; // Weather the lattice is ordered or not
  int N = L * L;                // No. of spin particles

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

  mt19937 generator;
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator.seed(seed);

  double E = 0;
  double M = 0;
  double sumE = 0;
  double sumEE = 0;
  double sumM = 0;
  double sumMM = 0;
  double heatCap = 0;
  double X = 0;

  mat lattice = initialize_lattice(L, T, E, M, ordered);

  ofstream outfile;
  string outfileName;

  if (ordered)
  {
    outfileName = "ordered_" + to_string(L) +
                  "by" + to_string(L) + "_lattice.txt";
  }

  if (!ordered)
  {
    outfileName = "unordered_" + to_string(L) +
                   "by" + to_string(L) + "_lattice.txt";
  }

  outfile.open(outfileName, ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "<e>"
          << setw(width) << "<|m|>"
          << setw(width) << "C_V"
          << setw(width) << "X"
          << endl;

  double avgEps, avgEps_sqrd, avgM, avgM_sqrd;

  for (int n = 1; n <= nCycles; n++)
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

    outfile << setw(width) << n
            << setw(width) << avgEps
            << setw(width) << avgM
            << setw(width) << N * (avgEps_sqrd - avgEps * avgEps) / (T * T)
            << setw(width) << N * (avgM_sqrd - avgM * avgM) / T
            << endl;
  }

  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "<e>"
          << setw(width) << "<|m|>"
          << setw(width) << "C_V"
          << setw(width) << "X"
          << endl;

  outfile.close();

  if (N == 4)
  {
    cout << "Computed values after " << nCycles << " cycles:" << endl;
    cout << "<e>: " << avgEps << endl;
    cout << "<|m|>: " << avgM << endl;
    cout << "C_V: " << N * (avgEps_sqrd - avgEps * avgEps) / (T * T) << endl;
    cout << "X: " << N * (avgM_sqrd - avgM * avgM) / T << endl;
  }

  return 0;
}
