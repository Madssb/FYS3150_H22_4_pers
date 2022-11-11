/*
Main code for simulating the Ising model
*/

#include "lattice.hpp"
#include "mcmc.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace arma;

int main(int argc, const char* argv[])
{
  int L = atoi(argv[2]);        // Lattice length
  int nCycles = atoi(argv[3]);  // No. of MCMC cycles
  double T = atof(argv[4]);     // Temperature [J/k]
  int N = L * L;                // No. of spin particles

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

  mat lattice = initialize_lattice(L, T, E, M);

  ofstream outfile;
  outfile.open("energies.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "e"
          << setw(width) << "<e>"
          << setw(width) << "<e^2>"
          << setw(width) << "<|m|>"
          << setw(width) << "<m^2>"
          << endl;

  for (int n = 0; n < nCycles; n++)
  {
    mcmc(lattice, generator, L, E, M, T);

    sumE += E;
    sumEE += (E * E);
    sumM += fabs(M);
    sumMM += (M * M);

    outfile << setw(width) << n + 1
            << setw(width) << E / N
            << setw(width) << sumE / nCycles / N / N
            << setw(width) << sumEE / nCycles / N / N
            << setw(width) << sumM / nCycles / N / N
            << setw(width) << sumMM / nCycles / N / N
            << endl;
  }

  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "e"
          << setw(width) << "<e>"
          << setw(width) << "<e^2>"
          << setw(width) << "<|m|>"
          << setw(width) << "<m^2>"
          << endl;

  outfile.close();

  return 0;
}
