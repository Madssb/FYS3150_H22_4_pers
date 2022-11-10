/*
In this program we run sumulations on a NxN lattice containing spin particles
that have either spin up (+1) or spin down (-1).
*/

#include "lattice.hpp"
#include "mcmc.hpp"
#include "energy_magnetization.hpp"
#include "expectation_values.hpp"
#include "heat_cap_susceptibility.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <string>

using namespace std;
using namespace arma;

// Getting seed from chrono
unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();

// Declearing RNG
std::mt19937 generator;

int main(int argc, const char* argv[])
{
  // Setting seed to RNG
  generator.seed(seed);

  int L = atoi(argv[2]);        // Lattice length
  int nCycles = atoi(argv[3]);  // No. of MCMC cycles
  double T = atof(argv[4]);     // Temperature [J/k]
  int N = L * L;                // No. of spin particles

  // Containing the energies and  energies per spins
  vector<double> energy;
  vector<double> epsilon;
  vector<double> magnetization;
  vector<double> magnPerSpin;
  vector<double> expEpsilon;
  vector<double> expEpsilon_squared;
  vector<double> expMagn;
  vector<double> expMagn_squared;
  vector<double> heatCap;
  vector<double> suscept;

  mat lattice = create_lattice(L);

  double initEnergy = totE(lattice, L);  // Initial energy
  double initEpsilon = 1. / N * initEnergy;
  double initMagn = magnetizationLattice(lattice);
  double initMagnPerSpin = 1. / N * initMagn;

  energy.push_back(initEnergy);
  epsilon.push_back(initEpsilon);
  magnetization.push_back(initMagn);
  magnPerSpin.push_back(initMagnPerSpin);

  for (int n = 0; n < nCycles; n++)
  {
    mcmc(lattice, L, T, generator);

    int latticeEnergy = totE(lattice, L);
    energy.push_back(latticeEnergy);
    epsilon.push_back(1. / N * latticeEnergy);

    double magnLattice = magnetizationLattice(lattice);
    magnetization.push_back(magnLattice);
    magnPerSpin.push_back(1. / N * magnLattice);
  }
  expected_values(expEpsilon, expEpsilon_squared, epsilon, magnPerSpin, expMagn,
                  expMagn_squared, initMagnPerSpin, initEpsilon);

  heatCap_suscept(expEpsilon, expEpsilon_squared, expMagn, expMagn_squared,
                  heatCap, suscept, N, T);

  // Writing to file
  int width = 15;

  ofstream outfile;
  outfile.open("energies.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "e(s)"
          << setw(width) << "<e>"
          << setw(width) << "<e^2>"
          << setw(width) << "C_V"
          << setw(width) << "m(s)"
          << setw(width) << "<|m|>"
          << setw(width) << "<m^2>"
          << setw(width) << "X"
          << endl;

  for (int i = 0; i < epsilon.size(); i++)
  {
    outfile << setw(width) << i
            << setw(width) << epsilon[i]
            << setw(width) << expEpsilon[i]
            << setw(width) << expEpsilon_squared[i]
            << setw(width) << heatCap[i]
            << setw(width) << magnPerSpin[i]
            << setw(width) << expMagn[i]
            << setw(width) << expMagn_squared[i]
            << setw(width) << suscept[i]
            << endl;
  }
  
  outfile << "#" << setw(width - 1) << "nCycles"
          << setw(width) << "e(s)"
          << setw(width) << "<e>"
          << setw(width) << "<e^2>"
          << setw(width) << "C_V"
          << setw(width) << "m(s)"
          << setw(width) << "<|m|>"
          << setw(width) << "<m^2>"
          << setw(width) << "X"
          << endl;

  outfile.close();

  return 0;
}
