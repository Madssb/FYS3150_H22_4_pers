/*
In this program we run sumulations on a NxN lattice containing spin particles
that have either spin up (+1) or spin down (-1).
*/

#include "lattice.hpp"
#include "mcmc.hpp"
#include "energy_magnetization.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

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
  int nsteps = atoi(argv[3]);   // No. of MCMC cycles
  double T = atof(argv[4]);     // Temperature [J/k]
  int N = L * L;                // No. of spin particles

  // Containing the energies and  energies per spins
  vector<double> energy;
  vector<double> epsilon;
  vector<double> magnetization;
  vector<double> magnPerSpin;

  mat lattice = create_lattice(L);

  double initEnergy = totE(lattice, L);  // Initial energy
  double initEpsilon = 1. / N * initEnergy;
  double initMagn = magnetizationLattice(lattice);
  double initMagnPerSpin = 1. / N * initMagn;

  energy.push_back(initEnergy);
  epsilon.push_back(initEpsilon);
  magnetization.push_back(initMagn);
  magnPerSpin. push_back(initMagnPerSpin);

  for (int n = 0; n < nsteps; n++)
  {
    int deltaE = mcmc(lattice, L, T, generator);
    energy.push_back(energy[n] + deltaE);
    epsilon.push_back(epsilon[n] + 1. / N * deltaE);
    magnetization.push_back(magnetizationLattice(lattice));
    magnPerSpin.push_back(1. / N * magnetizationLattice(lattice));
  }
  // Writing to file
  int width = 15;

  ofstream outfile;
  outfile.open("energies.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width - 1) << "e(s)"
          << setw(width) << "E(s)"
          << setw(width) << "M(s)"
          << setw(width) << "m(s)"
          << endl;

  for (int i = 0; i < epsilon.size(); i++)
  {
    outfile << setw(width) << epsilon[i]
            << setw(width) << energy[i]
            << setw(width) << magnetization[i]
            << setw(width) << magnPerSpin[i]
            << endl;
  }

  outfile.close();

  return 0;
}
