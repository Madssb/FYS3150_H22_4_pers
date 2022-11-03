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

int main()
{
  // Containing the energies per spins
  vector<double> epsilon;

  int L = 20;        // Lattice length
  int N = L * L;    // No. of spin particles
  double T = 1.;    // Temperature [J/k]
  int nsteps = 100;  // No. of MCMC cycles to perform

  mat lattice = create_lattice(L);

  for (int n = 0; n < nsteps; n++)
  {
    mcmc(lattice, L, T);
    int tempTotE = totE(lattice, L);
    epsilon.push_back(1. / N * tempTotE);   // NB: Be aware of integer division!
  }
  // Writing to file
  int width = 15;

  ofstream outfile;
  outfile.open("energies.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width - 1) << "e(s)" << endl;

  for (int i = 0; i < epsilon.size(); i++)
  {
    outfile << setw(width) << epsilon[i] << endl;
  }

  outfile.close();

  return 0;
}
