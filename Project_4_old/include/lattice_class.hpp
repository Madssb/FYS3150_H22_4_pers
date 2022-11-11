/*
Class containing the lattice
*/

#ifndef __Lattice_Class_hpp__
#define __Lattice_Class_hpp__

#include <string>
#include <armadillo>
#include <vector>

class Lattice
{
  public:
    // Member variables
    int L;                    // Length of lattice
    int N;                    // No. of spins = LxL
    double T;                 // Temperature [J/k]
    int nCycles;              // No. of MC cycles to perform
    int E;                    // Energy [J]
    int M;                    // Magnetization [?]
    double eps;               // Energy per spin [J]
    double m;                 // Magnetization per spin [?]
    double heatCap;           // Specific heat capacity [?]
    double suscept;           // Susceptibility [?]

  // Constructor
  Lattice(int L_in, double T_in, int nCycles_in);

  // Creating spins in lattice
  arma::mat create_spins();

  // Calculate energy of spin (i,j)
  double energySpin_ij(int i, int j);

  // Calculate total energy and magnetization of lattice
  void energy_magnetization();

  // Metropolis algorithm
  void metropolis();

  // MCMC cycle
  void mcmc();
};

#endif
