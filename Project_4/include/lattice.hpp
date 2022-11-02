/*
Header file containing declarations of functions to initialize a LxL lattice
containing random spins (+-1), and doing calculations on it
*/

#ifndef __Lattice_hpp__
#define __Lattice_hpp__
#include <armadillo>

// Creating LxL lattice
arma::mat create_lattice(int L);

// Compute the energy of spin (i, j)
int energySpin_ij(arma::mat& lattice, int L, int i, int j);

// Compute the magnetization of lattice
int magnetizationLattice(arma::mat& lattice, int L);

#endif
