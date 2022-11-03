/*
Header file for the functions to calculate the energy and magnetization of
the lattice
*/

#ifndef __Energy_Magnetization_hpp__
#define __Energy_Magnetization_hpp__

#include <armadillo>

// Calculate the total energy of the lattice
int totE(arma::mat& lattice, int L);

// Compute the magnetization of lattice
int magnetizationLattice(arma::mat& lattice);

#endif
