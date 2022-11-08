/*
Header file for the functions to calculate the expected values
*/

#ifndef __Expectation_Values_hpp__
#define __Expectation_Values_hpp__
#include <vector>

std::vector<double> expected_energy(std::vector<double> epsilonEnergies, double T);

std::vector<double> expected_magnetization(std::vector<double> magnetization, double T);

#endif
