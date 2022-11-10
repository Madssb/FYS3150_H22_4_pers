/*
Header file fot the functions computing the specific heat capacity and
susceptibility
*/

#ifndef __Heat_Cap_Susceptibility_hpp__
#define __Heat_Cap_Susceptibility_hpp__

#include<vector>

void heatCap_suscept(std::vector<double>& expEps, std::vector<double>& expEps_sqrd,
                     std::vector<double>& expMagn, std::vector<double>& expMagn_sqrd,
                     std::vector<double>& heatCap, std::vector<double>& suscept,
                     int noSpins, double T);

#endif
