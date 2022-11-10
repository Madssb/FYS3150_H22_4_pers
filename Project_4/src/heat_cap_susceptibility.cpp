/*
Definition of the function to compute the specific heat cap. and susceptibility
*/

#include "heat_cap_susceptibility.hpp"

void heatCap_suscept(std::vector<double>& expEps, std::vector<double>& expEps_sqrd,
                     std::vector<double>& expMagn, std::vector<double>& expMagn_sqrd,
                     std::vector<double>& heatCap, std::vector<double>& suscept,
                     int noSpins, double T)
{
  int N = expEps.size();

  double heatCapFactor = noSpins / (T * T);
  double susceptFactor = noSpins / T;

  double epsVar, magnVar;

  heatCap.push_back(heatCapFactor * (expEps_sqrd[0] - (expEps[0] * expEps[0])));
  suscept.push_back(susceptFactor * (expMagn_sqrd[0] - (expMagn[0] * expMagn[0])));

  for (int n = 1; n < N; n++)
  {
    epsVar = expEps_sqrd[n] - (expEps[n] * expEps[n]);
    magnVar = expMagn_sqrd[n] - (expMagn[n] * expMagn[n]);

    heatCap.push_back(heatCapFactor * epsVar);
    suscept.push_back(susceptFactor * magnVar);
  }
}
