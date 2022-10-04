// Main code of project 3
#include "particle.hpp"
#include "penning_trap.hpp"

using namespace std;

int main()
{
  /// Defining constants
  const double T = 9.64852558 * 10;   // [u/(mu*s)e ]
  const double V = 9.64852558 * 1e7;  // [u(mu*m)^2/(mu*s)^2e]
  const double B_0 = 1. * T;          // Magnetic field constant
  const double V_0 = 10. * V;         // Electric field constant
  const double d = 1e4;               // [mu*m]

  PenningTrap trap = PenningTrap(B_0, V_0, d);




  return 0;
}
