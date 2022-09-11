/*
Draft code for solving problem 2 within project 1
*/
// std::cout lives here
#include <iostream>
// std::exp lives here
#include <cmath>
// std::vector lives here
#include <vector>
// std::scientific and .precision() lives here
#include <iomanip>
/*
Analytic solution function:
u(x) = 1 - (1 - e^{-10})x - e^{-10x}
*/
double analytic_sol(double x)
{
  double res = 1 - (1 - std::exp(-10)) * x - std::exp(-10 * x);
  return res;
}

int main()
{
  int n_points = 100;                 // Number of points
  double step_length = 1. / n_points; // Step length

  std::vector<double> x(n_points + 1, 0); // Vector to contain x-values

  // int j = 0;  // Counter
  for (int i = 0; i <= n_points; i++) // Loop over and fill x-vector
  {
    x[i + 1] = x[i] + step_length;

    // std::cout.precision(3);
    // std::cout << std::scientific << x[j] << std::endl;   // Print out the vector to see if it works
    // j += 1;
  }

  std::vector<double> u(n_points + 1, 0); // Vector to contain analytic solution
  std::cout << "\nx-values\t"
            << "u-values\n " << std::endl;
  for (int k = 0; k <= n_points; k++) // Loop over and solve for u with analytic formula
  {
    u[k] = analytic_sol(x[k]);
    //amount of decimals are fixed
    std::cout.precision(3);
    //data is stored in scientific notation
    std::cout << std::scientific << x[k] << '\t' << u[k] << std::endl;
  }

  return 0;
}
