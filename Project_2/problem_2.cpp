#include <armadillo>
#include <cmath> //std::pow lives here
int N = 6;
//finding the step length
double x_hat_max = 1
double x_hat_min = 0
float h = (x_hat_max - x_hat_min)/N;
//initializing tridiagonal matrix
arma::mat tridiag(N, N);
//defining the diagonal of the matrix
double a = -1/std::pow(h,2)
for (i = 0; i < 6, i++)
{
    tridiag(i,i) = a
}