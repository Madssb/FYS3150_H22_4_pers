#include <armadillo>
// std::pow and std::cos lives here
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    // initializing tridiagonal matrix
    int N = 6;
    int n_rows = N;
    int n_cols = N;
    arma::mat tridiag(n_rows, n_cols);

    // finding the step length h
    double x_hat_max = 1.;
    double x_hat_min = 0.;
    double h = (x_hat_max - x_hat_min) / N;
    std::cout << "step length h = " << h <<std::endl;
    // defining the diagonal of tridiagonal matrix
    double a = -1 / std::pow(h, 2);
    std::cout << "diagonal elements a = " << a << std::endl;
    for (int i = 0; i < 6; i++)
    {
        tridiag(i, i) = a;
    }
    // definding the supra and sub-diagonals
    double d = 2 / std::pow(h, 2);
    std::cout << "sub/supra diagonal elements d = " << d << std::endl;
    for (int i = 0; i < 5; i++)
    {
        // defining the subdiagonal
        tridiag(i + 1, i) = d;
        // defining the supradiagonal
        tridiag(i, i + 1) = d;
    }
    std::cout << "tridiagonal matrix:\n"
              << tridiag << std::endl;

    // finding the numerical eigenvalues and eigvectors
    arma::vec num_eigval;
    arma::mat num_eigvec;
    arma::eig_sym(num_eigval, num_eigvec, tridiag);
    std::cout << "numerical eigvals:\n"
              << num_eigval << std::endl;
    std::cout << "numerical eigvecs:\n"
              << num_eigvec << std::endl;

    // finding analytical eigenvalues
    arma::vec anal_eigval = arma::vec(N - 1);
    double pi = 3.1415926535897;
    for (int i = 1; i < 6; i++)
    {
        anal_eigval(i - 1) = d + 2 * a * std::cos((i * pi) / (N + 1));
    }

    std::cout << "analytical eigvals:\n"
              << anal_eigval << std::endl;
}