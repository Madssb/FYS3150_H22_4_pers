#include <armadillo>
// std::pow and std::cos lives here
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    int N = 6;
    int n_rows = N;
    int n_cols = N;
    arma::mat tridiagonal_matrix(n_rows, n_cols);

    // finding the step length h
    double x_hat_max = 1.;
    double x_hat_min = 0.;
    double h = (x_hat_max - x_hat_min) / N;

    // defining the diagonal of tridiagonal matrix

    double d = 2 / std::pow(h, 2);
    for (int i = 0; i < 6; i++)
    {
        tridiagonal_matrix(i, i) = d;
    }
    // definding the supra and sub-diagonals
    double a = -1 / std::pow(h, 2);

    for (int i = 0; i < 5; i++)
    {
        // defining the subdiagonal
        tridiagonal_matrix(i + 1, i) = a;
        // defining the supradiagonal
        tridiagonal_matrix(i, i + 1) = a;
    }

    // finding the numerical eigenvalues and eigenvectors
    arma::vec num_eigvals;
    arma::mat num_eigvecs;
    arma::eig_sym(num_eigvals, num_eigvecs, tridiagonal_matrix);

    // normalizing eigenvectors in terms of units p=1
    arma::mat num_eigvecs_norm = arma::normalise(num_eigvecs, 1);

    // finding analytical eigenvalues
    arma::vec anal_eigvals = arma::vec(N - 1);
    double pi = 3.1415926535897;
    for (int i = 1; i < 6; i++)
    {
        anal_eigvals(i - 1) = d + 2 * a * std::cos((i * pi) / (N + 1));
    }

    // finding analytical eigenvectors
    arma::mat anal_eigvecs(n_rows, n_cols);
    for (int col = 0; col < N; col++)
    {
        for (int row = 0; row < N; row++)
        {
            std::cout << row << std::endl;
            anal_eigvecs(row, col) = std::sin((col+1) * (row+1) * pi / (N + 1));
        }
    }

    // normalize analytical eigenvectors
    arma::mat anal_eigvecs_norm = arma::normalise(anal_eigvecs, 1);

    //
    // section for all output
    //

    std::cout << "tridiagonal matrix:\n"
              << tridiagonal_matrix << std::endl;

    std::cout << "step length h = " << h << std::endl;
    std::cout << "sub/supra diagonal elements a = " << a << std::endl;
    std::cout << "diagonal elements d = " << d << std::endl;

    std::cout << "analytic eigvenvectors (normalized):\n"
              << anal_eigvecs_norm << std::endl;

    std::cout << "numerical eigenvectors (normalized):\n"
              << num_eigvecs_norm << std::endl;

    std::cout << "numerical eigvecs:\n"
              << num_eigvecs << std::endl;

    std::cout << "analytic eigvecs:\n"
              << anal_eigvecs << std::endl;

    std::cout << "analytical eigvals:\n"
              << anal_eigvals << std::endl;

    std::cout << "numerical eigvals:\n"
              << num_eigvals << std::endl;
}