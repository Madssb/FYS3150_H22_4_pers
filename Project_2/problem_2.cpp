//
// This program generates a tridiagonal matrix A, computes the eigenvalues and eigenvectors by analytical means, 
// and by more numerical means, using armadillo properties. 
// Both are compared to find the relative error of armadillo-generated eigenvalues and eigenvectors.
//


#include <armadillo>
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    //
    // initializing the tridiagonal (a,d,a) matrix A
    //

    constexpr int N = 6;
    constexpr int n_rows = N;
    constexpr int n_cols = N;
    arma::mat tridiagonal_matrix(n_rows, n_cols);

    // finding the step length h
    constexpr double x_hat_max = 1.;
    constexpr double x_hat_min = 0.;
    constexpr int n_steps = N + 1;
    constexpr double h = (x_hat_max - x_hat_min) / n_steps;

    // defining the main, sub and supra-diagonals
    constexpr double d = 2 / std::pow(h, 2);
    constexpr double a = -1 / std::pow(h, 2);
    tridiagonal_matrix.diag(1).fill(a);
    tridiagonal_matrix.diag(0).fill(d);
    tridiagonal_matrix.diag(-1).fill(a);

    //
    // computing the eigenvalues and eigenvectors for A by numerical means
    //

    arma::vec num_eigvals;
    arma::mat num_eigvecs;
    arma::eig_sym(num_eigvals, num_eigvecs, tridiagonal_matrix);

    // normalizing eigenvectors in terms of units p=1
    arma::mat num_eigvecs_norm = arma::normalise(num_eigvecs);

    //
    // solving for the eigenvalues and eigenvalues for A by analytic means
    //

    arma::vec anal_eigvals = arma::vec(N);
    constexpr double pi = 3.1415926535897;
    for (int i = 0; i < 6; i++)
    {
        anal_eigvals(i) = d + 2 * a * std::cos(((i + 1) * pi) / (N + 1));
    }

    // finding analytical eigenvectors
    arma::mat anal_eigvecs(n_rows, n_cols);
    for (int col = 0; col < N; col++)
    {
        for (int row = 0; row < N; row++)
        {
            anal_eigvecs(row, col) = std::sin((col + 1) * (row + 1) * pi / (N + 1));
        }
    }

    // normalize analytical eigenvectors
    arma::mat anal_eigvecs_norm = arma::normalise(anal_eigvecs);

    // manual rescaling such that signs match for the numeric and analytic eigenvectors
    std::vector<int> indices_for_rows_that_must_be_inverted = {1, 5};
    for (int index : indices_for_rows_that_must_be_inverted)
    {
        anal_eigvecs_norm.col(index) *= -1;
    }

    //
    // finding relative errors for the computed eigenvalues & eigenvectors
    //

    arma::vec rel_err_eigvals = (num_eigvals - anal_eigvals) / anal_eigvals;
    arma::mat rel_err_eigvecs = (num_eigvecs_norm - anal_eigvecs_norm) / anal_eigvecs_norm;

    //
    // section for output of all results
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

    std::cout << "numerical eigenvectors:\n"
              << num_eigvecs << std::endl;

    std::cout << "analytic eigvectors:\n"
              << anal_eigvecs << std::endl;

    std::cout << "analytical eigenvalues:\n"
              << anal_eigvals << std::endl;

    std::cout << "numerical eigenvalues:\n"
              << num_eigvals << std::endl;

    std::cout << "relative error for eigenvalues:\n"
              << rel_err_eigvals << std::endl;

    std::cout << "relative error for eigenvectors:\n"
              << rel_err_eigvecs << std::endl;
}