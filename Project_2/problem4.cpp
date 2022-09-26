//
//Program implements Jacobi's rotation algorithm
//


#include <armadillo>
#include <iostream>
#include <cmath> 


double jacobi_rotate(const arma::mat symmetric_matrix)
{
    arma::mat A_M = symmetric_matrix;
    arma::mat A_Mp1;
    arma::mat R0(5, 5 , fill::eye);
    arma::mat S_m;

    //k = k, l = l
    double &max_element_offdiag = max_element_offdiag(A0,k,l);

    while(a_kl > tolerance)
    {
        double tau = (A_M(k,k) - A_M(l,l)) / (2*A_M(k,l)) 
        double tan;
        if (tau>0)
            tan = -tau + std::sqrt(1 + std::pow(tau,2));
        else
            tan  = -tau - std::sqrt(1 + std::pow(tau,2));

        double cos = 1/std::sqrt(1+ std::pow(tan,2));
        double sin = cos*tan;
        S_m.eye();
        S_m(k,k) = cos;
        S_m(l,k) = -sin;
        S_m(k,l) = sin;
        S_m(l,l) = cos;
        // evauluating A = STAS 
    }

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat &symmetric_matrix, arma::mat &rotation_matrix, int k, int l)
{
    double tau = (symmetric_matrix(l,l) - symmetric_matrix(k,k))/2*symmetric_matrix(k,l)
    double tan;

    if (tau>0)
    {
        tan = -tau + std::sqrt(1 + std::pow(tau,2));
    }
    else
    {
        tan  = -tau - std::sqrt(1 + std::pow(tau,2));
    }
    double cos = 1/std::sqrt(1+ std::pow(tan,2));
    double sin = cos*tan;    
    double temp_ik;
    //transforming symmetric matrix (STAS)
    symmetric_matrix(k,k) = symmetric_matrix(k,k)*std::pow(cos,2) - 2*symmetric_matrix(k,l)*cos + symmetric_matrix(k,k)*std::pow(sin,2);
    symmetric_matrix(l,l) = symmetric_matrix(k,k)*std::pow(cos,2) + 2*symmetric_matrix(k,l)*cos + symmetric_matrix(k,k)*std::pow(sin,2);
    symmetric_matrix(k,l) = 0;
    symmetric_matrix(l,k) = 0;
    for(int i = 0; i < symmetric_matrix.n_cols; i++)
    {
        if(i != k && i != l)
        {
            temp_ik = symmetric_matrix(i,k);
            symmetric_matrix(i,k) = symmetric_matrix(i,k)*cos - symmetric_matrix(i,l)*sin;
            symmetric_matrix(k,i) = symmetric_matrix(i,k);
            symmetric_matrix(i,l) = symmetric_matrix(i,l)*cos + temp_ik*sin;
            symmetric_matric(l,i) = symmetric_matrix(i,l);
        }
    }
    rotation_matrix(r)
}


int main()
{
    return 0;
}
