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
    const double tolerance = 1e-12;
    double k;
    double l;
    //k = k, l = l
    double &a_kl = max_element_offdiag(A0,k,l)

    while(a_kl > tolerance)
    {
        double tau = (A_M(k,k) - A_M(l,l)) / (2*A_M(k,l)) 
        double t;
        if (tau>0)
            tan = -tau + std::sqrt(1 + std::pow(tau,2));
        else (tau<0)
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

}


int main()
{
}
