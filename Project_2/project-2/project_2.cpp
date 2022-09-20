#include<iostream>
#include "header.hpp"

int main(){

  double h = 1./6;
  double a = -1./(h*h);
  double d = 2./(h*h);

  ana_vs_arma_test(d, a, 6.);

  max_offdiag_symm_test();

  return 0;
}
