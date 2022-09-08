# include <iostream>
# include <cmath>
# include <vector>
/*
Problem 7 on FYS3150 project 1. Implementing the Thomas-algorithm as general as possible.
*/

// Function that prints all elements in vector
void print_all(std::vector<double> y){
  for (double x : y){
    std::cout << x << " ";
  }
  std::cout << "" << std::endl;
}

int main(){
  // The dimension of matrix, i.e B: n x n
  int n = 4;
  // Initial values for vector we are solving, end points
  double v0 = 1.;
  double vN = 2.;

  /*
  Defining and initializing all vectors needed to solve. Use n elements for
  a and c (upper and lower diagonals of the tridiagonal matrix) as the 0'th element are thrown away when indexing in the loop. Therefore
  we need the n'th element, and the 0'th element is never used. This is to
  keep the indexing convention of the algorithm.
  */
  std::vector<double> a(n, -1.);
  std::vector<double> b(n, 2.);
  std::vector<double> c(n, -1.);
  std::vector<double> g = {1., 2., 3., 4.};
  std::vector<double> bt(n);
  std::vector<double> gt(n);
  std::vector<double> v(n);
  std::vector<double> vt(n+2);
  vt[0] = v0;
  vt[n+1] = vN;
  bt[0] = b[0];
  gt[0] = g[0];

  // Loop solving for bt and gt using Thomas algorithm
  for (int i=1; i<n; i++){
    bt[i] = b[i] - a[i]/bt[i-1]*c[i-1];
    gt[i] = g[i] - a[i]/bt[i-1]*gt[i-1];
  }
  // Initializing v, where last element is gt / bt
  v[n-1] = gt[n-1] / bt[n-1];
  // Solving v using Thomas algorithm, backwars substitution
  for (int i=n-2; i>=0; i--){
    v[i] = (gt[i] - c[i]*v[i+1]) / bt[i];
  }
  // Appending solutions from v in our original vector. This has size n+2 as our inital conditions is added.
  for (int i=1; i<=n; i++){
    vt[i] = v[i-1];
  }

// Printing all vector-elements to check if everything went okay
print_all(a);
print_all(b);
print_all(c);
print_all(g);
print_all(bt);
print_all(gt);
print_all(v);
print_all(vt);

  return 0;
}
