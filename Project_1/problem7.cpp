# include <iostream>
# include <vector>
# include <cmath>
#include <string>
#include <fstream>
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

std::vector<double> linspace(double start, double end, int points){
  std::vector<double> linvec(points);
  for (int i=0; i<points; i++){
    linvec[i] = start + i*(end - start) / (points-1);
  }
  return linvec;
}

std::vector<double> source_term(std::vector<double> x){
  std::vector<double> f(x.size());
  for (double i=0; i<x.size(); i++){
    f[i] = 100 * std::exp(-10*x[i]);
  }
  return f;
}

int main(){
  // The dimension of matrix, i.e B: n x n
  int n = 10000;

  /*
  Defining and initializing all vectors needed to solve. Use n elements for
  a and c (upper and lower diagonals of the tridiagonal matrix) as the 0'th element are thrown away when indexing in the loop. Therefore
  we need the n'th element, and the 0'th element is never used. This is to
  keep the indexing convention of the algorithm.
  */
  std::vector<double> x = linspace(0., 1., n+2);
  double h = x[1] - x[0];
  std::vector<double> f = source_term(x);
  std::vector<double> a(n, -1.);
  std::vector<double> b(n, 2.);
  std::vector<double> c(n, -1.);
  std::vector<double> g(n);
  std::vector<double> bt(n);
  std::vector<double> gt(n);
  std::vector<double> v(n);
  std::vector<double> vt(n+2);
  // Initial values for vector we are solving, end points
  double v0 = 0.;
  double vN = 0.;
  // Initializing g, the vector containing our solution term
  g[0] = h*h*f[1] + v0;
  g[n-1] = h*h*f[n-2] + vN;
  for (int i=1; i<g.size()-1; i++){
    g[i] = h*h * f[i];
  }
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
// print_all(a);
// print_all(b);
// print_all(c);
// print_all(g);
// print_all(bt);
// print_all(gt);
// print_all(v);
// print_all(vt);

// Set a filename
std::string filename = "problem7_data.txt";

// Create and open the output file. Or, technically, create
// an "output file stream" (type std::ofstream) and connect
// it to our filename.
std::ofstream ofile;
ofile.open(filename);

// Send some text to this output file
for (int i=0; i<vt.size(); i++){
  ofile << x[i] << " " << vt[i] << std::endl;
}

// Close the output file
ofile.close();

  return 0;
}
