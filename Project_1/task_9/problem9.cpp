#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
/*
Problem 9 on FYS3150 project 1. Implementing the Thomas-algorithm as special as possible for the
tridiaginal matrix A with signature (-1, 2, -1).
Then
bt_i = 2 - 1 / bt_(i-1)  since a, c = -1, a*c = 1, b_i = 2.
gt_i = g_i + gt_(i-1)/bt_(i-1) since a = -1
v_i = (gt_i + v_(i+1)) / bt_i since c = -1
*/

// Function that prints all elements in vector
void print_all(std::vector<double> y)
{
  for (double x : y)
  {
    std::cout << x << " ";
  }
  std::cout << "" << std::endl;
  return;
}

// Function that acts like Python linspace
std::vector<double> linspace(double start, double end, int points)
{
  std::vector<double> linvec(points);
  for (int i = 0; i < points; i++)
  {
    linvec[i] = start + i * (end - start) / (points - 1);
  }
  return linvec;
}

std::vector<double> source_term(std::vector<double> x)
{
  std::vector<double> f(x.size());
  for (double i = 0; i < x.size(); i++)
  {
    f[i] = 100 * std::exp(-10 * x[i]);
  }
  return f;
}

int main()
{
  int k = 6;
  std::vector<int> n_vec(k);
  n_vec[0] = 10;
  for (int j = 0; j < k - 1; j++)
  {
    n_vec[j + 1] = n_vec[j] * 10;
  }
  for (int N_step : n_vec)
  {

    auto t1 = std::chrono::high_resolution_clock::now();

    int n_step = N_step;
    // The dimension of matrix, i.e B: n x n. Is number of points (n_step + 1) minus 2 => n_step - 1
    int n = n_step - 1;

    /*
    Defining and initializing all vectors needed to solve. Remember that a, b, c in the special case have signature (-1, 2, -1),
    but generally can be any arbitrary vector a, b, c with arbitrary entries.
    Now we alter the algorithm to be specific for a = c = -1, and b = 2 as it will reduce FLOPs.
    */
    std::vector<double> x = linspace(0., 1., n_step + 1);
    double h = x[1] - x[0];
    std::vector<double> f = source_term(x);
    std::vector<double> g(n);
    std::vector<double> bt(n);
    std::vector<double> gt(n);
    std::vector<double> v(n);
    std::vector<double> vt(n_step + 1);
    // Initial values for vector we are solving, end points
    double v0 = 0.;
    double vN = 0.;
    // Initializing g, the vector containing our solution term
    g[0] = h * h * f[1] + v0;
    g[n - 1] = h * h * f[n - 2] + vN;
    for (int i = 1; i < g.size() - 1; i++)
    {
      g[i] = h * h * f[i];
    }
    vt[0] = v0;
    vt[n + 1] = vN;
    bt[0] = 2;
    gt[0] = g[0];

    // Loop solving for bt and gt using Thomas algorithm
    for (int i = 1; i < bt.size(); i++)
    {
      /*
      Here, b_i = 2 and a_i*c_(i-1) = 1, so we reduce the first line by 1 FLOP.
      In the second line we have that a_i = -1, so the sign is flipped. This allow us to
      move gt_(i-1) on top of bt_(i-1), reducing the line by 1 FLOP.
      */
      bt[i] = 2 - 1 / bt[i - 1];
      gt[i] = g[i] + gt[i - 1] / bt[i - 1];
    }
    // Initializing v, where last element is gt / bt
    v[n - 1] = gt[n - 1] / bt[n - 1];
    // Solving v using Thomas algorithm, backwards substitution
    for (int i = v.size() - 2; i >= 0; i--)
    {
      /*
      Here c_(i-1) = -1, flipping the sign. This reduces the line by 1 FLOP
      */
      v[i] = (gt[i] + v[i + 1]) / bt[i];
    }
    // Appending solutions from v in our original vector. This has size n+2 as our inital conditions is added.
    for (int i = 1; i <= v.size(); i++)
    {
      vt[i] = v[i - 1];
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
    std::cout << "n_step = " << n_step << " used " << duration_seconds << " seconds." << std::endl;

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
    std::string filename = "problem9_data_n=";
    std::string txt = ".txt";
    filename += std::to_string(n_step);
    filename += txt;

    // Create and open the output file. Or, technically, create
    // an "output file stream" (type std::ofstream) and connect
    // it to our filename.
    std::ofstream ofile;
    ofile.open(filename);

    // Send some text to this output file
    for (int i = 0; i < vt.size(); i++)
    {
      ofile << x[i] << " " << vt[i] << std::endl;
    }

    // Close the output file
    ofile.close();
  }
  return 0;
}
