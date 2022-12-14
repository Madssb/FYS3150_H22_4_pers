/*
This program will calculate the machine precision.
*/

#include <iostream>

using namespace std;

void machine_epsilon(double n)
{
  int k = 0;
  double eps;

  while ((1 + n) != 1)
  {
    eps = n;
    n /= 2;
    k += 1;
  }

  cout << "Found machine epsilon to be " << eps
       << " after " << k << " iterations."
       << endl;
}

int main()
{
  machine_epsilon(.5);

  return 0;
}
