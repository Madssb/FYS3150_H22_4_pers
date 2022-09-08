# include <iostream>
# include <cmath>
# include <vector>
// # define print(x) (std:: cout << x << std::endl)

void print_all(std::vector<double> y){
  for (double x : y){
    std::cout << x << " ";
  }
  std::cout << "" << std::endl;
}

int main(){
  int n = 4;
  double v0 = 1.;
  double vN = 2.;

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

  for (int i=1; i<n; i++){
    bt[i] = b[i] - a[i]/bt[i-1]*c[i-1];
    gt[i] = g[i] - a[i]/bt[i-1]*gt[i-1];
  }

  v[n-1] = gt[n-1] / bt[n-1];
  for (int i=n-2; i>=0; i--){
    v[i] = (gt[i] - c[i]*v[i+1]) / bt[i];
  }

  for (int i=1; i<=n; i++){
    vt[i] = v[i-1];
  }

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
