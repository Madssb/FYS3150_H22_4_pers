# Project 4
##### Notes:
- The `<ordered/unordered>` should be passed as either `ordered` or `unordered`.
- Temperature must be passed in unit Kelvin

##### Parallel coding
To build:
```
g++ -O2 ising.cpp src/lattice.cpp src/mcmc.cpp -I include -std=c++11 -Xpreprocessor -fopenmp -o ising.exe -larmadillo -lomp
```
Set no. of threads:
```
export OMP_NUM_THREADS=<no. of threads>
```
To run:
```
./ising.exe 7 <no. of threads> <lattice length> <cycles per thread> <temperature> <ordered/unordered> <filename.txt>
```
---
##### Serial coding:
First compile:
```
g++ -c -O2 ising.cpp src/lattice.cpp src/mcmc.cpp -I include -std=c++11
```
Then link:
```
g++ ising.o lattice.o mcmc.o -o ising.exe -larmadillo
```
Run:
```
./ising.exe 7 1 <lattice length> <cycles per thread> <temperature> <ordered/unordered> <filename.txt>
```
