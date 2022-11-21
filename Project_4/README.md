# The Ising model
In this project we study spin particles in a lattice. We compute expectation values for different size lattices and compare to analytical values. The crux of this project is the study of phase transition in the model.

We also take advantage of parallel computing, more specificly OpenMP.

Programs
---
`ising.cpp`
- The main program of this project
- Takes a total of 8 input arguments. The first two are the run command `./ising.exe` and the number `8`. The following 6 arguments are:

  - Number of threads
  - Lattice length (must be an even number)
  - Number of cycles to perform
  - Temperature in Kelvin
  - Weather the spins should be randomly set to +-1
    - Pass `ordered` to set all spins to +1
    - Pass `unordered` to randomly generate spins
  - The output filename (should be a .txt file)
  - Weather to study phase transition
    - Pass arbitrary word to run the program outisde of the phase transition block (a word should be passed, e.g. `normal`)
    - Pass `phasetransition` to enter the phase transition simulation. The program will ask for 3 values:
      - Lower boundary of the temperature to perform the simulation on
      - Upper boundary
      - Number of temperature elements to simulate on
      - At this point the program will print to the temperature range and the step length to screen

##### Building

This program can be built to run in either parallel or serial.

To build in **parallel**:
```
g++ -O2 ising.cpp src/lattice.cpp src/mcmc.cpp -I include -std=c++11 -Xpreprocessor -fopenmp -o ising.exe -larmadillo -lomp
```
Set number of threads for OpenMP to use (optional):
```
export OMP_NUM_THREADS=<no. of threads>
```
Run as describes above.

To build in **serial**:

First compile:
```
g++ -c -O2 ising.cpp src/lattice.cpp src/mcmc.cpp -I include -std=c++11
```
Then link:
```
g++ ising.o lattice.o mcmc.o -o ising.exe -larmadillo
```
Run as described above.

---
`plotting.py`
- Visualization program for this project 
