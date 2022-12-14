Solving the time-dependent Schrödinger equation
---
This folder contains the programs used to solve the time-dependent Schrödinger equation, using the Crank-Nicolson scheme. We simulate a wave packet hitting a single, double, and triple-slit and measure the probability after impact.

**How to use**

The program `potential.py` will create the needed potentials for the simulations. It can create a single, double, and triple-slit potential barrier. After setting the arguments as wanted, run with `python3 potential.py`.

The program `main.cpp` is built with
```
g++ -O2 main.cpp src/matrix.cpp -I include -std=c++11 -o main.exe -larmadillo
```
In order to run the program, one has to have predefined configuration files as arguments. One of these has to be a potential file (Created with e.g. the `pptential.py` program), and the other a `.txt` file with the following setup:

- Spatial step length
- Time step length
- Total simulation time
- Center coordinate of initial wave packet in x-direction
- Width of initial wave packet in x-direction
- Momentum in x-direction
- Center coordinate of initial wave packet in y-direction
- Width of initial wave packet in y-direction
- Momentum in y-direction

The file must be written on a **single row** in this exact order.

Run the program with
```
./main.cpp 4 <configuration file> <potential file> <out-file name>
```
**Note** that the config and potential files must be passed with their respective file extension. The `<out-file name>` is passed without extension and will be set to `.bin`.

For visualisation, the program `plotting.py` can be used. The structure of this folder is set so that **all programs should be run from this folder** (i.e. from the folder `Project_5`). The figures generated from running
```
python3 plotting/plotting.py
```
from this folder will be saved in the folder `Project_5/figures/`.
