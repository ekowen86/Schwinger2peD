# (2+1)D-Schwinger Model

C++ code for a 2D Schwinger Model, with the option to add a 3rd dimension
for the gauge field with open boundary conditions and unit valued links in
the extra dimension. Fermions are restricted to the central slice.

## Features

The code is 'built for comfort, not for speed'. No attempts have been made to
optimise the code for performance.

### Actions

We use the Wilson gauge action throughout. At the present time, we offer

   2. 2 flavour Wilson

fermion actions.

### Simulation

We use Leapfrog HMC throughout the code. You have the option to perform dynamic
or quenched simulations, dictated from the command line. An example is given
in the `launcher.sh` file loacted in `wilson/`

### Measurements

At the preset time, only gauge measurements are implemented.

   1. Average plaquette
   2. Wilson loops (for Creutz ratios)
   3. Topological charge
   4. Pion correlation function

### Usage

Usual CMake build routine:

```
mkdir build
cd build
cmake <path/to/src>
make
```

To run the executable, use the provided `launcher.sh` file, edit as desired.

Happy Hacking!