// testIO.cpp

#include <random>
#include "schwinger2peD_internal.h"
#include "utils.h"
#include "hmc.h"
#include "measurements.h"
#include "io.h"

using namespace std;

int main(int argc, char **argv) {

    param_t p;
    mt19937 gen(1234);
    p.gen = &gen;

    // Lattice size
    p.Nx = 32;
    p.Ny = 32;
    p.Nz = 1;
    p.zCenter = (p.Nz - 1) / 2;

    p.beta = 0.0;
    p.betaZ = 0.0;
    p.m = 0.1;

    p.n_step = 20;
    p.tau = 0.5;

    p.dynamic = true;
    p.lockedZ = true;
    p.max_iter_cg = 10000;
    p.eps = 1e-16;

    string filenameBin = "testIO.bin";
    string filenameRng = "testIO.rng";

    printf("Generating initial lattice configuration...\n");
    field3D<Complex> gauge3D(p);
    field<Complex> gauge2D(p);
    hotStart(&gauge3D);

    // do some hmc steps
    leapfrogHMC HMCStep(p);
    for (int i = 0; i < 20; i++) HMCStep.hmc(gauge3D, true);

    // save the lattice and rng
    printf("Writing initial lattice and rng configuration...\n");
    writeGaugeBinary(gauge3D, filenameBin);
    writeRngState(gauge3D, filenameRng);

    extract2DSlice(&gauge2D, &gauge3D, p.zCenter);
    Complex plaq0 = measPlaq(&gauge2D);
    printf("Average plaquette = %.12lf + %.12lf i\n", real(plaq0), imag(plaq0));

    // do some more hmc steps
    printf("Updating and measuring original lattice...\n");
    for (int i = 0; i < 20; i++) HMCStep.hmc(gauge3D, true);
    extract2DSlice(&gauge2D, &gauge3D, p.zCenter);
    Complex plaq1 = measPlaq(&gauge2D);
    printf("Average plaquette = %.12lf + %.12lf i\n", real(plaq1), imag(plaq1));

    // save the lattice and rng
    printf("Restoring initial lattice and rng configuration...\n");
    readGaugeBinary(gauge3D, filenameBin);
    readRngState(gauge3D, filenameRng);

    extract2DSlice(&gauge2D, &gauge3D, p.zCenter);
    Complex plaq2 = measPlaq(&gauge2D);
    printf("Average plaquette = %.12lf + %.12lf i\n", real(plaq2), imag(plaq2));

    // do some more hmc steps
    printf("Updating and measuring restored lattice...\n");
    for (int i = 0; i < 20; i++) HMCStep.hmc(gauge3D, true);
    extract2DSlice(&gauge2D, &gauge3D, p.zCenter);
    Complex plaq3 = measPlaq(&gauge2D);
    printf("Average plaquette = %.12lf + %.12lf i\n", real(plaq3), imag(plaq3));

    printf("Error = %.12lf + %.12f i\n", real(plaq3 - plaq1), imag(plaq3 - plaq1));

    remove(filenameBin.c_str());
    remove(filenameRng.c_str());

    return 0;
}
