// testWilsonFlow.cpp

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
    p.Nx = 48;
    p.Ny = 48;
    p.Nz = 1;
    p.zCenter = (p.Nz - 1) / 2;

    p.beta = 4.0;
    p.betaZ = 0.0;
    p.m = 0.1;

    p.n_step = 30;
    p.tau = 1.0;

    p.dynamic = false;
    p.lockedZ = true;
    p.max_iter_cg = 10000;
    p.eps = 1e-16;

    printf("Generating initial lattice configuration...\n");
    field3D<Complex> gauge3D(p);
    field<Complex> gauge2D(p);
    hotStart(gauge3D);

    // do some hmc steps
    leapfrogHMC HMCStep(p);
    for (int i = 0; i < 200; i++) HMCStep.hmc(gauge3D, true);
    extract2DSlice(gauge2D, gauge3D, p.zCenter);

    // measure field strength as a function of wilson flow time
    double dt = 0.02;
    for (int i = 0; i <= 1000; i++) {
        double t = double(i) * dt;
        if (i != 0) wilsonFlow(gauge2D, dt);
        double fs = measFieldStrength(gauge2D) * 8.0;
        printf("%.02f %.12f %.12f\n", t, fs, fs * t * t);
    }

    return 0;
}
