//
// schwinger2peD.cpp
// program to generate lattice configurations
//

#include <random>
#include "schwinger2peD_internal.h"
#include "utils.h"
#include "hmc.h"
#include "io.h"
#include "measurements.h"

using namespace std;

int quickTopoCharge(field3D<Complex>& gauge);
void writeTopoCharge(field3D<Complex>& gauge, int n);

int main(int argc, char **argv) {

    if (argc != 12) {
        printf("Error: Invalid command line arguments.\n");
        return 1;
    }

    int Nx = atoi(argv[1]); printf("Nx: %d\n", Nx);
    int Ny = atoi(argv[2]); printf("Ny: %d\n", Ny);
    int Nz = atoi(argv[3]); printf("Nz: %d\n", Nz);
    double beta = stod(argv[4]); printf("beta: %.3f\n", beta);
    double eps3 = stod(argv[5]); printf("eps3: %.3f\n", eps3);
    double m = stod(argv[6]); printf("m: %.3f\n", m);

    // number of thermalization sweeps
    int nTherm = atoi(argv[7]); printf("nTherm: %d\n", nTherm);
    // number of trajectories per checkpoint
    int nSkip = atoi(argv[8]); printf("nSkip: %d\n", nSkip);
    // total number of trajectories (after thermalization)
    int nTraj = atoi(argv[9]); printf("nTraj: %d\n", nTraj);

    param_t p;

    p.beta = beta;
    p.betaZ = beta * eps3;
    p.m = m;

    p.n_step = atoi(argv[10]); printf("hmcSteps: %d\n", p.n_step);
    p.tau = stod(argv[11]); printf("hmcTau: %.3f\n", p.tau);

    p.dynamic = true;
    p.lockedZ = true;
    p.max_iter_cg = 1000;
    p.eps = 1e-16;

    // Lattice size
    p.Nx = Nx;
    p.Ny = Ny;
    p.Nz = Nz;
    p.zCenter = (p.Nz - 1) / 2;
    p.gen = mt19937(1234);

    field<Complex> gauge2D(p);
    field3D<Complex> gauge3D(p);

    gaussStart(&gauge3D); // hot start

    leapfrogHMC HMCStep(p);
    int accept = 0;
    double dH_sum = 0.0;
    double exp_dH_sum = 0.0;
    int nStuck = 0;
    int prevTopoCharge = 0;

    // thermalize
    for (int n = 0; n <= (nTherm + nTraj); n++) {

        // no hmc for initial trajectory
        // no metropolis for first half of thermalization
        if (n) {
            accept += HMCStep.hmc(&gauge3D, n < (nTherm / 2));
            dH_sum += HMCStep.dH;
            exp_dH_sum += HMCStep.exp_dH;
            int topoCharge = quickTopoCharge(gauge3D);
            if (topoCharge == prevTopoCharge) nStuck++;
        }

        // write topological charge every trajectory
        // writeTopoCharge(gauge3D, n);

        if (n % nSkip) continue;

        // write checkpoint
        char ckpointPath[50];
        sprintf(ckpointPath, "ckpoint/ckpoint.%d", n);
        writeGaugeBinary(gauge3D, ckpointPath);

        extract2DSlice(&gauge2D, &gauge3D, p.zCenter);
        measVacuumTrace(&gauge2D, n);
        measPionCorrelation(&gauge2D, n);
        writeTopoCharge(gauge3D, n);

        printf("%06d", n);
        for (int z = 0; z < p.Nz; z++) {
            // print plaquette on each slice
            extract2DSlice(&gauge2D, &gauge3D, z);
            printf(" %.8f", real(measPlaq(&gauge2D)));
        }
        printf(" %.2f", double(nStuck) / double(nSkip)); nStuck = 0;
        printf(" %+.8f", dH_sum / double(nSkip)); dH_sum = 0.0;
        printf(" %.8f", exp_dH_sum / double(nSkip)); exp_dH_sum = 0.0;
        printf(" %.2f", double(accept) / double(nSkip)); accept = 0;
        printf("\n");
    }

    return 0;
}

int quickTopoCharge(field3D<Complex>& gauge) {
    field<Complex> coolGauge(gauge.p);
    extract2DSlice(&coolGauge, &gauge, gauge.p.zCenter);
    coolLink(coolGauge, 10);
    return int(measTopCharge(&coolGauge));
}

void writeTopoCharge(field3D<Complex>& gauge, int n) {

    // open data file
    char path[50];
    sprintf(path, "topo_charge/topo_charge.%d", n);
    ofstream file;
    file.open(path, ofstream::out | ofstream::trunc);

    // extract each 2D slice into an array
    vector<field<Complex>> coolGauge;
    for (int z = 0; z < gauge.p.Nz; z++) {
        coolGauge.push_back(field<Complex>(gauge.p));
        extract2DSlice(&coolGauge[z], &gauge, z);
    }

    double topoCharge[gauge.p.Nz];
    for (int nCool = 0; nCool <= 10; nCool++) {

        // smear the gauge field and measure topological charge
        for (int z = 0; z < gauge.p.Nz; z++) {
            if (nCool) coolLink(coolGauge[z]);
            topoCharge[z] = measTopCharge(&coolGauge[z]);
        }

        // number of cooling steps
        file << fixed << nCool;

        for (int z = 0; z < gauge.p.Nz; z++) {
            // round topological charge to the nearest integer
            file << " " << int(round(topoCharge[z]));
        }
        file << scientific << setprecision(12);
        for (int z = 0; z < gauge.p.Nz; z++) {
            // exact measured value of topological charge
            file << " " << topoCharge[z];
        }
        file << endl;
    }
}
