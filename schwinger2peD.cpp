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

// wilson loops take a long time to measure
#define MEASURE_WILSON_LOOPS 0

using namespace std;

int quickTopoCharge(field3D<Complex>& gauge);
void writeTopoCharge(field3D<Complex>& gauge, int n);
void doWilsonFlow(const field<Complex>& gauge, int n);

int main(int argc, char **argv) {

    if (argc < 12) {
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

    p.dynamic = false;
    p.lockedZ = true;
    p.max_iter_cg = 10000;
    p.eps = 1e-16;

    // Lattice size
    p.Nx = Nx;
    p.Ny = Ny;
    p.Nz = Nz;
    p.zCenter = (p.Nz - 1) / 2;
    mt19937 gen(1234);
    p.gen = &gen;

    field<Complex> gauge2D(p);
    field3D<Complex> gauge3D(p);

    int nCkpoint = 0;
    if (argc == 12) {
        hotStart(gauge3D); // hot start
    } else {
        nCkpoint = atoi(argv[12]);
        printf("ckPoint: %d\n", nCkpoint);
        nTherm = 0; // no thermalization for checkpoint start
        char ckpointPath[50];
        sprintf(ckpointPath, "ckpoint/ckpoint.%d", nCkpoint);
        readGaugeBinary(gauge3D, ckpointPath);
        char rngPath[50];
        sprintf(rngPath, "ckpoint/rng.%d", nCkpoint);
        readRngState(gauge3D, rngPath);
    }

    leapfrogHMC HMCStep(p);
    int accept = 0;
    double dH_sum = 0.0;
    double exp_dH_sum = 0.0;
    int nStuck = 0;
    int prevTopoCharge = quickTopoCharge(gauge3D);

    for (int n = nCkpoint; n <= (nTherm + nTraj + nCkpoint); n++) {

        // no hmc for initial configuration
        if (n != nCkpoint) {
            // no metropolis for first half of thermalization
            accept += HMCStep.hmc(gauge3D, n < (nTherm / 2));
            dH_sum += HMCStep.dH;
            exp_dH_sum += HMCStep.exp_dH;
            int topoCharge = quickTopoCharge(gauge3D);
            if (topoCharge == prevTopoCharge) nStuck++;
            prevTopoCharge = topoCharge;
        } else if (nCkpoint) {
            // don't write initial output for checkpoint start
            continue;
        }

        writeTopoCharge(gauge3D, n);

        if (n % nSkip) continue;

        // write checkpoint
        char ckpointPath[50];
        sprintf(ckpointPath, "ckpoint/ckpoint.%d", n);
        writeGaugeBinary(gauge3D, ckpointPath);
        char rngPath[50];
        sprintf(rngPath, "ckpoint/rng.%d", n);
        writeRngState(gauge3D, rngPath);

        extract2DSlice(gauge2D, gauge3D, p.zCenter);
        measChiralCond(gauge2D, n);
        measPionCorrelation(gauge2D, n);
        doWilsonFlow(gauge2D, n);

        // write trajectory data to file and console
        FILE* file = fopen("plaq.dat", "a");
        fprintf(file, "%06d", n);
        printf("%06d", n);
        for (int z = 0; z < p.Nz; z++) {
            // print plaquette on each slice
            extract2DSlice(gauge2D, gauge3D, z);
            double plaq = real(measPlaq(gauge2D));
            fprintf(file, " %.12e", plaq);
            printf(" %.8f", plaq);
        }
        fprintf(file, " %.2f", double(nStuck) / double(nSkip));
        printf(" %.2f", double(nStuck) / double(nSkip));
        nStuck = 0;

        fprintf(file, " %+.12e", dH_sum / double(nSkip));
        printf(" %+.8f", dH_sum / double(nSkip));
        dH_sum = 0.0;

        fprintf(file, " %.12e", exp_dH_sum / double(nSkip));
        printf(" %.8f", exp_dH_sum / double(nSkip));
        exp_dH_sum = 0.0;

        fprintf(file, " %.2f", double(accept) / double(nSkip));
        printf(" %.2f", double(accept) / double(nSkip));
        accept = 0;

        fprintf(file, "\n");
        printf("\n");
        fclose(file);
    }

    return 0;
}

int quickTopoCharge(field3D<Complex>& gauge) {
    field<Complex> coolGauge(gauge.p);
    extract2DSlice(coolGauge, gauge, gauge.p.zCenter);
    coolLink(coolGauge, 10);
    return int(round(measTopCharge(coolGauge)));
}

void writeTopoCharge(field3D<Complex>& gauge, int n) {

    // open data file
    char path[50];
    FILE* file = fopen("topo_charge.dat", "a");
    fprintf(file, "%06d", n);

    // extract each 2D slice into an array
    vector<field<Complex>> coolGauge;
    for (int z = 0; z < gauge.p.Nz; z++) {
        coolGauge.push_back(field<Complex>(gauge.p));
        extract2DSlice(coolGauge[z], gauge, z);
    }

    double topoCharge[gauge.p.Nz];
    int topoChargeInt[gauge.p.Nz];

    // smear the gauge field and measure topological charge
    for (int z = 0; z < gauge.p.Nz; z++) {
        coolLink(coolGauge[z], 10);
        topoCharge[z] = measTopCharge(coolGauge[z]);
        topoChargeInt[z] = int(round(topoCharge[z]));
    }

    for (int z = 0; z < gauge.p.Nz; z++) {
        // round topological charge to the nearest integer
        fprintf(file, " %+d", topoChargeInt[z]);
    }
    for (int z = 0; z < gauge.p.Nz; z++) {
        // exact measured value of topological charge
        fprintf(file, " %.16e", topoCharge[z]);
    }
    fprintf(file, "\n");
    fclose(file);
}

void doWilsonFlow(const field<Complex>& gauge, int n) {

    // copy the gauge field
    field<Complex> gauge_wf(gauge.p);
    gauge_wf.copy(gauge);

    double fs;
    int loopMax = min(gauge.p.Nx, gauge.p.Ny) / 2;
    int nLoops = loopMax * 3 + 1;
    std::vector<double> wLoops(nLoops);

    double t_wf_max = 30.0;
    double dt = 0.02;
    int n_wf = int(ceil(t_wf_max / dt));

    char path[50];
    FILE* file;

    for (int i = 0; i <= n_wf; i++) {
        double t = i * dt;
        if (i) wilsonFlow(gauge_wf, dt);

        // measure field strength and wilson loops
        double fs = measFieldStrength(gauge_wf);

        // write field strength to file
        sprintf(path, "wf/field_strength.%d", n);
        file = fopen(path, "a");
        fprintf(file, "%.02lf", t);
        fprintf(file, " %.12e", fs);
        fprintf(file, "\n");
        fclose(file);

#if MEASURE_WILSON_LOOPS
        // measure wilson loops in this order
        // (1,1), (1,2)
        // (2,1), (2,2), (2,3)
        //        (3,2), (3,3), (3,4) ...
        for (int l = 1; l <= loopMax; l++) {
            if (l != 1) {
                wLoops[l * 3 - 1] = real(measWilsonLoop(gauge_wf, l, l - 1));
            }
            wLoops[l * 3] = real(measWilsonLoop(gauge_wf, l, l));
            if (l != loopMax) {
                wLoops[l * 3 + 1] = real(measWilsonLoop(gauge_wf, l, l + 1));
            }
        }

        // write wilson loops to file
        sprintf(path, "wf/wilson_loops.%d", n);
        file = fopen(path, "a");
        fprintf(file, "%.02lf", t);
        for (int l = 0; l < nLoops; l++) {
            fprintf(file, " %.12e", wLoops[l]);
        }
        fprintf(file, "\n");
        fclose(file);
#endif // MEASURE_WILSON_LOOPS
    }
}
