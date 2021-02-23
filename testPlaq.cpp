// testPlaq.cpp

#include <random>
#include "schwinger2peD_internal.h"
#include "utils.h"
#include "hmc.h"
#include "io.h"
#include "measurements.h"

using namespace std;

int main(int argc, char **argv) {

    if (argc != 13) {
        printf("Error: Invalid command line arguments.\n");
        return 1;
    }

    int Nx = atoi(argv[1]); printf("Nx: %d\n", Nx);
    int Ny = atoi(argv[2]); printf("Ny: %d\n", Ny);
    int Nz = atoi(argv[3]); printf("Nz: %d\n", Nz);
    double betaMin = stod(argv[4]); printf("betaMin: %.3f\n", betaMin);
    double betaMax = stod(argv[5]); printf("betaMax: %.3f\n", betaMax);
    double betaInc = stod(argv[6]); printf("betaInc: %.3f\n", betaInc);
    double eps3 = stod(argv[7]); printf("eps3: %.3f\n", eps3);

    // number of thermalization sweeps
    int nTherm = atoi(argv[8]); printf("nTherm: %d\n", nTherm);
    // number of trajectories per checkpoint
    int nSkip = atoi(argv[9]); printf("nSkip: %d\n", nSkip);
    // total number of trajectories (after thermalization)
    int nTraj = atoi(argv[10]); printf("nTraj: %d\n", nTraj);

    param_t p;

    p.beta = betaMin;
    p.betaZ = betaMin * eps3;

    p.n_step = atoi(argv[11]); printf("hmcSteps: %d\n", p.n_step);
    p.tau = stod(argv[12]); printf("hmcTau: %.3f\n", p.tau);

    p.dynamic = false;
    p.max_iter_cg = 10000;
    p.eps = 1e-16;

    // Lattice size
    p.Nx = Nx;
    p.Ny = Ny;
    p.Nz = Nz;
    p.zCenter = (p.Nz - 1) / 2;
    mt19937 gen(1234);
    p.gen = &gen;


    for (double beta = betaMin; beta <= betaMax; beta += betaInc) {
        p.beta = beta;
        p.betaZ = beta * eps3;
        printf("beta: %.3lf\n", p.beta);
        printf("betaZ: %.3lf\n", p.betaZ);

        field<Complex> gauge2D(p);
        field3D<Complex> gauge3D(p);

        leapfrogHMC HMCStep(p);
        int accept = 0;
        double dH_sum = 0.0;
        double exp_dH_sum = 0.0;

        hotStart(gauge3D); // hot start

        for (int n = 0; n <= (nTherm + nTraj); n++) {

            // no hmc for initial configuration
            if (n != 0) {
                // no metropolis for first half of thermalization
                accept += HMCStep.hmc(gauge3D, n < (nTherm / 2));
                dH_sum += HMCStep.dH;
                exp_dH_sum += HMCStep.exp_dH;
            }

            if (n % nSkip) continue;

            // write trajectory data to file and console
            char path[50];
            sprintf(path, "plaq/plaq_%d_%d_%d_%d.dat",
                Nx, Nz, int(round(p.beta * 1000)), int(round(eps3 * 1000)));
            FILE* file = fopen(path, "a");
            fprintf(file, "%06d", n);
            printf("%06d", n);

            // print plaquette
            extract2DSlice(gauge2D, gauge3D, p.zCenter);
            double plaq = real(measPlaq(gauge2D));
            fprintf(file, " %.12e", plaq);
            printf(" %.8f", plaq);

            // print 1x1 wilson loop
            double w11 = real(measWilsonLoop(gauge2D, 1, 1));
            fprintf(file, " %.12e", w11);
            printf(" %.8f", w11);

            // print 1x2 and 2x1 wilson loops
            double w12 = real(measWilsonLoop(gauge2D, 1, 2));
            fprintf(file, " %.12e", w12);
            printf(" %.8f", w12);
            double w21 = real(measWilsonLoop(gauge2D, 2, 1));
            fprintf(file, " %.12e", w21);
            printf(" %.8f", w21);

            // print 2x2 wilson loop
            double w22 = real(measWilsonLoop(gauge2D, 2, 2));
            fprintf(file, " %.12e", w22);
            printf(" %.8f", w22);

            // print hmc statistics
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
    }

    return 0;
}
