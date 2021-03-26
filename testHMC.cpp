// testHMC.cpp

#include <random>
#include "schwinger2peD_internal.h"
#include "utils.h"
#include "hmc.h"
#include "measurements.h"
#include "dirac_op.h"

using namespace std;

int main(int argc, char **argv) {

    param_t p;
    mt19937 gen(1234);
    p.gen = &gen;

    // Lattice size
    p.Nx = 32;
    p.Ny = 32;
    p.Nz = 5;
    p.zCenter = (p.Nz - 1) / 2;
    p.beta = 4.0;
    p.betaZ = 0.5;
    p.m = 0.1;
    p.dynamic = true;

    p.n_step = 20;
    p.tau = 1.0;

    leapfrogHMC HMCStep(p);

    hotStart(HMCStep.gauge3D);
    gaussReal(HMCStep.mom3D);

    // Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
    gaussComplex(HMCStep.chi);

    // Create pseudo fermion field, phi = D * chi
    extract2DSlice(HMCStep.gauge2D, HMCStep.gauge3D, p.zCenter);
    g3Dpsi(HMCStep.phi, HMCStep.chi, HMCStep.gauge2D);
    // blas::caxpy(-I * sqrt(p.musq), HMCStep.chi.data, HMCStep.phi.data);

    // measure initial hamiltonian
    double H0 = 0.0;
    H0 += blas::norm2(HMCStep.mom3D.data) * 0.5;
    H0 += measGaugeAction(HMCStep.gauge3D);
    if (p.dynamic) {
        H0 += real(blas::cDotProd(HMCStep.chi.data, HMCStep.chi.data));
    }

    // do forward trajectory
    HMCStep.trajectory();

    // measure final hamiltonian
    double H1 = 0.0;
    H1 += blas::norm2(HMCStep.mom3D.data) * 0.5;
    H1 += measGaugeAction(HMCStep.gauge3D);
    if (p.dynamic) {
        extract2DSlice(HMCStep.gauge2D, HMCStep.gauge3D, p.zCenter);
        cg(HMCStep.chi.data, HMCStep.phi.data, HMCStep.gauge2D, &_DdagDpsi);
        H1 += real(blas::cDotProd(HMCStep.chi.data, HMCStep.phi.data));
    }

    // do backward trajectory
    p.tau = -p.tau;
    blas::ax(-1.0, HMCStep.mom3D.data);
    HMCStep.trajectory();

    // measure (new) initial hamiltonian
    double H2 = 0.0;
    H2 += blas::norm2(HMCStep.mom3D.data) * 0.5;
    H2 += measGaugeAction(HMCStep.gauge3D);
    if (p.dynamic) {
        extract2DSlice(HMCStep.gauge2D, HMCStep.gauge3D, p.zCenter);
        cg(HMCStep.chi.data, HMCStep.phi.data, HMCStep.gauge2D, &_DdagDpsi);
        H2 += real(blas::cDotProd(HMCStep.chi.data, HMCStep.phi.data));
    }

    printf("H0 = %.12lf\n", H0);
    printf("H1 = %.12lf\n", H1);
    printf("H2 = %.12lf\n", H2);
    printf("percent error = %.12e\n", abs(H2 - H0) / H0);
    printf("\n");

    return 0;
}
