#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"
#include "inverters.h"

class leapfrogHMC {

public:

    // The inverter
    // inverterCG inv;

    // Objects for fermion force and guess tracking
    field<Complex> *phip;
    field<Complex> *g3Dphi;

    field3D<double> *mom3D;
    field<double> *mom2D;

    // 3D gauge field
    field3D<Complex> gauge3D_hmc;

    // 2D gauge field for center slice
    field<Complex> gauge2D;

    // Fermion fields are 2D
    field<Complex> *phi;
    field<Complex> *chi;

    // int hmc_count = 0;
    double exp_dH = 1.0;
    double dH = 0.0;
    double oldH = 0.0;
    double newH = 0.0;

    leapfrogHMC(param_t param);
    ~leapfrogHMC();
    int hmc(field3D<Complex> *gauge, bool noMetropolis = false);
    void trajectory(field3D<double> *mom, field3D<Complex> *gauge, field<Complex> *phi);
    void forceU(field3D<double> *fU, const field3D<Complex> *gauge);
    int forceD(field<double> *fD, field<Complex> *gauge, field<Complex> *phi);
    void update_mom(field3D<double> *fU, field<double> *fD, field3D<double> *mom, double dtau);
    void update_gauge(field3D<Complex> *gauge, field3D<double> *mom, double dtau);
};
