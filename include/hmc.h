#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"
#include "inverters.h"

class leapfrogHMC {

public:

    // Objects for fermion force and guess tracking
    field<Complex> phip;
    field<Complex> g3Dphi;

    field3D<double> mom3D;

    // 3D gauge field
    field3D<Complex> gauge3D;

    // 2D gauge field for center slice
    field<Complex> gauge2D;

    // Fermion fields are 2D
    field<Complex> phi;
    field<Complex> chi;

    // force fields
    field3D<double> fU; // gauge force
    field<double> fD; // fermion force

    double exp_dH = 1.0;
    double dH = 0.0;
    double oldH = 0.0;
    double newH = 0.0;

    leapfrogHMC(param_t p);
    int hmc(field3D<Complex>& oldGauge, bool noMetropolis = false);
    void trajectory();
    void forceU();
    int forceD();
    void update_mom(double dtau);
    void update_gauge(double dtau);
};
