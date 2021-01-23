#pragma once

#include "schwinger2peD_internal.h"
#include "utils.h"
#include "dirac_op.h"
#include "measurements.h"

class leapfrogHMC {
  
private:

  // The inverter
  inverterCG *inv;

  // Objects for fermion force and guess tracking
  field<Complex> *phip;
  field<Complex> *g3Dphi;

  field3D<double> *mom3D;
  field<double> *mom2D;
  field3D<Complex> *gauge3D_old;
  field<Complex> *gauge2D;

  // Fermion fields are 2D
  field<Complex> *phi;
  field<Complex> *chi;
  
public:

  int hmc_count = 0;
  double exp_dH_ave = 0.0;
  double dH_ave = 0.0;
  
  leapfrogHMC(param_t param);
  int hmc(field3D<Complex> *gauge, int iter);
  void trajectory(field3D<double> *mom, field3D<Complex> *gauge, field<Complex> *phi, int iter);
  void forceU(field3D<double> *fU, const field3D<Complex> *gauge);
  int forceD(field<double> *fD, field<Complex> *gauge, field<Complex> *phi, int iter);
  void update_mom(field3D<double> *fU, field<double> *fD, field3D<double> *mom, double dtau);
  void update_gauge(field3D<Complex> *gauge, field3D<double> *mom, double dtau);
  
  ~leapfrogHMC();
  
};
