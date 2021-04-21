#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"
// #include "utils.h"
// #include "blas.h"
// #include "inverters.h"

//-----------------------------------------------------------------------------------
// 2 Dimensional routines
//-----------------------------------------------------------------------------------

Complex measWilsonLoop(const field<Complex>& gauge, int a, int b);
// void measWilsonLoops(field<Complex>& gauge, double plaq, int iter);
void measPionCorrelation(field<Complex>& gauge, int iter);
double measGaugeAction(field3D<Complex>& gauge);
Complex measPlaq(field<Complex>& gauge);
double measFieldStrength(field<Complex>& gauge);
void drawInstantons(field<Complex>& gauge);
double measTopCharge(field<Complex>& gauge);
double measTopChargeBig(field<Complex>& gauge, int N);
void measChiralCond(field<Complex>& gauge, int iter);
