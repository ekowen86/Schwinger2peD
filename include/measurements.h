#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"
// #include "utils.h"
// #include "blas.h"
// #include "inverters.h"

//-----------------------------------------------------------------------------------
// 2 Dimensional routines
//-----------------------------------------------------------------------------------

//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
void measWilsonLoops(field<Complex>& gauge, double plaq, int iter);

//Pion correlation function
//                              |----------------|
//                              |        |-------|---------|
//  < pi(x) | pi(0) > = < ReTr[dn*(x) g3 up(x) | dn(0) g3 up*(0)] >
//                    = < ReTr(g3 Gd[0,x] g3 Gu[x,0]) >
//                    = < ReTr(G*[x,0] G[x,0]) >
//
// using g3 G[x,0] g3 = G*[x,0] and Gu[x,0] \eq Gd[x,0]
//
// if H = Hdag, Tr(H * Hdag) = Sum_{n,m} (H_{n,m}) * (H_{n,m})^*,
// i.e., the sum of the modulus squared of each element
void measPionCorrelation(field<Complex>& gauge, int iter);
double measGaugeAction(field3D<Complex>& gauge);
Complex measPlaq(field<Complex>& gauge);
double measTopCharge(field<Complex>& gauge);
void measChiralCond(field<Complex>& gauge, int iter);
