#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"

//===============================================================
// CG solutions to Apsi = b
// see http://en.wikipedia.org/wiki/Conjugate_gradient_method
//===============================================================
//       x: The solution
//       b: The RHS vector
//      x0: An initial guess
//   gauge: The gauge field defining the operator
//   param: The parameter container

// Wilson g5Dg5D matrix inverter
//---------------------------------------------------------------

typedef vector<Complex> cvec;
typedef void (*op)(cvec&, const cvec&, const field<Complex>&);
int cg(cvec& P, const cvec& Q, const field<Complex>& gauge, op f);
