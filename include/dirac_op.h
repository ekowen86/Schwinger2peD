#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"

void Dpsi(field<Complex> *psi2, const field<Complex> *psi1, const field<Complex> *gauge);
void g3Dpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);
void Ddagpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);
void DdagDpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge);
void g3psi(field<Complex> *out, const field<Complex> *in);
void g2psi(field<Complex> *out, const field<Complex> *in);
void g1psi(field<Complex> *out, const field<Complex> *in);

void _DdagDpsi(vector<Complex>& out, const vector<Complex>& in, const field<Complex>& gauge);
void _DdagDpsiImp(vector<Complex>& out, const vector<Complex>& in, const field<Complex>& gauge);
