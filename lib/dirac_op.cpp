#include "dirac_op.h"
#include "blas.h"


void Dpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge) {

    double m0 = gauge->p.m;
    double r = 1.0;
    double constant = 2 * r + m0;

    //Sum over 0,1 directions.
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
// #pragma omp parallel for
    for (int x = 0; x < Nx; x++) {
        int xp1 = (x+1)%Nx;
        int xm1 = (x-1+Nx)%Nx;
        for (int y = 0; y < Ny; y++) {
            int yp1 = (y+1)%Ny;
            int ym1 = (y-1+Ny)%Ny;

            Complex gauge_px = gauge->read(x,y,0);
            Complex gauge_py = gauge->read(x,y,1);
            Complex gauge_mx = conj(gauge->read(xm1,y,0));
            Complex gauge_my = conj(gauge->read(x,ym1,1));
            Complex in0 = in->read(x,y,0);
            Complex in0_px = in->read(xp1,y,0);
            Complex in0_py = in->read(x,yp1,0);
            Complex in0_mx = in->read(xm1,y,0);
            Complex in0_my = in->read(x,ym1,0);
            Complex in1 = in->read(x,y,1);
            Complex in1_px = in->read(xp1,y,1);
            Complex in1_py = in->read(x,yp1,1);
            Complex in1_mx = in->read(xm1,y,1);
            Complex in1_my = in->read(x,ym1,1);

            // antiperiodic boundary conditions in y direction
            if (yp1 == 0) {
                in0_py *= -1.0;
                in1_py *= -1.0;
            }

            if (y == 0) {
                in0_my *= -1.0;
                in1_my *= -1.0;
            }

            //upper
            out->write(x,y,0, constant * in0 -
            0.5*(gauge_px * (r*in0_px -   in1_px) +
                 gauge_py * (r*in0_py + I*in1_py) +
                 gauge_mx * (r*in0_mx +   in1_mx) +
                 gauge_my * (r*in0_my - I*in1_my)));

            //lower
            out->write(x,y,1, constant * in1 -
            0.5*(gauge_px * (  -in0_px + r*in1_px) +
                 gauge_py * (-I*in0_py + r*in1_py) +
                 gauge_mx * (   in0_mx + r*in1_mx) +
                 gauge_my * ( I*in0_my + r*in1_my)));
        }
    }
}

void g3psi(field<Complex> *out, const field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
// #pragma omp parallel for
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0,  in->read(x,y,0));
      out->write(x,y,1, -in->read(x,y,1));
    }
}

void g2psi(field<Complex> *out, const field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
// #pragma omp parallel for
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0, -I*in->read(x,y,1));
      out->write(x,y,1,  I*in->read(x,y,0));
    }
}

void g1psi(field<Complex> *out, field<Complex> *in){

  int Nx = in->p.Nx;
  int Ny = in->p.Ny;
// #pragma omp parallel for
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++) {
      out->write(x,y,0, in->read(x,y,1));
      out->write(x,y,1, in->read(x,y,0));
    }
}

void g3Dpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge){

  field<Complex> *temp = new field<Complex>(in->p);
  Dpsi(temp, in, gauge);
  g3psi(out, temp);
  delete temp;
}

void Ddagpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge){

  field<Complex> *temp = new field<Complex>(in->p);
  g3psi(out, in);
  Dpsi(temp, out, gauge);
  g3psi(out, temp);
  delete temp;
}


void DdagDpsi(field<Complex> *out, const field<Complex> *in, const field<Complex> *gauge) {

  field<Complex> *temp = new field<Complex>(in->p);
  Dpsi(temp, in, gauge);
  g3psi(out, temp);
  Dpsi(temp, out, gauge);
  g3psi(out, temp);
  delete temp;
}

void _DdagDpsi(vector<Complex>& out, const vector<Complex>& in, const field<Complex>& gauge) {
    field<Complex> outField(gauge.p);
    field<Complex> inField(gauge.p);
    blas::copy(outField.data, out);
    blas::copy(inField.data, in);
    DdagDpsi(&outField, &inField, &gauge);
    blas::copy(out, outField.data);
}

void _DdagDpsiImp(vector<Complex>& out, const vector<Complex>& in, const field<Complex>& gauge) {
    field<Complex> outField(gauge.p);
    field<Complex> inField(gauge.p);
    blas::copy(outField.data, out);
    blas::copy(inField.data, in);
    DdagDpsi(&outField, &inField, &gauge);
    blas::caxpy(gauge.p.musq, in, outField.data, out);
}
