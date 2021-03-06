#include "inverters.h"
#include "dirac_op.h"
#include "utils.h"
#include "measurements.h"
#include "hmc.h"

//2D HMC Routines
//---------------------------------------------------------------------

leapfrogHMC::leapfrogHMC(param_t p) :
    gauge3D(p),
    gauge2D(p),
    phip(p),
    g3Dphi(p),
    mom3D(p),
    phi(p),
    chi(p),
    fU(p),
    fD(p) {

};

int leapfrogHMC::hmc(field3D<Complex>& oldGauge, bool noMetropolis) {

    int accept = 0;

    gauge3D.copy(oldGauge);
    gaussReal(mom3D);

    if (gauge3D.p.dynamic) {
        //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
        gaussComplex(chi);
        //Create pseudo fermion field, phi = D * chi
        extract2DSlice(gauge2D, gauge3D, gauge3D.p.zCenter);
        g3Dpsi(phi, chi, gauge2D);
        blas::caxpy(-I * sqrt(gauge3D.p.musq), chi.data, phi.data);
    }

    double oldH = 0.0;
    oldH += blas::norm2(mom3D.data) * 0.5;
    oldH += measGaugeAction(gauge3D);
    if (gauge3D.p.dynamic) {
        oldH += real(blas::cDotProd(chi.data, chi.data));
    }

    trajectory();

    double newH = 0.0;
    newH += blas::norm2(mom3D.data) * 0.5;
    newH += measGaugeAction(gauge3D);
    if (gauge3D.p.dynamic) {
        extract2DSlice(gauge2D, gauge3D, gauge3D.p.zCenter);
        cg(chi.data, phi.data, gauge2D, &_DdagDpsiImp);
        newH += real(blas::cDotProd(chi.data, phi.data));
    }

    dH = newH - oldH;
    exp_dH = exp(-dH);

    if (dH <= 0 || gauge3D.rand_double() < exp_dH) accept = 1;
    if (accept || noMetropolis) oldGauge.copy(gauge3D);

    return accept;
}

void leapfrogHMC::trajectory() {

    // time step
    double dtau = gauge3D.p.tau / gauge3D.p.n_step;

    //Initial half step.
    //P_{1/2} = P_0 - dtau/2 * (fU - fD)
    forceU();
    // ave_iter += forceD(fD, phi, gauge2D);
    forceD();
    update_mom(0.5 * dtau);

    // intermediate steps
    for (int k = 1; k < gauge3D.p.n_step; k++) {

      //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
      update_gauge(dtau);

      //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
      forceU();
      forceD();

      update_mom(dtau);
    }

    //Final half step.
    //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
    update_gauge(dtau);

    //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
    forceU();
    forceD();
    update_mom(0.5 * dtau);
}

// void leapfrogHMC::trajectory() {
//
//     double dtau1 = gauge3D.p.tau / 6.0;
//
//     update_mom_fermion(0.5 * dtau1);
//
//     for (int k = 1; k < 6; k++) {
//         trajectory2();
//         update_mom_fermion(dtau1);
//     }
//
//     trajectory2();
//     update_mom_fermion(0.5 * dtau1);
// }
//
// void leapfrogHMC::trajectory2() {
//
//     double dtau2 = gauge3D.p.tau / (12.0 * 6.0);
//
//     update_mom_gauge(0.5 * dtau2);
//
//     for (int k = 1; k < 12; k++) {
//
//         update_gauge(dtau2);
//         update_mom_gauge(dtau2);
//     }
//
//     update_gauge(dtau2);
//     update_mom_gauge(0.5 * dtau2);
// }

void leapfrogHMC::forceU() {

    int xp1, xm1, yp1, ym1, zp1, zm1;

    int Nx = gauge3D.p.Nx;
    int Ny = gauge3D.p.Ny;
    int Nz = gauge3D.p.Nz;
    double betaZ = gauge3D.p.betaZ;
    double f;
    Complex u1, u2, u3, u4;

    for (int x = 0; x < Nx; x++) {
        xp1 = (x + 1) % Nx;
        xm1 = (x - 1 + Nx) % Nx;
        for (int y = 0; y < Ny; y++) {
            yp1 = (y + 1) % Ny;
            ym1 = (y - 1 + Ny) % Ny;
            for (int z = 0; z < Nz; z++) {
                zp1 = (z + 1) % Nz;
                zm1 = (z - 1 + Nz) % Nz;

                double beta = gauge3D.p.beta;
                if (gauge3D.p.linearBeta && (z < gauge3D.p.zCenter)) {
                    // vary the coupling linearly to zero away from the center slice
                    beta *= double(z) / double(gauge3D.p.zCenter);
                } else if (gauge3D.p.linearBeta && (z > gauge3D.p.zCenter)) {
                    beta *= double(gauge3D.p.Nz - z - 1) / double(gauge3D.p.zCenter);
                }

                //X dir
                //-------
                f = 0.0;

                // +x, +y, -x, -y
                u1 = gauge3D.read(x,y,z,0);
                u2 = gauge3D.read(xp1,y,z,1);
                u3 = conj(gauge3D.read(x,yp1,z,0));
                u4 = conj(gauge3D.read(x,y,z,1));
                f += beta * imag(u1 * u2 * u3 * u4);

                // -y, +x, +y, -x
                u1 = conj(gauge3D.read(x,ym1,z,1));
                u2 = gauge3D.read(x,ym1,z,0);
                u3 = gauge3D.read(xp1,ym1,z,1);
                u4 = conj(gauge3D.read(x,y,z,0));
                f -= beta * imag(u1 * u2 * u3 * u4);

                if (z != Nz-1) {
                    // +x, +z, -x, -z
                    u1 = gauge3D.read(x,y,z,0);
                    u2 = conj(gauge3D.read(x,y,zp1,0));
                    f += betaZ * imag(u1 * u2);
                }

                if (z != 0) {
                    // -z, +x, +z, -x
                    u1 = gauge3D.read(x,y,zm1,0);
                    u2 = conj(gauge3D.read(x,y,z,0));
                    f -= betaZ * imag(u1 * u2);
                }

                fU.write(x, y, z, 0, f);

                //Y dir
                //------
                f = 0.0;

                // +y, -x, -y, +x
                u1 = gauge3D.read(x,y,z,1);
                u2 = conj(gauge3D.read(xm1,yp1,z,0));
                u3 = conj(gauge3D.read(xm1,y,z,1));
                u4 = gauge3D.read(xm1,y,z,0);
                f += beta * imag(u1 * u2 * u3 * u4);

                // +x, +y, -x, -y
                u1 = gauge3D.read(x,y,z,0);
                u2 = gauge3D.read(xp1,y,z,1);
                u3 = conj(gauge3D.read(x,yp1,z,0));
                u4 = conj(gauge3D.read(x,y,z,1));
                f -= beta * imag(u1 * u2 * u3 * u4);

                if (z != Nz-1) {
                    // y, z, -y, -z
                    u1 = gauge3D.read(x,y,z,1);
                    u2 = conj(gauge3D.read(x,y,zp1,1));
                    f += betaZ * imag(u1 * u2);
                }

                if (z != 0) {
                    // -z, +y, +z, -y
                    u1 = gauge3D.read(x,y,zm1,1);
                    u2 = conj(gauge3D.read(x,y,z,1));
                    f -= betaZ * imag(u1 * u2);
                }

                fU.write(x, y, z, 1, f);
            }
        }
    }
}

void leapfrogHMC::update_mom(double dtau) {

    forceU();
    blas::axpy(-dtau, fU.data, mom3D.data);

    //Update from the 2D fermion field
    if (gauge3D.p.dynamic) {
        int Nx = gauge3D.p.Nx;
        int Ny = gauge3D.p.Ny;
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                for (int mu = 0; mu < 2; mu++) {
                    double temp = mom3D.read(x,y,gauge3D.p.zCenter,mu);
                    temp += fD.read(x,y,mu) * dtau;
                    mom3D.write(x,y,gauge3D.p.zCenter,mu, temp);
                }
            }
        }
    }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void leapfrogHMC::update_mom_gauge(double dtau) {

    forceU();
    blas::axpy(-dtau, fU.data, mom3D.data);
}

void leapfrogHMC::update_mom_fermion(double dtau) {

    //Update from the 2D fermion field
    forceD();
    if (gauge3D.p.dynamic) {
        int Nx = gauge3D.p.Nx;
        int Ny = gauge3D.p.Ny;
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                for (int mu = 0; mu < 2; mu++) {
                    double temp = mom3D.read(x,y,gauge3D.p.zCenter,mu);
                    temp += fD.read(x,y,mu) * dtau;
                    mom3D.write(x,y,gauge3D.p.zCenter,mu, temp);
                }
            }
        }
    }
}


//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void leapfrogHMC::update_gauge(double dtau) {

    for (int i = 0; i < gauge3D.data.size(); i++) {
        gauge3D.data[i] *= polar(1.0, mom3D.data[i] * dtau);
    }
}


// Optimise this to operate only on a single parity of sites.
int leapfrogHMC::forceD() {
    extract2DSlice(gauge2D, gauge3D, gauge3D.p.zCenter);
    int cg_iter = 0;
    if(gauge3D.p.dynamic == true) {

        blas::zero(fD.data);

        //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
        // phip = (D^-1 * Ddag^-1)
        //  phi = (D^-1 * g3 * D^-1 g3) phi.
        cg_iter += cg(phip.data, phi.data, gauge2D, &_DdagDpsiImp);
        //g3Dphi = g3D * phip
        g3Dpsi(g3Dphi, phip, gauge2D);

        double r = 1.0;
        int Nx = gauge2D.p.Nx;
        int Ny = gauge2D.p.Ny;
        // #pragma omp parallel for
        for(int x=0; x<Nx; x++) {
            int xp1 = (x+1)%Nx;
            for(int y=0; y<Ny; y++) {
                int yp1 = (y+1)%Ny;

                Complex gauge_px = gauge2D.read(x,y,0);
                Complex gauge_py = gauge2D.read(x,y,1);
                Complex phip0 = conj(phip.read(x,y,0));
                Complex phip0_px = conj(phip.read(xp1,y,0));
                Complex phip0_py = conj(phip.read(x,yp1,0));
                Complex phip1 = conj(phip.read(x,y,1));
                Complex phip1_px = conj(phip.read(xp1,y,1));
                Complex phip1_py = conj(phip.read(x,yp1,1));
                Complex g3Dphi0 = g3Dphi.read(x,y,0);
                Complex g3Dphi0_px = g3Dphi.read(xp1,y,0);
                Complex g3Dphi0_py = g3Dphi.read(x,yp1,0);
                Complex g3Dphi1 = g3Dphi.read(x,y,1);
                Complex g3Dphi1_px = g3Dphi.read(xp1,y,1);
                Complex g3Dphi1_py = g3Dphi.read(x,yp1,1);

                // antiperiodic boundary conditions in y direction
                if (yp1 == 0) {
                    phip0_py *= -1.0;
                    phip1_py *= -1.0;
                    g3Dphi0_py *= -1.0;
                    g3Dphi1_py *= -1.0;
                }

                double temp = 0.0;
                //mu = 0
                //upper
                // | r  1 |
                // | 1  r |
                //lower
                // | r -1 |
                // | 1 -r |
                temp = real(I*((conj(gauge_px) *
                    (phip0_px * (r*g3Dphi0 +   g3Dphi1) -
                     phip1_px * (  g3Dphi0 + r*g3Dphi1)))
                - (gauge_px *
                    (phip0 * (r*g3Dphi0_px -   g3Dphi1_px) +
                     phip1 * (  g3Dphi0_px - r*g3Dphi1_px)))));

                fD.write(x,y,0,temp);

                //mu = 1
                //upper
                // | r -i |
                // | i  r |
                //lower
                // | r  i |
                // | i -r |
                temp = real(I*((conj(gauge_py) *
                    (phip0_py * (r*g3Dphi0 - I*g3Dphi1) -
                     phip1_py * (I*g3Dphi0 + r*g3Dphi1)))
                - (gauge_py *
                    (phip0 * (r*g3Dphi0_py + I*g3Dphi1_py) +
                     phip1 * (I*g3Dphi0_py - r*g3Dphi1_py)))));

                fD.write(x,y,1,temp);
            }
        }
    }
    return cg_iter;
}

//----------------------------------------------------------------------------------
