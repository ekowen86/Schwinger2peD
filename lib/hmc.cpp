#include "inverters.h"
#include "dirac_op.h"
#include "utils.h"
#include "measurements.h"
#include "hmc.h"

//2D HMC Routines
//---------------------------------------------------------------------

leapfrogHMC::leapfrogHMC(param_t param) : gauge3D_hmc(param), gauge2D(param) {
  phip = new field<Complex>(param);
  g3Dphi = new field<Complex>(param);

  mom3D = new field3D<double>(param);
  mom2D = new field<double>(param);

  // Fermion fields are 2D
  phi = new field<Complex>(param);
  chi = new field<Complex>(param);

};

leapfrogHMC::~leapfrogHMC() {

  delete phip;
  delete g3Dphi;

  delete mom3D;
  delete mom2D;

  delete phi;
  delete chi;

};

int leapfrogHMC::hmc(field3D<Complex> *gauge3D, bool noMetropolis) {

    int accept = 0;

    gauge3D_hmc.copy(gauge3D);
    gaussReal(mom3D);

    if (gauge3D_hmc.p.dynamic) {
        //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
        gaussComplex(chi);
        //Create pseudo fermion field, phi = D * chi
        extract2DSlice(&gauge2D, &gauge3D_hmc, gauge3D_hmc.p.zCenter);
        g3Dpsi(phi, chi, &gauge2D);
        blas::caxpy(-I * sqrt(gauge3D->p.musq), chi->data, phi->data);
    }

    double oldH = 0.0;
    oldH += blas::norm2(mom3D->data) * 0.5;
    oldH += measGaugeAction(&gauge3D_hmc);
    if (gauge3D->p.dynamic) {
        oldH += real(blas::cDotProd(chi->data, chi->data));
    }

    trajectory(mom3D, &gauge3D_hmc, phi);

    double newH = 0.0;
    newH += blas::norm2(mom3D->data) * 0.5;
    newH += measGaugeAction(&gauge3D_hmc);
    if (gauge3D->p.dynamic) {
        extract2DSlice(&gauge2D, &gauge3D_hmc, gauge3D_hmc.p.zCenter);
        cg(chi->data, phi->data, gauge2D, &_DdagDpsiImp);
        newH += real(blas::cDotProd(chi->data, phi->data));
    }

    dH = newH - oldH;
    exp_dH = exp(-dH);

    if (dH <= 0 || gauge3D_hmc.rand_double() < exp_dH) accept = 1;
    if (accept || noMetropolis) gauge3D->copy(&gauge3D_hmc);

    return accept;
}

void leapfrogHMC::trajectory(field3D<double> *mom3D, field3D<Complex> *gauge3D, field<Complex> *phi){

  // time step
  double dtau = gauge3D->p.tau / gauge3D->p.n_step;

  // gauge force (U field)
  field3D<double> fU(gauge3D->p);

  // fermion force (D operator)
  field<double> fD(gauge3D->p);

  // Start HMC trajectory
  //----------------------------------------------------------
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU - fD)
  forceU(&fU, gauge3D);
  extract2DSlice(&gauge2D, gauge3D, gauge3D->p.zCenter);
  // ave_iter += forceD(fD, phi, gauge2D);
  forceD(&fD, phi, &gauge2D);
  update_mom(&fU, &fD, mom3D, 0.5*dtau);

  // intermediate steps
  for (int k = 1; k < gauge3D->p.n_step; k++) {

    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge3D, mom3D, dtau);

    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(&fU, gauge3D);
    extract2DSlice(&gauge2D, gauge3D, gauge3D->p.zCenter);
    forceD(&fD, phi, &gauge2D);

    update_mom(&fU, &fD, mom3D, dtau);
  }

  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge3D, mom3D, dtau);

  //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
  forceU(&fU, gauge3D);
  extract2DSlice(&gauge2D, gauge3D, gauge3D->p.zCenter);
  forceD(&fD, phi, &gauge2D);
  update_mom(&fU, &fD, mom3D, 0.5*dtau);

  // HMC trajectory complete
  //----------------------------------------------------------
}

void leapfrogHMC::forceU(field3D<double> *fU, const field3D<Complex> *gauge3D) {

  Complex plaq = 0.0, plaq0 = 0.0;
  double temp = 0.0;
  int xp1, xm1, yp1, ym1, zp1, zm1;

  int Nx = gauge3D->p.Nx;
  int Ny = gauge3D->p.Ny;
  int Nz = gauge3D->p.Nz;
  double beta = gauge3D->p.beta;
  double betaz = gauge3D->p.betaZ;

  for(int x=0; x<Nx; x++) {
    xp1 = (x+1)%Nx;
    xm1 = (x-1+Nx)%Nx;
    for(int y=0; y<Ny; y++) {
      yp1 = (y+1)%Ny;
      ym1 = (y-1+Ny)%Ny;
      for(int z=0; z<Nz; z++) {
	zp1 = (z+1)%Nz;
	zm1 = (z-1+Nz)%Nz;

	//X dir
	//-------
	// +x, +y, -x, -y
	plaq0 = gauge3D->read(x,y,z,0) * gauge3D->read(xp1,y,z,1) *
	  conj(gauge3D->read(x,yp1,z,0)) * conj(gauge3D->read(x,y,z,1));

	temp = beta*(plaq0.imag());

	// -y, +x, +y, -x
	plaq = conj(gauge3D->read(x,ym1,z,1)) * gauge3D->read(x,ym1,z,0) *
	  gauge3D->read(xp1,ym1,z,1) * conj(gauge3D->read(x,y,z,0));

	temp -= beta*(plaq.imag());

	if(z != Nz-1) {
	  // +x, +z, -x, -z
	  plaq = gauge3D->read(x,y,z,0) * cUnit * conj(gauge3D->read(x,y,zp1,0)) * cUnit;
	  temp += betaz*(plaq.imag());
	}

	if(z != 0) {
	  // -z, +x, +z, -x
	  plaq = cUnit * gauge3D->read(x,y,zm1,0) * cUnit * conj(gauge3D->read(x,y,z,0));
	  temp -= betaz*(plaq.imag());
	}

	fU->write(x,y,z,0, temp);

	temp = 0.0;

	//Y dir
	//------
	// +y, -x, -y, +x
	plaq = gauge3D->read(x,y,z,1) * conj(gauge3D->read(xm1,yp1,z,0)) * conj(gauge3D->read(xm1,y,z,1)) * gauge3D->read(xm1,y,z,0);

	temp += beta*(plaq.imag());

	//This plaquette was aleady computed. We want the conjugate.
	temp -= beta*(plaq0.imag());

	if(z != Nz-1) {
	  // y, z, -y, -z
	  plaq = gauge3D->read(x,y,z,1) * cUnit * conj(gauge3D->read(x,y,zp1,1)) * cUnit;
	  temp += betaz*(plaq.imag());
	}

	if(z != 0) {
	  // -z, +y, +z, -y
	  plaq = cUnit * gauge3D->read(x,y,zm1,1) * cUnit * conj(gauge3D->read(x,y,z,1));
	  temp -= betaz*(plaq.imag());
	}

	fU->write(x,y,z,1, temp);
      }
    }
  }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void leapfrogHMC::update_mom(field3D<double> *fU, field<double> *fD, field3D<double> *mom, double dtau){

    //Always update from the gauge fields
    blas::axpy(-dtau, fU->data, mom->data);

    //Update from the 2D fermion field
    if (fU->p.dynamic) {
        for (int x = 0; x < fU->p.Nx; x++) {
            for (int y = 0; y < fU->p.Ny; y++) {
                for (int mu = 0; mu < 2; mu++) {
                    double temp = mom->read(x,y,fU->p.zCenter,mu);
                    temp += fD->read(x,y,mu) * dtau;
                    mom->write(x,y,fU->p.zCenter,mu, temp);
                }
            }
        }
    }
}


//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void leapfrogHMC::update_gauge(field3D<Complex> *gauge, field3D<double> *mom, double dtau){

    for (int x = 0; x < gauge->p.Nx; x++) {
        for (int y = 0; y < gauge->p.Ny; y++) {
            for (int z = 0; z < gauge->p.Nz; z++) {
                for (int mu = 0; mu < 2; mu++) {
                    Complex temp = gauge->read(x,y,z,mu);
                    temp *= polar(1.0, mom->read(x,y,z,mu) * dtau);
                    gauge->write(x,y,z,mu, temp);
                }
            }
        }
    }
}


// Optimise this to operate only on a single parity of sites.
int leapfrogHMC::forceD(field<double> *fD, field<Complex> *phi, field<Complex> *gauge) {
    int cg_iter = 0;
    if(gauge->p.dynamic == true) {

        blas::zero(fD->data);

        //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
        // phip = (D^-1 * Ddag^-1)
        //  phi = (D^-1 * g3 * D^-1 g3) phi.
        cg_iter += cg(phip->data, phi->data, *gauge, &_DdagDpsiImp);
        //g3Dphi = g3D * phip
        g3Dpsi(g3Dphi, phip, gauge);

        double r = 1.0;
        int Nx = gauge->p.Nx;
        int Ny = gauge->p.Ny;
        // #pragma omp parallel for
        for(int x=0; x<Nx; x++) {
            int xp1 = (x+1)%Nx;
            for(int y=0; y<Ny; y++) {
                int yp1 = (y+1)%Ny;

                Complex gauge_px = gauge->read(x,y,0);
                Complex gauge_py = gauge->read(x,y,1);
                Complex phip0 = conj(phip->read(x,y,0));
                Complex phip0_px = conj(phip->read(xp1,y,0));
                Complex phip0_py = conj(phip->read(x,yp1,0));
                Complex phip1 = conj(phip->read(x,y,1));
                Complex phip1_px = conj(phip->read(xp1,y,1));
                Complex phip1_py = conj(phip->read(x,yp1,1));
                Complex g3Dphi0 = g3Dphi->read(x,y,0);
                Complex g3Dphi0_px = g3Dphi->read(xp1,y,0);
                Complex g3Dphi0_py = g3Dphi->read(x,yp1,0);
                Complex g3Dphi1 = g3Dphi->read(x,y,1);
                Complex g3Dphi1_px = g3Dphi->read(xp1,y,1);
                Complex g3Dphi1_py = g3Dphi->read(x,yp1,1);

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

                fD->write(x,y,0,temp);

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

                fD->write(x,y,1,temp);
            }
        }
    }
    return cg_iter;
}

//----------------------------------------------------------------------------------
