#include "hmc.h"
#include "inverters.h"

//2D HMC Routines
//---------------------------------------------------------------------

leapfrogHMC::leapfrogHMC(param_t param){
  inv = new inverterCG(param);
  phip = new field<Complex>(param);
  g3Dphi = new field<Complex>(param);

  mom3D = new field3D<double>(param);
  mom2D = new field<double>(param);
  gauge3D_old = new field3D<Complex>(param);
  gauge2D = new field<Complex>(param);  

  // Fermion fields are 2D
  phi = new field<Complex>(param);
  chi = new field<Complex>(param);
  
};

int leapfrogHMC::hmc(field3D<Complex> *gauge3D, int iter) {

  int accept = 0;
  double H = 0.0, H_old = 0.0;
  
  gauge3D_old->copy(gauge3D);

  // init mom[LX][LY][D]  <mom^2> = 1;
  gaussReal(mom3D); 
  
  if(gauge3D->p.dynamic) {    
    //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
    gaussComplex(chi);
    //Create pseudo fermion field, phi = D * chi
    extract2DSlice(gauge2D, gauge3D, (gauge3D->p.Nz - 1)/2);
    g3Dpsi(phi, chi, gauge2D);
  }

  if (iter >= gauge3D->p.therm) H_old = measAction(mom3D, gauge3D, chi, false);
  trajectory(mom3D, gauge3D, phi, iter);
  if (iter >= gauge3D->p.therm) H = measAction(mom3D, gauge3D, phi, true);

  if (iter >= 2*gauge3D->p.therm) {      
    hmc_count++;
    exp_dH_ave += exp(-(H-H_old));
    dH_ave += (H-H_old);
  }

  // Metropolis accept/reject step
  if (iter >= gauge3D->p.therm) {    
    if ( drand48() > exp(-(H-H_old)) ) gauge3D->copy(gauge3D_old);
    else accept = 1;
  }

  return accept;
}

void leapfrogHMC::trajectory(field3D<double> *mom3D, field3D<Complex> *gauge3D, field<Complex> *phi, int iter){
  
  double dtau = gauge3D->p.tau/gauge3D->p.n_step;
  double H = 0.0;
  field<Complex> *gauge2D = new field<Complex>(gauge3D->p);
  
  double ave_iter = 0;
  
  //gauge force (U field)
  field3D<double> *fU = new field3D<double>(gauge3D->p);
  //fermion force (D operator)
  field<double> *fD = new field<double>(gauge3D->p);

  // Start HMC trajectory
  //----------------------------------------------------------
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU - fD)
  forceU(fU, gauge3D);
  extract2DSlice(gauge2D, gauge3D, (gauge3D->p.Nz - 1)/2);  
  ave_iter += forceD(fD, phi, gauge2D, iter);
  update_mom(fU, fD, mom3D, 0.5*dtau);  

  for(int k=1; k<gauge3D->p.n_step; k++) {


    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge3D, mom3D, dtau);

    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge3D);
    extract2DSlice(gauge2D, gauge3D, (gauge3D->p.Nz - 1)/2);
    ave_iter += forceD(fD, phi, gauge2D, iter);

    update_mom(fU, fD, mom3D, dtau);
  }
  
  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge3D, mom3D, dtau);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
  forceU(fU, gauge3D);
  extract2DSlice(gauge2D, gauge3D, (gauge3D->p.Nz - 1)/2);  
  ave_iter += forceD(fD, phi, gauge2D, iter);
  update_mom(fU, fD, mom3D, 0.5*dtau);

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
  
  int Nx = fU->p.Nx;
  int Ny = fU->p.Ny;
  int Nz = fU->p.Nz;
  double temp = 0.0;

  //Always update from the gauge fields
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int z=0; z<Nz; z++)
	for(int mu=0; mu<2; mu++) {
	  temp = mom->read(x,y,z,mu) - fU->read(x,y,z,mu)*dtau;
	  mom->write(x,y,z,mu, temp);
	}

  //Update from the 2D fermion field
  if(fU->p.dynamic == true) {
    for(int x=0; x<Nx; x++)
      for(int y=0; y<Ny; y++)
	for(int mu=0; mu<2; mu++) {
	  temp = mom->read(x,y,(Nz-1)/2,mu) + fD->read(x,y,mu)*dtau;
	  mom->write(x,y,(Nz-1)/2,mu, temp);
	}
  }
}


//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void leapfrogHMC::update_gauge(field3D<Complex> *gauge, field3D<double> *mom, double dtau){
  
  int Nx = gauge->p.Nx;
  int Ny = gauge->p.Ny;
  int Nz = gauge->p.Nz;
  Complex temp = 0.0;  
  for(int x=0; x<Nx; x++)
    for(int y=0; y<Ny; y++)
      for(int z=0; z<Nz; z++)
	for(int mu=0; mu<2; mu++) {
	  temp = gauge->read(x,y,z,mu) * polar(1.0, mom->read(x,y,z,mu) * dtau);
	  gauge->write(x,y,z,mu, temp);
	  //cout << x << " " << y << " " << z << " " << mu << endl;
	}
}


// Optimise this to operate only on a single parity of sites.
int leapfrogHMC::forceD(field<double> *fD, field<Complex> *phi, field<Complex> *gauge,
			int iter)
{
  int cg_iter = 0;
  if(gauge->p.dynamic == true) {

    blas::zero(fD->data);
    
    //Ainvpsi inverts using the DdagD (g3Dg3D) operator, returns
    // phip = (D^-1 * Ddag^-1)
    //  phi = (D^-1 * g3 * D^-1 g3) phi.
    cg_iter += inv->solve(phip, phi, gauge);
    //g3Dphi = g3D * phip
    g3Dpsi(g3Dphi, phip, gauge);
    
    double r = 1.0;
    int Nx = gauge->p.Nx;
    int Ny = gauge->p.Ny;
#pragma omp parallel for
    for(int x=0; x<Nx; x++) {
      int xp1 = (x+1)%Nx;
      for(int y=0; y<Ny; y++) {
	int yp1 = (y+1)%Ny;

	double temp = 0.0;
	//mu = 0
	//upper
	// | r  1 | 
	// | 1  r |
	//lower
	// | r -1 |
	// | 1 -r |	
	temp = real(I*((conj(gauge->read(x,y,0)) *
			 (conj(phip->read(xp1,y,0)) * (r*g3Dphi->read(x,y,0) +   g3Dphi->read(x,y,1)) -
			  conj(phip->read(xp1,y,1)) * (  g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			-
			(gauge->read(x,y,0) *
			 (conj(phip->read(x,y,0)) * (r*g3Dphi->read(xp1,y,0) -   g3Dphi->read(xp1,y,1)) +
			  conj(phip->read(x,y,1)) * (  g3Dphi->read(xp1,y,0) - r*g3Dphi->read(xp1,y,1))))
			)
		     );
	
	fD->write(x,y,0,temp);
	
	//mu = 1
	//upper
	// | r -i | 
	// | i  r |
	//lower
	// | r  i |
	// | i -r |
	temp = real(I*((conj(gauge->read(x,y,1)) *
			 (conj(phip->read(x,yp1,0)) * (r*g3Dphi->read(x,y,0) - I*g3Dphi->read(x,y,1)) -
			  conj(phip->read(x,yp1,1)) * (I*g3Dphi->read(x,y,0) + r*g3Dphi->read(x,y,1))))
			-			       
			(gauge->read(x,y,1) *
			 (conj(phip->read(x,y,0)) * (r*g3Dphi->read(x,yp1,0) + I*g3Dphi->read(x,yp1,1)) +
			  conj(phip->read(x,y,1)) * (I*g3Dphi->read(x,yp1,0) - r*g3Dphi->read(x,yp1,1))))
			)
		     );
	
	fD->write(x,y,1,temp);
      }
    }
  }
  return cg_iter;
}

//----------------------------------------------------------------------------------

leapfrogHMC::~leapfrogHMC() {
  
  delete inv;
  delete phip;
  delete g3Dphi;

  delete mom3D;
  delete mom2D;
  delete gauge3D_old;
  delete gauge2D;

  delete phi;
  delete chi;
  
};
