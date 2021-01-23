#include "schwinger2peD_internal.h"
#include "utils.h"
#include "hmc.h"
#include "io.h"

int main(int argc, char **argv) {

  struct timeval start, end, total_start, total_end;
  gettimeofday(&total_start, NULL);
  double t_hmc =  0.0;
  double t_meas = 0;
  double t_total =  0.0;
  
  param_t p;

  p.beta = atof(argv[1]);
  p.betaZ = atof(argv[2]); 
  p.iter_hmc = atoi(argv[3]);
  p.therm = atoi(argv[4]);
  p.skip = atoi(argv[5]);
  p.chkpt = atoi(argv[6]);
  p.checkpoint_start = atoi(argv[7]);  
  p.n_step = atoi(argv[8]);
  p.tau = atof(argv[9]);
  
  p.smear_iter = atoi(argv[10]);
  p.alpha = atof(argv[11]);  
  p.seed = atoi(argv[12]);
  p.dynamic = (atoi(argv[13]) == 0 ? false : true);
  p.lockedZ = (atoi(argv[14]) == 0 ? false : true);
  p.m = atof(argv[15]);
  p.max_iter_cg = atoi(argv[16]);
  p.eps = atof(argv[17]);
  
  //Measurements
  p.meas_wl = (atoi(argv[18]) == 0 ? false : true);
  p.meas_pc = (atoi(argv[19]) == 0 ? false : true);
  
  // Lattice size 
  p.Nx = atoi(argv[20]);
  p.Ny = atoi(argv[21]);
  p.Nz = atoi(argv[22]);
  
  if(p.loop_max > std::min(p.Nx/2, p.Ny/2)) {
    cout << "Warning: requested Wilson loop max " << p.loop_max << " greater than ";
    cout << min(p.Nx/2, p.Ny/2) << ", truncating." << endl;
    p.loop_max = std::min(p.Nx/2, p.Ny/2);
  }
  
  //Pseudo RNG seed
  srand48((long)p.seed);
  
  //Topology
  double top = 0.0;
  std::vector<int> top_int(p.Nz, 0);
  std::vector<int> top_old(p.Nz, 0);
  std::vector<int> top_stuck(p.Nz, 0);

  int histL = 101;
  std::vector<std::vector<int>> histQ(p.Nz, std::vector<int>(histL, 0));
  std::vector<double> plaq(p.Nz, 0);
  std::vector<double> plaqSum(p.Nz, 0);
  std::vector<int> index(p.Nz, 0);
      
  int count = 0;
  string name;

  int accept;
  int accepted = 0;
  char fname[256];
  FILE *fp;

  printParams(p);
  
  field<Complex> *gauge2D = new field<Complex>(p);
  field3D<Complex> *gauge3D = new field3D<Complex>(p);
  
  gaussStart(gauge3D);  // hot start
  
  //Start simulation
  double time0 = -((double)clock());
  int iter_offset = 0;
  int iter = 0;
  cout << setprecision(16);

  leapfrogHMC *HMCStep = new leapfrogHMC(p);
  
  if(p.checkpoint_start > 0) {
    
    //Read in gauge field if requested
    //---------------------------------------------------------------------
    name = "gauge/gauge";
    constructName(name, p);
    name += "_traj" + to_string(p.checkpoint_start) + ".dat";	
    readGauge(gauge3D, name);
    iter_offset = p.checkpoint_start;    
  } else {

    //Thermalise from random start
    //---------------------------------------------------------------------

    gettimeofday(&start, NULL);  
    for(iter=0; iter<2*p.therm; iter++){      
      //Perform HMC step
      accept = HMCStep->hmc(gauge3D, iter);
      gettimeofday(&total_end, NULL);  
      t_total = ((total_end.tv_sec  - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
      cout << t_total << " " << endl;                     //Time
    }
    gettimeofday(&end, NULL);  
    t_hmc += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
    iter_offset = 2*p.therm;
  }  
    
  //Begin thermalised trajectories
  //---------------------------------------------------------------------
  for(iter=iter_offset; iter<p.iter_hmc + iter_offset; iter++){

    for(int z=0; z<p.Nz; z++) {
      //Measure the topological charge at each step
      //---------------------------------------------------------------------
      extract2DSlice(gauge2D, gauge3D, z);
      top = measTopCharge(gauge2D);
      top_int[z] = round(top);
      name = "data/top/top_charge_Lz" + to_string(z);
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d %d\n", iter, top_int[z]);
      fclose(fp);
      
      index[z] = top_int[z] + (histL-1)/2;
      histQ[z][index[z]]++;
      if(top_old[z] == top_int[z]) top_stuck[z]++;
      top_old[z] = top_int[z];      
    }
    
    //Perform Measurements
    //---------------------------------------------------------------------
    if((iter)%p.skip == 0 && iter > 2*p.therm) {
      
      count++; //Number of measurements taken
      
      //Plaquette action
      for(int z=0; z<p.Nz; z++) {
	extract2DSlice(gauge2D, gauge3D, z);
	plaq[z] = measPlaq(gauge2D).real();
	plaqSum[z] += plaq[z];
      }
      
      //Dump simulation data to stdout
      gettimeofday(&total_end, NULL);  
      t_total = ((total_end.tv_sec  - total_start.tv_sec) * 1000000u + total_end.tv_usec - total_start.tv_usec) / 1.e6;
      cout << fixed << setprecision(16) << iter << " ";   //Iteration
      cout << t_total << " ";                             //Time
      cout << plaqSum[(p.Nz-1)/2]/count << " ";           //Action
      cout << (double)top_stuck[(p.Nz-1)/2]/(count*p.skip) << " " ;   //P(stuck)
      cout << HMCStep->exp_dH_ave/(count*p.skip) << " ";  //Average exp(-dH)
      cout << HMCStep->dH_ave/(count*p.skip) << " ";      //Average dH
      cout << (double)accepted/(count*p.skip) << " ";     //Acceptance
      for(int z=0; z<p.Nz; z++) printf("%+.2d ", top_int[z]);  //T charge
      cout << endl;
      
      for(int z=0; z<p.Nz; z++) {
	
	//Dump simulation data to file
	//I cannot make bricks without clay!
	name = "data/data/data_Lz" + to_string(z); 
	constructName(name, p);
	name += ".dat";	
	sprintf(fname, "%s", name.c_str());	
	fp = fopen(fname, "a");	
	fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e %.16e %d\n",
		iter,
		t_total,
		plaqSum[z]/count,
		(double)top_stuck[z]/(accepted),
		HMCStep->exp_dH_ave/(count*p.skip),
		HMCStep->dH_ave/(count*p.skip),
		(double)accepted/(count*p.skip),
		top_int[z]);
	fclose(fp);

      
	//Update topoligical charge histogram
	name = "data/top/top_hist_Lz" + to_string(z);
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "w");
	for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[z][i]);
	fclose(fp);
      }

      //Physical observables
      //-------------------------------------------------------------      
      //Polyakov Loops      
      //if(p.meas_pl) measPolyakovLoops(gauge, iter, p);
      
      //Creu<tz Ratios (for string tension)
      if(p.meas_wl) measWilsonLoops(gauge2D, plaq[(p.Nz-1)/2], iter);
      //Pion Correlation
      if(p.meas_pc) measPionCorrelation(gauge2D, iter);
      //Vacuum Trace
      //if(p.meas_vt) measVacuumTrace(gauge, top_old, iter, p);
      //-------------------------------------------------------------
#if 0
      
      //Test deflation routines
      //-------------------------------------------------------------
      if(gauge->p.deflate) {

	// Construct objects for an eigensolver
	//-------------------------------------
	eig_param_t eig_param;
	std::vector<field<Complex>*> kSpace;
	std::vector<Complex> evals;	
	prepareKrylovSpace(kSpace, evals, eig_param, gauge->p);
	
	// Compute a deflation space using IRAM
	iram(gauge, kSpace, evals, eig_param);
		
	// Test the block compression
	//---------------------------
	int Nx = gauge->p.Nx;
	int Ny = gauge->p.Ny;
	int Ns = 2;
	int blk_scheme[2] = {gauge->p.block_scheme[0], gauge->p.block_scheme[1]};
	int x_block_size = Nx/blk_scheme[0];
	int y_block_size = Ny/blk_scheme[1];
	int n_blocks = blk_scheme[0]*blk_scheme[1];
	int blk_size = Ns * x_block_size * y_block_size; // Complex elems per block
	int n_low = gauge->p.n_low;
	int n_conv = eig_param.n_conv;
	int n_deflate = eig_param.n_deflate;
	
	// Object to hold the block orthonormal low mode space
	std::vector<std::vector<Complex>> block_data_ortho(n_blocks, std::vector<Complex> (n_low * blk_size, 0.0));
	// Object to hold the projection coeffiecients of the high modes on the ow space
	std::vector<std::vector<Complex>> block_coef(n_blocks, std::vector<Complex> (n_low * n_conv, 0.0));

	gettimeofday(&start, NULL);

	// Krylov space
	std::vector<field<Complex>*> kSpace_recon(n_conv);
	for(int i=0; i<n_conv; i++) kSpace_recon[i] = new field<Complex>(gauge->p);
	// eigenvalues
	std::vector<Complex> evals_recon(n_conv);
	// Compress kSpace into block_data_ortho and block_coeffs...
	blockCompress(kSpace, block_data_ortho, block_coef, blk_scheme, n_low, n_conv);
	// ...then expand to into kSpace_recon test the quality
	blockExpand(kSpace_recon, block_data_ortho, block_coef, blk_scheme, n_low, n_conv);
	gettimeofday(&end, NULL);  
	
	// Compute the eigenvalues and residua using the reconstructed kSpace
	std::vector<double> resid(n_conv, 0.0);
	computeEvals(gauge, kSpace_recon, resid, evals_recon, n_conv);
	
	cout << "Compare eigenvalues and residua: " << endl;	
	for(int i=0; i<eig_param.n_conv; i++) printf("%d: %e %e \n", i, abs(evals[i].real() - evals_recon[i].real())/evals[i].real(), resid[i]);
	double delta_eval = 0.0;
	double delta_resid = 0.0;  
	for(int i=0; i<n_conv; i++) {
	  delta_eval += abs(evals[i].real() - evals_recon[i].real())/evals[i].real();
	  delta_resid += resid[i];
	}
	printf("<delta eval> = %e\n", delta_eval/n_conv);
	printf("<delta resid> = %e\n", delta_resid/n_conv);
	cout << endl;
	
	// Check compression ratio
	int pre = n_conv * 2 * Nx * Ny;
	int post= n_blocks * n_low * (blk_size + n_conv);
	cout << "Algorithmic compression: " << endl;
	cout << "Complex(double) elems pre = " << pre << " Complex(double) elems post = " << post << endl;
	cout << "Ratio1: " << (100.0 * post)/pre << "% of original data " << endl;
	cout << "Ratio2: " << 100*((1.0 * n_low)/n_conv + (1.0*n_low*n_blocks)/(Ns*Nx*Ny))<< "% of original data " << endl;
	cout << n_low << " low eigenvectors used " << endl;
	cout << (n_conv - n_low) << " high eigenvectors reconstructed " << endl;
	cout << "Compress/decompress time = " << ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 << endl;
	
	// Test deflation with reconstructed space
	//----------------------------------------
	field<Complex> *src = new field<Complex>(gauge->p);
	field<Complex> *sol = new field<Complex>(gauge->p);
	field<Complex> *check = new field<Complex>(gauge->p);

	// Populate src with rands
	gaussComplex(src);

	// Create inverter
	inverterCG *inv = new inverterCG(gauge->p);
	
	// Inversion with no deflation
	int undef_iter = inv->solve(sol, src, gauge);
	
	// Inversion with deflation
	blas::zero(sol->data);
	int def_iter = inv->solve(sol, src, kSpace, evals, gauge);
	DdagDpsi(check, sol, gauge);
	blas::axpy(-1.0, src->data, check->data);
	
	cout << "Deflation efficacy: " << endl;
	cout << "Undeflated CG iter = " << undef_iter << endl;
	cout << "Deflated CG iter   = " << def_iter << endl;
	cout << "Solution fidelity  = " << std::scientific << blas::norm2(check->data) << endl;
      }
      //-------------------------------------------------------------      
#endif
    }
    
    //Perform HMC step
    gettimeofday(&start, NULL);
    accept = HMCStep->hmc(gauge3D, iter);
    accepted += accept;
    gettimeofday(&end, NULL);  
    t_hmc += ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    
    //Checkpoint the gauge field?
    if((iter+1)%p.chkpt == 0) {	  
      name = "gauge/gauge";
      constructName(name, p);
      name += "_traj" + to_string(iter+1) + ".dat";
      writeGauge(gauge3D, name);
    }
  }

  delete gauge2D;
  delete gauge3D;
  delete HMCStep;
  
}
//-------------------------------------------------------------------------------
