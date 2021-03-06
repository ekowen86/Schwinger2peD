#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <sys/time.h>
#include <random>
#include <algorithm>

using namespace std;

typedef complex<double> Complex;

#define PI 3.141592653589793
#define TWO_PI 6.283185307179586
#define I Complex(0,1.0)
#define cUnit Complex(1.0,0)

using namespace std;

typedef struct {

  //HMC
  int n_step = 25;
  double tau = 1.0;
  int iter_hmc = 1000;
  int therm = 50;
  int skip = 25;
  int chkpt = 100;
  int checkpoint_start = 0;
  int max_iter_cg = 1000;
  double eps = 1e-12;

  int seed = 1234;
  bool verbosity = true;

  //physics
  int Nx = 16;
  int Ny = 16;
  int Nz = 5;
  int zCenter = 2;
  double beta = 3.0;
  double betaZ = 1.0;
  double m = -0.06;
  double musq = 0.0;
  bool dynamic = true;
  bool linearBeta = true;

  //Smearing
  double alpha = 0.5;
  int smear_iter = 1;

  //Wilson loop and Polyakov loop max size.
  int loop_max = 16;

  //Measurements
  bool meas_wl = false; //Wilson loop and Creutz ratios
  bool meas_pc = false; //Pion
  bool meas_vt = false; //Vacuum Trace

  std::mt19937* gen; // random number generator

} param_t;
