#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"
#include "blas.h"
#include "hdf5.h"

void printParams(param_t p);
void constructName(string &name, param_t p);
void writeGauge(field3D<Complex> *gauge, string name);
void readGauge(field3D<Complex> *gauge, string name);

void gaussStart(field3D<Complex> *gauge);
void coldStart(field3D<Complex> *gauge);

void gaussStart(field<Complex> *gauge);
void coldStart(field<Complex> *gauge);

// Normalized gaussian exp(-phi*phi/2) and  <phi|phi> = 1
void gaussReal(field<double> *field);
void gaussReal(field3D<double> *field);

//normalized gaussian exp[ - eta*eta/2]  <eta|eta> = 1;
void gaussComplex(field<Complex> *field);

//APE smearing: project back on U(1)
// staple x is 0th, y is 1st.
void smearLink(field<Complex> *smeared, field<Complex> *gauge);

// Extract 2D lattice slice
void extract2DSlice(field<Complex> *gauge,
		    field3D<Complex> *gauge3D,
		    int slice);



void measBlockColinearity(std::vector<field<Complex> *> kSpace, int blockScheme[2], int nLow);

void blockCompress(std::vector<field<Complex> *> &kSpace,
		   std::vector<std::vector<Complex>> &block_data_ortho,
		   std::vector<std::vector<Complex>> &block_coef,
		   int blockScheme[2], int n_low, int n_conv);

void blockExpand(std::vector<field<Complex> *> &kSpace,
		 std::vector<std::vector<Complex>> &block_data_ortho,
		 std::vector<std::vector<Complex>> &block_coef,
		 int blockScheme[2], int n_low, int n_conv);

void readVectorToBlock(std::vector<field<Complex> *> &kSpace,
		       std::vector<std::vector<Complex>> &block_data,
		       int blockScheme[2], int iEig);

void hdf5Example();
