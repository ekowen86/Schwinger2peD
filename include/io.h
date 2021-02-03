#pragma once

#include "schwinger2peD_internal.h"
#include "lattice.h"

void writeGaugeText(field3D<Complex> *gauge, string name);
void readGaugeText(field3D<Complex> *gauge, string name);
void writeGaugeBinary(field3D<Complex>& gauge, string name);
void readGaugeBinary(field3D<Complex>& gauge, string name);
void writeRngState(field3D<Complex>& gauge, string name);
void readRngState(field3D<Complex>& gauge, string name);
