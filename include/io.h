#pragma once

#include "schwinger2peD_internal.h"
#include "blas.h"

using namespace std;

void writeGauge(field<Complex> *gauge, string name);
void readGauge(field<Complex> *gauge, string name);
