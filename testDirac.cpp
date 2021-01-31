// testDirac.cpp

#include <random>
#include "schwinger2peD_internal.h"
#include "utils.h"
#include "blas.h"
#include "dirac_op.h"

using namespace std;

int main(int argc, char **argv) {

    param_t p;
    p.gen = mt19937(1234);

    int sizeX[] = {8,16,32,64,128,256,16,64};
    int sizeY[] = {8,16,32,64,128,256,48,512};
    double mass[] = {0.1,0.3,-0.1,0.0,0.2,0.4,0.2,-0.2};

    for (int i = 0; i < 8; i++) {

        // Lattice size
        p.Nx = sizeX[i]; printf("Nx = %d\n", p.Nx);
        p.Ny = sizeY[i]; printf("Ny = %d\n", p.Ny);
        p.m = mass[i]; printf("m = %lf\n", p.m);

        // generate a random gauge field
        field<Complex> gauge(p); gaussStart(&gauge);

        // generate two random spinor fields
        field<Complex> X(p); gaussComplex(&X);
        field<Complex> Y(p); gaussComplex(&Y);

        // check that < g3DX | Y > == < X | g3DY > to check hermiticity of D
        field<Complex> g3DX(p); g3Dpsi(&g3DX, &X, &gauge);
        field<Complex> g3DY(p); g3Dpsi(&g3DY, &Y, &gauge);

        Complex c1 = blas::cDotProd(g3DX.data, Y.data);
        Complex c2 = blas::cDotProd(X.data, g3DY.data);

        printf("c1 = %.12e + %.12e i\n", real(c1), imag(c1));
        printf("c2 = %.12e + %.12e i\n", real(c2), imag(c2));
        printf("error = %.12e + %.12e i\n", real(c2) - real(c1), imag(c2) - imag(c1));
        printf("\n");
    }

    return 0;
}
