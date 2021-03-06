#include "measurements.h"
#include "utils.h"
#include "dirac_op.h"
#include "inverters.h"

Complex measWilsonLoop(field<Complex>& gauge, int a, int b) {
    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;

    Complex wSum = Complex(0.0, 0.0);

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {

            Complex w = Complex(1.0, 0.0);
            int x1 = x;
            int y1 = y;

            //Move in +x up to a
            for (int n = 0; n < a; n++) {
                w *= gauge.read(x1,y1,0);
                x1 = (x1 + 1) % Nx;
            }

            //Move in +y up to b
            for (int n = 0; n < b; n++) {
                w *= gauge.read(x1,y1,1);
                y1 = (y1 + 1) % Ny;
            }

            //Move in -x from a to 0
            for (int n = 0; n < a; n++) {
                x1 = (x1 - 1 + Nx) % Nx;
                w *= conj(gauge.read(x1,y1,0));
            }

            //Move in -y from b to 0
            for (int n = 0; n < b; n++) {
                y1 = (y1 - 1 + Ny) % Ny;
                w *= conj(gauge.read(x1,y1,1));
            }
            wSum += w;
        }
    }
    return wSum / double(Nx * Ny);
}

// void measWilsonLoops(field<Complex>& gauge, double plaq, int iter) {
//   int Nx = gauge.p.Nx;
//   int Ny = gauge.p.Ny;
//
//   std::vector<std::vector<Complex>> wLoops(Nx/2 ,std::vector<Complex> (Ny/2 ,0.0));
//   std::vector<double> sigma(Nx/2, 0.0);
//
//   Complex w;
//   int p1, p2, dx, dy, x, y;
//   double inv_Lsq = 1.0/(Nx*Ny);
//   int loop_max = gauge.p.loop_max;
//
//   //Smear the gauge field
//   field<Complex> smeared(gauge.p);
//   smearLink(smeared, gauge);
//
//   //Loop over all X side sizes of rectangle
// // #pragma omp parallel for
//   for(int Xrect=1; Xrect<loop_max; Xrect++) {
//
//     //Loop over all Y side sizes of rectangle
//     for(int Yrect=1; Yrect<loop_max; Yrect++) {
//
//       //Loop over all x,y starting points
//       for(x=0; x<Nx; x++)
//     for(y=0; y<Ny; y++){
//
//       Complex w = Complex(1.0,0.0);
//
//       //Move in +x up to p1.
//       for(int dx=0; dx<Xrect; dx++)     w *= smeared.read(x+dx,y,0);
//
//       //Move in +y up to p2 (p1 constant)
//       int p1 = (x + Xrect)%Nx;
//       for(int dy=0; dy<Yrect; dy++)     w *= smeared.read(p1,y+dy,1);
//
//       //Move in -x from p1 to (p2 constant)
//       int p2 = (y + Yrect)%Ny;
//       for(int dx=Xrect-1; dx>=0; dx--)  w *= conj(smeared.read(x+dx,p2,0));
//
//       //Move in -y from p2 to y
//       for(int dy=Yrect-1; dy>=0; dy--)  w *= conj(smeared.read(x,y+dy,1));
//       wLoops[Xrect][Yrect] += w*inv_Lsq;
//     }
//     }
//   }
//
//   //Compute string tension
// // #pragma	omp parallel for
//   for(int size=1; size<loop_max; size++) {
//     sigma[size]  = -log(abs((real(wLoops[size][size])/real(wLoops[size-1][size]))*
//         (real(wLoops[size-1][size-1])/real(wLoops[size][size-1]))));
//
//     sigma[size] += -log(abs((real(wLoops[size][size])/real(wLoops[size][size-1]))*
//         (real(wLoops[size-1][size-1])/real(wLoops[size-1][size]))));
//
//     sigma[size] *= 0.5;
//
//   }
//
//   string name;
//   FILE *fp;
//
//   name = "data/creutz/creutz";
//   constructName(name, gauge.p);
//   name += ".dat";
//   fp = fopen(name.c_str(), "a");
//   fprintf(fp, "%d %.16e ", iter, -log(abs(plaq)) );
//   for(int size=2; size<loop_max; size++)
//     fprintf(fp, "%.16e ", sigma[size]);
//   fprintf(fp, "\n");
//   fclose(fp);
//
//   for(int sizex=2; sizex<loop_max; sizex++)
//     for(int sizey=sizex-1; (sizey < loop_max && sizey <= sizex+1); sizey++) {
//       name = "data/rect/rectWL";
//       name += "_" + to_string(sizex) + "_" + to_string(sizey);
//       constructName(name, gauge.p);
//       name += ".dat";
//       fp = fopen(name.c_str(), "a");
//       fprintf(fp, "%d %.16e %.16e\n", iter, real(wLoops[sizex][sizey]), imag(wLoops[sizex][sizey]));
//       fclose(fp);
//     }
//
//   return;
// }

//Pion correlation function
//                              |----------------|
//                              |        |-------|---------|
//  < pi(x) | pi(0) > = < ReTr[dn*(x) g3 up(x) | dn(0) g3 up*(0)] >
//                    = < ReTr(g3 Gd[0,x] g3 Gu[x,0]) >
//                    = < ReTr(G*[x,0] G[x,0]) >
//
// using g3 G[x,0] g3 = G*[x,0] and Gu[x,0] \eq Gd[x,0]
//
// if H = Hdag, Tr(H * Hdag) = Sum_{n,m} (H_{n,m}) * (H_{n,m})^*,
// i.e., the sum of the modulus squared of each element

//void measWilsonLoops(field<Complex> *gauge, double plaq, int iter){
void measPionCorrelation(field<Complex>& gauge, int iter) {
    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;

    field<Complex> propUp(gauge.p);
    field<Complex> propDn(gauge.p);
    field<Complex> source(gauge.p);
    field<Complex> Dsource(gauge.p);

    //Up type prop
    blas::zero(source.data);
    blas::zero(Dsource.data);
    blas::zero(propUp.data);
    source.write(0, 0, 0, cUnit); // spinor index is 0

    // up -> (g3Dg3) * up *****
    // (g3Dg3D)^-1 * (g3Dg3) up = D^-1 * up *****

    g3psi(Dsource, source);
    g3Dpsi(source, Dsource, gauge);

    cg(propUp.data, source.data, gauge, &_DdagDpsi);

    //Down type prop
    blas::zero(source.data);
    blas::zero(Dsource.data);
    blas::zero(propDn.data);
    source.write(0, 0, 1, cUnit); // spinor index is 1

    // dn -> (g3Dg3) * dn
    g3psi(Dsource, source);
    g3Dpsi(source, Dsource, gauge);

    // (g3Dg3D)^-1 * (g3Dg3) dn = D^-1 * dn
    cg(propDn.data, source.data, gauge, &_DdagDpsi);

    double pion_corr[Ny + 1];
    // Let y be the 'time' dimension
    // double corr = 0.0, tmp = 0.0;
    for (int y = 0; y < Ny; y++) {
        //initialise
        pion_corr[y] = 0.0;
        //Loop over space and spin, fold propagator
        // corr = 0.0;
        for (int x = 0; x < Nx; x++) {
            pion_corr[y] += norm(propDn.read(x,y,0));
            pion_corr[y] += norm(propDn.read(x,y,1));
            pion_corr[y] += norm(propUp.read(x,y,0));
            pion_corr[y] += norm(propUp.read(x,y,1));
        }
    }

    // use boundary conditions to get correlator at Ny (simplifies folding)
    pion_corr[Ny] = pion_corr[0];

    FILE *fp = fopen("pion_corr.dat", "a");
    fprintf(fp, "%06d", iter);
    // fold correlators here
    for(int y = 0; y <= (Ny / 2); y++) {
        fprintf(fp, " %.16e", (pion_corr[y] + pion_corr[Ny - y]) * 0.5);
    }
    fprintf(fp, "\n");
    fclose(fp);
}

double measGaugeAction(field3D<Complex>& gauge) {

    double betaZ = gauge.p.betaZ;
    double action = 0.0;

    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;
    int Nz = gauge.p.Nz;
    Complex u1, u2, u3, u4;
    for (int x = 0; x < Nx; x++) {
        int xp1 = (x + 1) % Nx;
        for (int y = 0; y < Ny; y++) {
            int yp1 = (y + 1) % Ny;
            for (int z = 0; z < Nz; z++) {
                int zp1 = (z + 1) % Nz;

                double beta = gauge.p.beta;
                if (gauge.p.linearBeta && (z < gauge.p.zCenter)) {
                    beta *= double(z) / double(gauge.p.zCenter);
                } else if (gauge.p.linearBeta && (z > gauge.p.zCenter)) {
                    beta *= double(gauge.p.Nz - z - 1) / double(gauge.p.zCenter);
                }

                u1 = gauge.read(x,y,z,0);
                u2 = gauge.read(xp1,y,z,1);
                u3 = conj(gauge.read(x,yp1,z,0));
                u4 = conj(gauge.read(x,y,z,1));

                double plaq = real(u1 * u2 * u3 * u4);
                action += beta * (1.0 - plaq);

                //Compute extra dim contribution
                if(z != Nz - 1) {
                    //+x, +z, -x, -z
                    u1 = gauge.read(x,y,z,0);
                    u2 = conj(gauge.read(x,y,zp1,0));
                    plaq = real(u1 * u2);
                    action += betaZ * (1.0 - plaq);

                    //+y, +z, -y, -z
                    u1 = gauge.read(x,y,z,1);
                    u2 = conj(gauge.read(x,y,zp1,1));
                    plaq = real(u1 * u2);
                    action += betaZ * real(1.0 - plaq);
                }
            }
        }
    }
    return action;
}

Complex measPlaq(field<Complex>& gauge) {
    Complex plaq = 0.0;
    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;

// #pragma omp parallel for reduction(+:plaq)
    Complex u1, u2, u3, u4;
    for (int x = 0; x < Nx; x++) {
        int xp1 = (x + 1) % Nx;
        for (int y = 0; y < Ny; y++) {
            int yp1 = (y + 1) % Ny;
            // Anti-clockwise plaquette, starting at (x,y)
            u1 = gauge.read(x, y, 0);
            u2 = gauge.read(xp1, y, 1);
            u3 = conj(gauge.read(x, yp1, 0));
            u4 = conj(gauge.read(x, y, 1));
            plaq += u1 * u2 * u3 * u4;
        }
    }
    return plaq / double(Nx * Ny);
}

double measFieldStrength(field<Complex>& gauge) {
    double E = 0.0;
    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;

    Complex u1, u2, u3, u4;
    double F;
    for (int x = 0; x < Nx; x++) {
        int xp1 = (x + 1) % Nx;
        int xm1 = (x - 1 + Nx) % Nx;
        for (int y = 0; y < Ny; y++) {
            int yp1 = (y + 1) % Ny;
            int ym1 = (y - 1 + Ny) % Ny;

            u1 = gauge.read(x, y, 1);
            u1 *= gauge.read(x, yp1, 0);
            u1 *= conj(gauge.read(xp1, y, 1));
            u1 *= conj(gauge.read(x, y, 0));

            u2 = gauge.read(x, y, 0);
            u2 *= conj(gauge.read(xp1, ym1, 1));
            u2 *= conj(gauge.read(x, ym1, 0));
            u2 *= gauge.read(x, ym1, 1);

            u3 = conj(gauge.read(x, ym1, 1));
            u3 *= conj(gauge.read(xm1, ym1, 0));
            u3 *= gauge.read(xm1, ym1, 1);
            u3 *= gauge.read(xm1, y, 0);

            u4 = conj(gauge.read(xm1, y, 0));
            u4 *= gauge.read(xm1, y, 1);
            u4 *= gauge.read(xm1, yp1, 0);
            u4 *= conj(gauge.read(x, y, 1));

            F = imag(u1 + u2 + u3 + u4) * 0.25;
            E += (F * F);
        }
    }
    return E / double(Nx * Ny) * 0.5;
}

void drawInstantons(field<Complex>& gauge) {

    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;
    Complex u1, u2, u3, u4;
    double a1, a2, a3, a4;
    printf("\n\n\n\n\n\n\n");
    int q = 0;
    for (int x = 0; x < Nx; x++) {
        int xp1 = (x + 1) % Nx;
        for (int y = 0; y < Ny; y++) {
            int yp1 = (y + 1) % Ny;
            u1 = gauge.read(x, y, 0);
            u2 = gauge.read(xp1, y, 1);
            u3 = conj(gauge.read(x, yp1, 0));
            u4 = conj(gauge.read(x, y, 1));

            double w = arg(u1) + arg(u2) + arg(u3) + arg(u4);
            if (w <= -PI) {
                printf(" +");
                q++;
            } else if (w >= PI) {
                printf(" -");
                q--;
            } else {
                printf("  ");
            }
        }
        printf("\n");
    }
    printf("%+d\n", q);
}

double measTopCharge(field<Complex>& gauge) {

    double top = 0.0;
    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;
// #pragma omp parallel for reduction(+:top)
    Complex u1, u2, u3, u4;
    for (int x = 0; x < Nx; x++) {
        int xp1 = (x + 1) % Nx;
        for (int y = 0; y < Ny; y++) {
            int yp1 = (y + 1) % Ny;
            u1 = gauge.read(x, y, 0);
            u2 = gauge.read(xp1, y, 1);
            u3 = conj(gauge.read(x, yp1, 0));
            u4 = conj(gauge.read(x, y, 1));
            Complex w = u1 * u2 * u3 * u4;
            top += arg(w);  // -pi < arg(w) < pi  Geometric value is an integer.
        }
    }

    return top / TWO_PI;
}

void measChiralCond(field<Complex>& gauge, int iter) {

    int Nx = gauge.p.Nx;
    int Ny = gauge.p.Ny;

    field<Complex> propUp(gauge.p);
    field<Complex> propDn(gauge.p);
    field<Complex> source(gauge.p);
    field<Complex> Dsource(gauge.p);

    double cc = 0.0;

    vector<int> vPoints;

    int nPoints = max(Nx, Ny);

    for (int i = 0; i < nPoints; i++) {

        int n = 0;
        do {
            n = gauge.rand_int(0, Nx * Ny - 1);
        } while (find(vPoints.begin(), vPoints.end(), n) != vPoints.end());
        vPoints.push_back(n);

        int x = n % Nx;
        int y = n / Nx;

        // Up type prop
        blas::zero(source.data);
        blas::zero(Dsource.data);
        blas::zero(propUp.data);
        source.write(x, y, 0, cUnit);

        // up -> (g3Dg3) * up *****
        // (g3Dg3D)^-1 * (g3Dg3) up = D^-1 * up *****

        g3psi(Dsource, source);
        g3Dpsi(source, Dsource, gauge);

        cg(propUp.data, source.data, gauge, &_DdagDpsi);

        // Down type prop
        blas::zero(source.data);
        blas::zero(Dsource.data);
        blas::zero(propDn.data);
        source.write(x, y, 1, cUnit);

        // dn -> (g3Dg3) * dn
        g3psi(Dsource, source);
        g3Dpsi(source, Dsource, gauge);

        // (g3Dg3D)^-1 * (g3Dg3) dn = D^-1 * dn
        cg(propDn.data, source.data, gauge, &_DdagDpsi);

        cc += real(propUp.read(x,y,0) + propDn.read(x,y,1));
    }

    // for (int x = 0; x < Nx; x++) {
    //     for (int y = 0; y < Ny; y++) {
    //         //Up type prop
    //         blas::zero(source.data);
    //         blas::zero(Dsource.data);
    //         blas::zero(propUp.data);
    //         source.write(x, y, 0, cUnit);
    //
    //         // up -> (g3Dg3) * up *****
    //         // (g3Dg3D)^-1 * (g3Dg3) up = D^-1 * up *****
    //
    //         g3psi(Dsource, source);
    //         g3Dpsi(source, Dsource, gauge);
    //
    //         cg(propUp.data, source.data, gauge, &_DdagDpsi);
    //
    //         //Down type prop
    //         blas::zero(source.data);
    //         blas::zero(Dsource.data);
    //         blas::zero(propDn.data);
    //         source.write(x, y, 1, cUnit);
    //
    //         // dn -> (g3Dg3) * dn
    //         g3psi(Dsource, source);
    //         g3Dpsi(source, Dsource, gauge);
    //
    //         // (g3Dg3D)^-1 * (g3Dg3) dn = D^-1 * dn
    //         cg(propDn.data, source.data, gauge, &_DdagDpsi);
    //
    //         cc += real(propUp.read(x,y,0) + propDn.read(x,y,1));
    //     }
    // }

    FILE *fp = fopen("cc.dat", "a");
    fprintf(fp, "%06d", iter);
    fprintf(fp, " %.16e", cc / double(nPoints));
    fprintf(fp, "\n");
    fclose(fp);
}
