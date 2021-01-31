#include "inverters.h"
#include "dirac_op.h"

int cg(cvec& P, const cvec& Q, const field<Complex>& gauge, op f) {

    int success = 0;
    bool verbose = false;
    int size = Q.size();

    cvec r(size); // residue
    cvec p(size); // search vector
    cvec x(size); //
    cvec q2p(size); //

    double normsq, pro, err, alphaCG, betaCG;

    blas::zero(x);

    // compute initial residual
    blas::copy(r, Q);
    blas::copy(p, r);
    normsq = blas::norm2(Q);

    if(normsq == 0 || !isfinite(normsq)) {
        cout << "Error: CG source is zero or nan." << endl;
        exit(0);
    }

    double maxRelErr = gauge.p.eps * normsq;

    // Iterate until convergence
    int k;
    for (k = 0; k < gauge.p.max_iter_cg; k++) {

        // Compute f(p).
        f(q2p, p, gauge);

        pro = real(blas::cDotProd(p, q2p));

        alphaCG = normsq / pro;

        blas::axpy( alphaCG, p, x);
        blas::axpy(-alphaCG, q2p, r);

        // Exit if new residual is small enough
        err = blas::norm2(r);
        if (verbose) printf("CG iter %d, err = %g\n", k+1, err);
        if (err <= maxRelErr) {
            blas::copy(P, x);
            normsq = err;
            break;
        }

        // Update vec using new residual
        betaCG = err / normsq;
        blas::axpy(betaCG, p, r, p);
        normsq = err;

    } // End loop over k

    if (k == gauge.p.max_iter_cg) {
        // Failed convergence
        printf("CG: Failed to converge iter = %d, err = %.16e\n", k+1, normsq);
        success = 0;
    } else {
        // Convergence
        success = k+1;
    }

    if(verbose) {

        // Sanity
        cout << "source norm = " << blas::norm(Q) << endl;
        cout << "sol norm = " << blas::norm(P) << endl;

        cvec temp(size);
        f(temp, x, gauge);
        blas::axpy(-1.0, temp, Q, r);
        double truersq = blas::norm2(r);
        printf("CG: Converged iter = %d, err = %.16e, truersq = %.16e\n", k+1, normsq, truersq / blas::norm(Q));
    }

    return success;
}
