/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/14/2017
    File:   CPD.chpl

    Description:    This is a module file for things related to
                    performing the CPD.
*/

module CPD {
    use Base;
    use CSF;
    use Args;
    use Kruskal;
    use Matrices;
    use MTTKRP;
    use IO.FormattedIO;
    use Time;
    use ThreadInfo;
    use Barriers;
    use BLAS;

   /*****************************
    *
    *   Public Functions
    *
    ******************************/
    /*########################################################################
    #   Descriptipn:    Performs the CPD.
    #
    #   Parameters:     tensors (splatt_csf[]):     Tensor to factor in CSF
    #                   args (cpd_cmd_args):        Arguments
    #                   factored (splatt_kruskal):  Stores the result of CPD
    #
    #   Return:         int: SPLATT_SUCCESS on success
    ########################################################################*/
    proc splatt_cpd_als(const tensors, const args : cpd_cmd_args, factored : splatt_kruskal) : int
    {
        const nmodes : int = tensors[0].nmodes;
        const nfactors : int = args.decompRank;
        // Create factor matrices. One for each mode, plus 1 extra
        // for scratch space
        var mats : [0..(nmodes+1)-1] dense_matrix;
        for i in 0..(nmodes+1)-1 {
            mats[i] = new dense_matrix();
        }

        // Determine the maximum dimension of the tensor
        var maxdim : int = 0;
        for i in 0..nmodes-1 {
            maxdim = max(maxdim, tensors[0].dims[i]);
        }

        // Resize factor matrix domains and generate random values for
        // their data
        for m in 0..nmodes-1 {
            mats[m].matrix_domain = {0..tensors[0].dims[m]-1, 0..nfactors-1};
            mats[m].I = tensors[0].dims[m];
            mats[m].J = nfactors;
            mat_rand(mats[m], m);
            mats[m].vals_ref = c_ptrTo(mats[m].vals);
        }

        // Last matrix has maxdim rows. We store result of MTTKRP here
        mats[nmodes].matrix_domain = {0..maxdim, 0..nfactors-1};
        mats[nmodes].I = maxdim;
        mats[nmodes].J = nfactors;
        mats[nmodes].vals_ref = c_ptrTo(mats[nmodes].vals);

        // Perform iterations
        factored.fit = cpd_als_iterate(tensors, mats, factored.lambda_vals, nfactors, args);

        writef("Final fit: %0.5dr\n", factored.fit);

        // Store output
        factored.rank = nfactors;
        factored.nmodes = nmodes;
        factored.dims = tensors[0].dims;
        for m in 0..nmodes-1 {
            factored.factors[m] = mats[m];
        }

        return SPLATT_SUCCESS;
    }

    /*########################################################################
    #   Descriptipn:    Performs the iterations of CP-ALS
    #
    #   Parameters:     tensors (splatt_csf[]):     Tensor to factor in CSF
    #                   mats (dense_matrix[]):      Factor matrices
    #                   lambda_vals (real[]):       Output vector for scaling
    #                   nfactors (int):             Rank of factorization
    #                   args (cpd_cmd_args):        Arguments
    #
    #   Return:         real: Final fit of factorization
    ########################################################################*/
    proc cpd_als_iterate(const tensors, mats, lambda_vals, const nfactors, const args : cpd_cmd_args) : real
    {
        const nmodes : int = tensors[0].nmodes;
        const nthreads : int = args.numThreads;

        // Setup thread structures. + 64 bytes is to avoid false sharing.
        // SPLATT allocates this in terms of bytes but we're going to assume    
        // the items will always be 8 bytes.
        const lengths: [0..2] int = [((nmodes*nfactors*8)+64)/8, 0,((nmodes*nfactors*8)+64)/8];
        var thds = thd_init(nthreads, 3, lengths);

        // m1 is the output of the MTTKRP
        var m1 = mats[nmodes];

        // aTa is an array of dense_matrics. Each matrix in aTa
        // stores the result of A^T * A where A is the factor matrix
        // for some mode. We have a matrix for each mode, plus an extra
        // for scratch. The size of each matrix is nfactors X nfactors
        var aTa : [0..(nmodes+1)-1] dense_matrix;
        for m in 0..nmodes-1 {
            aTa[m] = new dense_matrix();
            aTa[m].matrix_domain = {0..nfactors-1, 0..nfactors-1};
            aTa[m].I = nfactors;
            aTa[m].J = nfactors;
            aTa[m].vals = 0;
            aTa[m].vals_ref = c_ptrTo(aTa[m].vals);
            // compute A^T*A
            mat_aTa(mats[m], aTa[m]);
        }
        aTa[nmodes] = new dense_matrix();
        aTa[nmodes].matrix_domain = {0..nfactors-1, 0..nfactors-1};
        aTa[nmodes].I = nfactors;
        aTa[nmodes].J = nfactors;
        aTa[nmodes].vals = 0;
        aTa[nmodes].vals_ref = c_ptrTo(aTa[nmodes].vals);

        // MTTKRP workspace
        var mttkrp_ws: splatt_mttkrp_ws = splatt_mttkrp_alloc_ws(tensors, nfactors, args);

        // Compute input tensor norm
        const ttnormsq : real = csf_frobsq(tensors);

        // Fit
        var fit : real = 0;
        var oldfit : real = 0;

        // Separate timer for the iterations
        var itertime : Timer;
        
        // Start CPD timer
        timers_g.timers["CPD"].start();

        // Perform iterations
        const niters = args.numIters;
        for it in 0..niters-1 {
            itertime.clear();
            itertime.start();
            for m in 0..nmodes-1 {
                mats[nmodes].I = tensors[0].dims[m];
                m1.I = mats[m].I;

                // M1 = X * (C o B)
                timers_g.timers["MTTKRP"].start();
                mttkrp_csf(tensors, mats, m, thds, mttkrp_ws, args);
                timers_g.timers["MTTKRP"].stop();

                //par_memcpy(mats[m].vals_ref, m1.vals_ref, m1.I * nfactors);
                //mat_solve_normals(m, nmodes, aTa, mats[m], args.regularization);

                /* M2 = (CtC * BtB *...) ^ -1 */
                calc_gram_inv(m, nmodes, aTa);
                /* A = M1 * M2 */
                c_memset(mats[m].vals_ref, 0, mats[m].I * mats[m].J * 8);
                timers_g.timers["MAT MULT"].start();
                gemm(m1.vals, aTa[nmodes].vals, mats[m].vals, 1, 0); // faster than dot
                timers_g.timers["MAT MULT"].stop();

                if it == 0 {
                    mat_normalize(mats[m], lambda_vals, MAT_NORM_2, thds);
                }
                else {
                    mat_normalize(mats[m], lambda_vals, MAT_NORM_MAX, thds);
                }
                
                // Update A^T*A
                mat_aTa(mats[m], aTa[m]);
            }
            fit = p_calc_fit(nmodes, thds, ttnormsq, lambda_vals, mats, m1, aTa);
            itertime.stop();

            writef("  its =  %i (%0.3drs)  fit = %0.5dr  delta = %+0.4er\n", it+1, itertime.elapsed(), fit, fit-oldfit);

            if fit == 1.0 || (it > 0 && abs(fit-oldfit) < args.tolerance) {
                break;
            }
            oldfit = fit;
        }
        timers_g.timers["CPD"].stop();

        return fit;
    }

    /*########################################################################
    #   Descriptipn:    Compute the fit of a Kruskal tensor, Z, to an input
    #                   tensor, X. This is computed via
    #                   1 - [sqrt(<X,X> + <Z,Z> - 2<X,Z>) / sqrt(<X,X>)]
    #
    #   Parameters:     nmodes (int):               Number of modes 
    #                   thds:                       Thread data structures
    #                   ttnormsq (real):            Norm squared of original input tensor
    #                   lambda_vals (real[]):       Vector of column norms
    #                   mats (dense_matrix[]):      Kruskal-tensor matrices
    #                   m1:                         Result of MTTKRP on last mode
    #                   aTa:                        Matrices for AtA, BtB, CtC, etc.
    #
    #   Return:         real: Inner product of two tensors, computed via
    #                   \lambda^T hadamard(mats[nmodes-1], m1) \lambda
    ########################################################################*/
    private proc p_calc_fit(const nmodes, const thds, const ttnormsq, const lambda_vals, mats, const m1, aTa) : real
    {
        timers_g.timers["FIT"].start();

        /* First get norm of new model: lambda^T * (hada aTa) * lambda */
        const norm_mats : real = p_kruskal_norm(nmodes, lambda_vals, aTa);

        /* Compute inner product of tensor with new model */
        const inner : real = p_tt_kruskal_inner(nmodes, thds, lambda_vals, mats, m1);

        //writef("norm_mats = %dr, inner = %dr\n", norm_mats, inner);

        /*
            We actually want sqrt(<X,X> + <Y,Y> - 2<X,Y>), but if the fit is perfect
            just make it 0
        */
        var residual : real = ttnormsq + norm_mats - (2 * inner);
        if residual > 0.0 {
            residual = sqrt(residual);
        }
        timers_g.timers["FIT"].stop();

        return 1 - (residual / sqrt(ttnormsq));
    }

    /*########################################################################
    #   Descriptipn:    Find the Frobenius norm squared of a Kruskal tensor.
    #                   This is equivalent to via computing <X,X>, the inner
    #                   product of X with itself. We find this via
    #                   \lambda^T (AtA * BtB *...)\lambda where * is the
    #                   Hadamard product.
    #
    #   Parameters:     nmodes (int):               Number of modes 
    #                   lambda_vals (real[]):       Vector of column norms
    #                   aTa:                        Matrices for AtA, BtB, CtC, etc.
    #
    #   Return:         real: The Frobenius norm of X, squared
    ########################################################################*/
    private proc p_kruskal_norm(const nmodes, const lambda_vals, aTa) : real
    {
        const rank = aTa[0].J;
        ref av = aTa[nmodes].vals;

        var norm_mats : real = 0;

        // The loops below assume that the data is stored in the upper
        // triangle

        /* use aTa[nmodes] as scratch space */
        for i in 0..rank-1 {
            for j in i..rank-1 {
                av[i,j] = 1.0;
            }
        }

        /* aTa[nmodes] = hada(aTa) */
        for m in 0..nmodes-1 {
            const ref atavals = aTa[m].vals;
            for i in 0..rank-1 {
                for j in i..rank-1 {
                    av[i,j] *= atavals[i,j];
                }
            }
        }

        /* now compute lambda^T * aTa[nmodes] * lambda */
        for i in 0..rank-1 {
            norm_mats += av[i,i] * lambda_vals[i] * lambda_vals[i];
            for j in i+1..rank-1 {
                norm_mats += av[i,j] * lambda_vals[i] * lambda_vals[j] * 2;
            }
        }

        return abs(norm_mats);
    }

    /*########################################################################
    #   Descriptipn:    Compute the inner product of a Kruskal tensor and an
    #                   unfactored tensor. Assumes that m1 contains the MTTKRP
    #                   result of the last mode.
    #
    #   Parameters:     nmodes (int):               Number of modes 
    #                   thds:                       Thread data structures
    #                   lambda_vals (real[]):       Vector of column norms
    #                   mats (dense_matrix[]):      Kruskal-tensor matrices
    #                   m1:                         Result of MTTKRP on last mode
    #
    #   Return:         real: The inner product of the two tensors, compute via
    #                   1^T hadamard(mats[nmodes-1], m1) \lambda
    ########################################################################*/
    private proc p_tt_kruskal_inner(const nmodes, const thds, const lambda_vals, mats, const m1) : real
    {
        const rank = mats[0].J;
        const lastm = nmodes-1;
        const dim = m1.I;

        const ref m0 = mats[lastm].vals;
        const ref mv = m1.vals;

        var myinner : real = 0;

        const b = new Barrier(numThreads_g);
        coforall tid in 0..numThreads_g-1 with (ref myinner) {
            ref accumF = thds[tid].scratch[0].buf;
            accumF = 0.0;
            // Divide up dim
            const I_per_thread = (dim + numThreads_g - 1) / numThreads_g;
            const I_begin = min(I_per_thread * tid, dim);
            const I_end = min(I_begin + I_per_thread, dim);
            for i in I_begin..I_end-1 {
                for r in 0..rank-1 {
                    accumF[r] += m0[i,r] * mv[i,r];
                }
            }
            b.barrier();

            /* accumulate everything into myinner */
            //TODO: This is pretty nice. In C, the omp parallel
            // section is set up with a +reduction on myinner, but 
            // there needs to be a for-loop that each thread executed to
            // accumulate into my inner. Here, we can use the built in
            // reduction.
            myinner += + reduce(accumF[0..rank-1] * lambda_vals); 
        }
        return myinner;
    }
}

