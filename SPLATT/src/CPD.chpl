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
    proc splatt_cpd_als(tensors, args : cpd_cmd_args, factored : splatt_kruskal) : int
    {
        var nmodes : int = tensors[0].nmodes;
        var nfactors : int = args.decompRank;
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
    proc cpd_als_iterate(tensors, mats, lambda_vals, nfactors, args : cpd_cmd_args) : real
    {
        var nmodes : int = tensors[0].nmodes;
        var nthreads : int = args.numThreads;

        // Setup thread structures. + 64 bytes is to avoid false sharing.
        // SPLATT allocates this in terms of bytes but we're going to assume    
        // the items will always be 8 bytes.
        var lengths: [0..2] int = [((nmodes*nfactors*8)+64)/8, 0,((nmodes*nfactors*8)+64)/8];
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
        var ttnormsq : real = csf_frobsq(tensors);

        // Separate timer for the iterations
        var itertime : Timer;
        
        // Start CPD timer
        timers_g.timers["CPD"].start();

        // Perform iterations
        var niters = args.numIters;
        for it in 0..niters-1 {
            itertime.clear();
            itertime.start();
            for m in 0..nmodes-1 {
                mats[nmodes].I = tensors[0].dims[m];
                m1.I = mats[m].I;

                // M1 = X * (C o B)
                writeln("## MODE: ", m, ", Calling mttkrp_csf");
                timers_g.timers["MTTKRP"].start();
                mttkrp_csf(tensors, mats, m, thds, mttkrp_ws, args);
                timers_g.timers["MTTKRP"].stop();
                writeln("M1:");
                for x in 0..9 {
                    write(mats[nmodes].vals_ref[x], " ");
                }
                writeln("");
                exit(-1);
            }
            itertime.stop();
        }

        timers_g.timers["CPD"].stop();

        return 0.0;
    }
}
