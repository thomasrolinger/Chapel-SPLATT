/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/14/2017
    File:   mats.chpl

    Description:    This is a module file for things related
                    to dense matrices that we use.
*/

module Matrices {
    use Base;
    use BLAS;
    use LAPACK;
    use ClassicLAPACK;
    use Barriers;
    use IO.FormattedIO;

    /*****************************
    *
    *   dense_matrix class
    *
    ******************************/
    class dense_matrix {
        var matrix_domain : domain(2) = {0..1, 0..1};
        // Holds the actual data
        var vals : [matrix_domain] real;
        // The number of rows (I) and cols (J). Not sure
        // if we'll need these since we can get that info from
        // the domain but it's nice to have.
        var I : int;
        var J : int;
        // This is a ptr to vals. For most cases, we don't need this
        // but when doing the MTTKRP, we need to get a reference to
        // one of the matrix's vals array and assign it to a private
        // copy.
        var vals_ref : c_ptr(real);
    }

    /*****************************
    *
    *   Private Functions
    *
    ******************************/
    /*########################################################################
    #   Description:    Form the Gram matrix from A^T * A
    #
    #   Parameters:     neq_matrix (real[]):    The matrix to fill
    #                   aTa (dense_matrix[]):   Individual Gram matrices
    #                   mode (int):             Which mode we are doing
    #                   nmodes (int):           Number of modes in tensor
    #                   reg (real):             Regularization param
    #
    #   Return:         None
    ########################################################################*/
    private proc p_form_gram(neq_matrix, aTa, mode, nmodes, reg)
    {
        /* nfactors */
        var N = aTa[0].J;

        /* 
            form upper triangual normal equations. We are ignoring
            the regularization parameter, so this loop is simplified
        */
        ref neqs = neq_matrix.vals;
            
        /* first initialize with 1s */
        forall (i,j) in neqs.domain {
            neqs(i,j) = 1.0;
        }
        /* now Hadamard product of all aTa matrices */
        for m in 0..nmodes-1 {
            if m == mode {
                continue;
            }
            var mat = aTa[m].vals;
            forall i in 0..N-1 {
                /* mat is symmetric but stored upper right triangular */
                /* copy upper triangle */
                for j in i..N-1 {
                    neqs[i,j] *= mat[i,j];
                }
            }
        } /* for each mode */

        /* now copy lower triangle */
        forall i in 0..N-1 {
            for j in 0..i-1 {
                neqs[i,j] = neqs[j,i];
            }
        }
    }

    /*########################################################################
    #   Description:    Calculates 2-norm
    #
    #   Parameters:     Stuff
    #
    #   Return:         None
    ########################################################################*/
    private proc p_mat_2norm(A, lambda_vals, thds)
    {
        var I = A.I;
        var J = A.J;
        ref vals = A.vals;

        var b = new Barrier(numThreads_g);

        coforall tid in 0..numThreads_g-1 {
            ref mylambda = thds[tid].scratch[0].buf;
            for j in 0..J-1 {
                mylambda[j] = 0;
            }   

            /*
                mylambda[j] will contain the sum of the square of all the values
                in the j-th column of vals
            */
            /*forall j in 0..J-1 {
                lambda_vals[j] = (+ reduce vals[0..I-1, j]**2);
            }*/

            /*
                With OMP, the coforall is a parallel section and this for-loop is just a
                omp for. The way it works is that each thread executes the statements within
                the parallel section and then the omp for iterations are divided amongst those
                threads. However in Chapel, it appears that the forall is NOT using the same pool
                of threads that were "created" by the coforall. Therefore, we must manually divide
                up the for loop iterations based on the TIDs.
            */
            // Divide up the outer loop iterations (rows)
            var I_per_thread = (I + numThreads_g - 1) / numThreads_g;
            var I_begin = min(I_per_thread * tid, I);
            var I_end = min(I_begin + I_per_thread, I);
            for i in I_begin..I_end-1 {
                for j in 0..J-1 {
                    mylambda[j] += vals(i,j) * vals(i,j);
                }
            }
            
            // reduction on partial sums
            thd_reduce(thds, 0, J, tid, b, REDUCE_SUM);
            
            if tid == 0 {
                lambda_vals = mylambda[0..J-1];
            }
        
            b.barrier();

            // Divide up J amongst threads
            var J_per_thread = (J + numThreads_g - 1) / numThreads_g;
            var J_begin = min(J_per_thread * tid, J);
            var J_end = min(J_begin + J_per_thread, J);
            for j in J_begin..J_end-1 {
                lambda_vals[j] = sqrt(lambda_vals[j]);
            }   
            //TODO: Does sqrt() work in parallel?

            /* do the normalization */
            for i in I_begin..I_end-1 {
                for j in 0..J-1 {
                    vals(i,j) /= lambda_vals[j];
                }
            }
        }
    }

    /*****************************
    *
    *   Public Functions
    *
    ******************************/
    /*########################################################################
    #   Description:    Computes A^T*A via BLAS syrk
    #
    #   Parameters:     A (dense_matrix):   The matrix to operate on
    #                   ret (dense_matrix): The matrix to store the result in
    #
    #   Return:         None
    ########################################################################*/
    proc mat_aTa(A : dense_matrix, ret : dense_matrix)
    {
        /*  
            In SPLATT, uplo is L and trans is N. I believe this is the case
            because SPLATT stores the dense matrices in row-major order.
            In Chapel, we are storing the matrices as multi-dimensional
            arrays.
        */  
        timers_g.timers["MAT A^TA"].start();
        var uplo = Uplo.Upper;
        var trans = Op.T;
        var order = Order.Row;
        var alpha : c_double = 1.0;
        var beta : c_double = 0.0;
        syrk(A.vals, ret.vals, alpha, beta, uplo, trans, order);
        timers_g.timers["MAT A^TA"].stop();
    }

    /*########################################################################
    #   Description:    Calculates (BtB * CtC *...)^-1 where * is the Hadamard
    #                   product. This is the Gram matrix of the CPD.
    #
    #   Parameters:     mode (int): Which mode we are operating on
    #                   nmodes (int):   Number of modes in tensor
    #                   aTa (dense_matrix[]):   Array of matrices that contains
    #                                           BtB, CtC, etc.
    #                   rhs (dense_matrix);     Factor matrix for this mode
    #                   reg (real):             Regularization value
    #
    #   Return:         None
    ########################################################################*/
    proc mat_solve_normals(mode, nmodes, aTa, rhs, reg)
    {
        timers_g.timers["INVERSE"].start();

        p_form_gram(aTa[nmodes], aTa, mode, nmodes, reg);
        
        var uplo = "U"; 
        ref neqs = aTa[nmodes].vals;

        /* Cholesky factorization */
        potrf(lapack_memory_order.row_major, uplo, neqs);
        
        // Solve against RHS 

        /*
          Not quite sure what's going on here. The standard defintion for potrs says
          that nrhs is equal to the number of columns in BX (i.e. rhs.vals). However,
          SPLATT sets nrhs to be the number of rows in BX. By using the built-in Chapel
          wrapper for potrs, we only solve the first N equations, where N is the number of
          columns in rhs.vals (i.e. number of factors). Furthermore, SPLATT sets ldb to be
          N. However, the only way I can get the same results as SPLATT is by calling the LAPACKE
          dpotrs function directly and setting nrhs to the number of rows in BX and setting ldb
          to the number of rows in BX. I am not sure why this is so different from what SPLATT
          is doing, since it is also directly calling LAPACK.
        */
        var N: c_int = neqs.domain.dim(1).size : c_int;
        var nrhs : c_int = rhs.vals.domain.dim(1).size : c_int;
        var lda : c_int = neqs.domain.dim(2).size : c_int;
        var ldb : c_int = rhs.vals.domain.dim(1).size : c_int;
        LAPACKE_dpotrs(lapack_memory_order.row_major, uplo, N, nrhs, neqs, lda, rhs.vals, ldb);
        timers_g.timers["INVERSE"].stop();
    }
    
    /*########################################################################
    #   Description:    Normalize the columns of A and return the norms in
    #                   lambda_vals. Supported norms are 2-norm and max-norm
    #
    #   Parameters:     A (dense_matrix):       The matrix to normalize
    #                   lambda_vals (reals[]:   Vector of columns norms
    #                   which (int):            Which norm to use
    #
    #   Return:         None
    ########################################################################*/
    proc mat_normalize(A, lambda_vals, which, thds)
    {
        timers_g.timers["MAT_NORM"].start();
        select (which) {
            when MAT_NORM_2 {
                //TODO: Implement
                p_mat_2norm(A, lambda_vals, thds);
            }
            when MAT_NORM_MAX {
                //TODO: Implement
                //p_mat_maxnorm(A, lambda_vals, thds);
            }
        }
        timers_g.timers["MAT_NORM"].stop();
    }

}
