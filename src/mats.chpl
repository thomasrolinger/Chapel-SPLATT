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
    use Norm;
    //use LAPACK;
    use LinearAlgebra;
    //use ClassicLAPACK;
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
    private proc p_form_gram(const neq_matrix, const aTa, const mode, const nmodes, const reg)
    {
        /* nfactors */
        const N = aTa[0].J;

        /* 
            form upper triangual normal equations. We are ignoring
            the regularization parameter, so this loop is simplified
        */
        const ref neqs = neq_matrix.vals;
          
        const b = new Barrier(numThreads_g);

        coforall tid in 0..numThreads_g-1 {
            /* first initialize with 1s */
            const I_per_thread = (N + numThreads_g - 1) / numThreads_g;
            const I_begin = min(I_per_thread * tid, N);
            const I_end = min(I_begin + I_per_thread, N);
            for i in I_begin..I_end-1 {
                for j in 0..N-1 {
                    neqs(i,j) = 1.0;
                }
            }
            /* now Hadamard product of all aTa matrices */
            for m in 0..nmodes-1 {
                if m == mode {
                    continue;
                }
                const ref mat = aTa[m].vals;
                for i in I_begin..I_end-1 {
                    /* mat is symmetric but stored upper right triangular */
                    /* copy upper triangle */
                    for j in i..N-1 {
                        neqs[i,j] *= mat[i,j];
                    }
                }
            } /* for each mode */

            b.barrier();

            /* now copy lower triangle */
            //for i in 0..N-1 {
            for i in I_begin..I_end-1 {
                for j in 0..i-1 {
                    neqs[i,j] = neqs[j,i];
                }
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
    private proc p_mat_2norm(const A, lambda_vals, const thds)
    {
        const I = A.I;
        const J = A.J;
        ref vals = A.vals;

        const b = new Barrier(numThreads_g);

        coforall tid in 0..numThreads_g-1 {
            ref mylambda = thds[tid].scratch[0].buf;
            mylambda = 0.0;

            /*
                mylambda[j] will contain the sum of the square of all the values
                in the j-th column of vals
                
                With OMP, the coforall is a parallel section and this for-loop is just a
                omp for. The way it works is that each thread executes the statements within
                the parallel section and then the omp for iterations are divided amongst those
                threads. However in Chapel, it appears that the forall is NOT using the same pool
                of threads that were "created" by the coforall. Therefore, we must manually divide
                up the for loop iterations based on the TIDs.
            */
            // Divide up the outer loop iterations (rows)
            const I_per_thread = (I + numThreads_g - 1) / numThreads_g;
            const I_begin = min(I_per_thread * tid, I);
            const I_end = min(I_begin + I_per_thread, I);
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
            /*const J_per_thread = (J + numThreads_g - 1) / numThreads_g;
            const J_begin = min(J_per_thread * tid, J);
            const J_end = min(J_begin + J_per_thread, J);
            for j in J_begin..J_end-1 {
                lambda_vals[j] = sqrt(lambda_vals[j]);
            }*/
            lambda_vals = sqrt(lambda_vals);   
            

            //TODO: in SPLATT< the above loop is a omp for within an omp parallel section.
            // There is an implicit barrier at the end of that loop that prevents all the
            // threads in the paralle section from continuing. Here, need to do the barrier
            // explicitly.
            b.barrier();
            /* do the normalization */
            for i in I_begin..I_end-1 {
                for j in 0..J-1 {
                    vals(i,j) /= lambda_vals[j];
                }
            }
        }
    }

    /*########################################################################
    #   Description:    Calculates max-norm
    #
    #   Parameters:     Stuff
    #
    #   Return:         None
    ########################################################################*/
    private proc p_mat_maxnorm(const A, lambda_vals, const thds)
    {
        const I = A.I;
        const J = A.J;
        ref vals = A.vals;

        const b = new Barrier(numThreads_g);

        coforall tid in 0..numThreads_g-1 {
            ref mylambda = thds[tid].scratch[0].buf;
            mylambda = 0.0;

            // Divide up the outer loop iterations (rows)
            const I_per_thread = (I + numThreads_g - 1) / numThreads_g;
            const I_begin = min(I_per_thread * tid, I);
            const I_end = min(I_begin + I_per_thread, I);
            for i in I_begin..I_end-1 {
                for j in 0..J-1 {
                    mylambda[j] = max(mylambda[j], vals(i,j));
                }
            }

            // reduction on partial maxes
            thd_reduce(thds, 0, J, tid, b, REDUCE_MAX);
    
            if tid == 0 {
                lambda_vals = mylambda[0..J-1];
            }

            b.barrier();

            // Divide up J amongst threads
            const J_per_thread = (J + numThreads_g - 1) / numThreads_g;
            const J_begin = min(J_per_thread * tid, J);
            const J_end = min(J_begin + J_per_thread, J);
            for j in J_begin..J_end-1 {
                lambda_vals[j] = max(lambda_vals[j], 1.0);
            }

            b.barrier();

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
    proc mat_aTa(const A : dense_matrix, const ret : dense_matrix)
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
    
        forall i in 1..ret.I-1 {
            for j in 0..i-1 {
                ret.vals(i,j) = ret.vals(j,i);
            }
        }

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
    /*proc mat_solve_normals(mode, nmodes, aTa, rhs, reg)
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
    }*/
    
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
    proc mat_normalize(A, lambda_vals, const which, const thds)
    {
        timers_g.timers["MAT NORM"].start();
        select (which) {
            when MAT_NORM_2 {
                p_mat_2norm(A, lambda_vals, thds);
            }
            when MAT_NORM_MAX {
                p_mat_maxnorm(A, lambda_vals, thds);
            }
        }
        timers_g.timers["MAT NORM"].stop();
    }

    /*############################################################################

    ##############################################################################*/
    proc calc_gram_inv(const mode, const nmodes, const aTa)
    {
        timers_g.timers["INVERSE"].start();

        const rank = aTa[0].J;
        ref av = aTa[nmodes].vals;
        
        //av = 1.0;  DOING THIS SEEMS TO BE SLIGHTLY SLOWER
        for x in 0..rank-1 {
            for y in 0..rank-1 {
                av[x,y] = 1.0;
            }
        }

        /* hadamard */
        for m in 1..nmodes-1 {
            const madjust = (mode + m) % nmodes;
            const ref vals = aTa[madjust].vals;
            //av *= vals; DOING THIS SEEMS TO BE SLIGHTLY SLOWER
            for x in 0..rank-1 {
                for y in 0..rank-1 {
                    av[x,y] *= vals[x,y];
                }
            }
        }

        /* M2 = M2^-1 */
        mat_syminv(aTa[nmodes]);

        timers_g.timers["INVERSE"].stop();
    }

    /*############################################################################

    ##############################################################################*/
    proc mat_syminv(const A) 
    {
        assert(A.I == A.J);

        const N = A.I;
        var L : dense_matrix = new dense_matrix();
        L.I = N;
        L.J = N;
        L.matrix_domain = {0..L.I-1, 0..L.J-1};
        L.vals = 0;

        /* Cholesky */
        mat_cholesky(A, L);
        
        /* setup identity matrix */
        A.vals = 0.0;
        for n in 0..N-1 {
            A.vals(n,n) = 1.0;
        }

        /* Solve L*Y = I */
        p_mat_forwardsolve(L, A);

        /* transpose L */
        L.vals = transpose(L.vals);

        /* Solve U*A = Y */
        p_mat_backwardsolve(L, A);
    }

    /*############################################################################

    ##############################################################################*/
    proc mat_cholesky(const A, const L)
    {
        /* check dimensions */
        assert(A.I == A.J);
        assert(A.I == L.J);
        assert(L.I == L.J);

        const N = A.I;
        const ref av = A.vals;
        ref lv = L.vals;

        for i in 0..N-1 {
            for j in 0..i {
                var inner = 0.0;
                for k in 0..j-1 {
                    inner += lv(i,k) * lv(j,k);
                }
                if i == j {
                    lv(i,j) = sqrt(av(i,i) - inner);
                }
                else {
                    lv(i,j) = 1.0 / lv(j,j) * (av(i,j) - inner);
                }
            }
        }
    }   

    private proc p_mat_forwardsolve(const L, const B)
    {
        const N = L.I;
        const ref lv = L.vals;
        ref bv = B.vals;

        /* first row of X is easy */
        //bv /= lv[0,0];
        for j in 0..N-1 {
            bv[0,j] /= lv[0,0];
        }
        
        /* now do forward sub */
        for i in 1..N-1 {
            for j in 0..i-1 {
                for f in 0..N-1 {
                    bv(i,f) -= lv(i,j) * bv(j,f);
                }
            }
            for f in 0..N-1 {
                bv(i,f) /= lv(i,i);
            }
        }
    } 

    private proc p_mat_backwardsolve(const U, const B)
    {
        const N = U.I;
        const ref rv = U.vals;
        ref bv = B.vals;

        /* last row of X is easy */
        for f in 0..N-1 {
            var i = N-1;
            bv(i,f) /= rv(i,i);
        }   
 
        /* now do backward sub */
        for row in 2..N {
            var i = N - row;
            for j in (i+1)..N-1 {
                for f in 0..N-1 {
                    bv(i,f) -= rv(i,j) * bv(j,f);
                }
            }
            for f in 0..N-1 {
                bv(i,f) /= rv(i,i);
            }    
        }
    }  

    proc mat_matmul(const A, const B, const C)
    {
        C.I = A.I;
        C.J = B.J;
        const ref av = A.vals;
        const ref bv = B.vals;
        const M = A.I;
        const N = B.J;
        const Na = A.J;
        ref cv = C.vals;
        param TILE = 16;
        forall i in 0..M-1 {
            for jt in 0..N-1 by TILE {
                for kt in 0..Na-1 by TILE {
                    const JSTOP = min(jt+TILE, N);
                    for j in jt..JSTOP-1 {
                        var accum : real = 0.0;
                        const KSTOP = min(kt+TILE, Na);
                        for k in kt..KSTOP-1 {
                            accum += av[i,k] * bv[k,j];
                        }
                        cv[i,j] += accum;
                    }
                }
            }
        }
    }
}

