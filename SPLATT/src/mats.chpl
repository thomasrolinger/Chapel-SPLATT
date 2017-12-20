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
    use BLAS.C_BLAS;

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
    *   Public Functions
    *
    ******************************/
    /*########################################################################
    #   Descriptipn:    Computes A^T*A via BLAS syrk
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
}
