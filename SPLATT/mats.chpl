/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/14/2017
    File:   mats.chpl

    Description:    This is a module file for things related
                    to dense matrices that we use.
*/

module Matrices {
    use BLAS;

    /*****************************
    *
    *   dense_matrix class
    *
    ******************************/
    class dense_matrix {
        var matrix_domain : domain(2) = {0..1, 0..1};
        // Holds the actual data (2D matrix)
        var vals : [matrix_domain] real;
        // The number of rows (I) and cols (J). Not sure
        // if we'll need these since we can get that info from
        // the domain but it's nice to have.
        var I : int;
        var J : int;
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
    proc mat_aTa(A, ret)
    {
        var alpha = 1.0;
        var beta = 0.0;
        var uplo = Uplo.Lower;
        var trans = Op.N;
        var order = Order.Row;
        syrk(A.vals, ret.vals, alpha, beta, uplo, trans, order);
    }
}
