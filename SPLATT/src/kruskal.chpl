/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/14/2017
    File:   kruskal.chpl

    Description:    This is a module file for the splatt_kruskal
                    struct/class and associated functions.
*/

module Kruskal {
    use Base;
    use Matrices;

    /*****************************
    *
    *   splatt_kruskal class
    *
    ******************************/
    /*
        Kruskal tensors are the output of the CPD. Each
        mode of the tensor is represented as a matrix with
        unit columns. Lambda is a vector whose entries scale
        the columns of the matrix factors.
    */
    class splatt_kruskal {
        var rank : int;     // The rank of the decomposition
        // The factor matrices, one for each mode. Each
        // factor matrix is defined above in the dense_matrix class
        var factors : [MAX_NMODES] dense_matrix; 
        var lambda_vals : [LAMBDA_d] real;
        var nmodes : int;
        var dims : [NUM_MODES_d] int; // Number of rows in each factor matrix
        var fit : real ; // quality of the CPD.
    }
}
