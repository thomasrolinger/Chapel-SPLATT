/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/07/2017
    File:   stats.chpl

    Description:    This is a module file that contains functions to
                    print various stats
*/

module Stats {
    use sptensor;
    use Util;
    /*########################################################################
    #   Descriptipn:    Output basic statistics about the sparse tensor tt.
    #
    #   Parameters:     tt (sptensor_t):    The tensor to get stats for
    #                   fname (string):     File name which tt came from
    #
    #   Return:         None
    ########################################################################*/
    proc stats_tt(tt : sptensor_t, fname : string)
    {
        writeln("Tensor information ---------------------------------------------");
        writeln("FILE=", fname);
        write("DIMS=", tt.dims[0]);
        for m in 1..tt.nmodes-1 {
            write("x", tt.dims[m]);
        }
        write(" NNZ=", tt.nnz, " DENSITY=", tt.tt_calcDensity());
        writeln("");
        writeln("COORD-STORAGE=", bytes_str(tt.nnz * ((8 * tt.nmodes) + 8)));
        writeln("");
        writeln("");
    }
}
