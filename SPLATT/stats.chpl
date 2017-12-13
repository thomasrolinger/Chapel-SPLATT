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
    use IO.FormattedIO;
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

    /*########################################################################
    #   Descriptipn:    Output work-related statistics before a CPD factorization.
    #                   This includes rank, #threads, tolerance, etc.
    #
    #   Parameters:     csf (splatt_csf[]): The CSF tensor we are factoring
    #                   opts (cpd_cmd_args):Arguments given to this program
    #
    #   Return:         None
    ########################################################################*/
    proc cpd_stats(csf, opts)
    {
        // Find total storage
        var fbytes = csf_storage(csf, opts);
        var mbytes : int = 0;
        for m in 0..csf[0].nmodes-1 {
            mbytes += csf[0].dims[m] * opts.decompRank * 8;
        }

        // header
        writeln("Factoring ------------------------------------------------------");
        writef("NFACTORS=%i MAXITS=%i TOL=%0.1er REG=%0.1er", opts.decompRank,  opts.numIters,
                                                            opts.tolerance, opts.regularization);
        writef(" SEED=%i THREADS=%i\n", opts.rndSeed, opts.numThreads);
        // CSF allocation
        write("CSF-ALLOC=");
        select(opts.numCSF) {
            when "one" {
                write("ONEMODE");
            }
            when "two" {
                write("TWOMODE");
            }
            when "all" {
                write("ALLMODE");
            }
        }
        // Tiling info
        write(" TILE=");
        if opts.tiling == 0 {
            writeln("NO");
        }
        else {
            writeln("DENSE TILE-DEPTH=1");
        }
        var fstorage = bytes_str(fbytes);
        var mstorage = bytes_str(mbytes);
        writeln("CSF-STORAGE=", fstorage, " FACTOR-STORAGE=", mstorage);
        writeln("");
    }
}
