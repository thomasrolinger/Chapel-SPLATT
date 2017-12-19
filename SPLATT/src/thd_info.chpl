/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/19/2017
    File:   thd_info.chpl

    Description:    This is a module file for thread info
*/

module ThreadInfo {
    /*****************************
    *
    *   scratch_buf class
    *
    ******************************/    
    class scratch_buf {
        var length : int;
        var buf : [0..length-1] real;
    }

    /*****************************
    *
    *   thd_info class
    *
    ******************************/
    class thd_info {
        // I believe the only types we use in scratch
        // are reals
        var nscratch : int;
        var scratch : [0..nscratch-1] scratch_buf;
    }

    /*****************************
    *
    *   Public Functions
    *
    ******************************/

    /*########################################################################
    #   Descriptipn:    Allocate and initialize a number of thd_info classes
    #
    #   Parameters:     nthreads:   Number of threads to allocate for
    #                   nscratch:   Number of scratch arrays to use
    #                   lengths:    Number of 'items' to allocate for each
    #                               scratch array. This is an array itself
    #
    #   Return:         Array of thd_info objects
    ########################################################################*/
    proc thd_init(nthreads, nscratch, lengths)
    {
        // For each thread, we allocate nscratch arrays that have the sizes
        // in lengths
        var thds : [0..nthreads-1] thd_info;
        for t in 0..nthreads-1 {
            thds[t] = new thd_info(nscratch=nscratch);
        }

        for s in 0..nscratch-1 {
            var len = lengths[s];
            for t in 0..nthreads-1 {
                thds[t].scratch[s] = new scratch_buf(length=len);
                thds[t].scratch[s].buf = 0;
            }
        }
        return thds;
    }
}
