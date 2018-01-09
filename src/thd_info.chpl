/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/19/2017
    File:   thd_info.chpl

    Description:    This is a module file for thread info
*/

module ThreadInfo {
    use Base;
    use Barriers;

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
    *   Private Functions
    *
    ******************************/
    
    /*########################################################################
    #   Description:    Performs a parallel sum reduction
    #
    #   Parameters:     Stuff
    #
    #   Return:         None
    ########################################################################*/
    private proc p_reduce_sum(const thds, const scratchid, const nelems, const tid, const b)
    {
        ref myvals = thds[tid].scratch[scratchid].buf;
        
        var half = numThreads_g/2;
        while half > 0 {
            if tid < half && tid + half < numThreads_g {
                const ref target = thds[tid+half].scratch[scratchid].buf;
                for i in 0..nelems-1 {
                    myvals[i] += target[i];
                }
            }
            b.barrier();
            if tid == 0 {
                /* check for odd number */
                if half > 1 && half % 2 == 1 {
                    const ref last = thds[half-1].scratch[scratchid].buf;
                    for i in 0..nelems-1 {
                        myvals[i] += last[i];
                    }
                }
            }
            /* next iteration */
            half /= 2;
        }
        /* account for odd thread at end */
        if tid == 0 {
            if numThreads_g % 2 == 1 {
                const ref last = thds[numThreads_g-1].scratch[scratchid].buf;
                for i in 0..nelems-1 {
                    myvals[i] += last[i];
                }
            }
        }
        b.barrier();
    }

    /*########################################################################
    #   Description:    Performs a parallel max reduction
    #
    #   Parameters:     Stuff
    #
    #   Return:         None
    ########################################################################*/
    private proc p_reduce_max(const thds, const scratchid, const nelems, const tid, const b)
    {
        ref myvals = thds[tid].scratch[scratchid].buf;

        var half = numThreads_g/2;
        while half > 0 {
            if tid < half && tid + half < numThreads_g {
                const ref target = thds[tid+half].scratch[scratchid].buf;
                for i in 0..nelems-1 {
                    myvals[i] = max(myvals[i], target[i]);
                }
            }
            b.barrier();
            if tid == 0 {
                /* check for odd number */
                if half > 1 && half % 2 == 1 {
                    const ref last = thds[half-1].scratch[scratchid].buf;
                    for i in 0..nelems-1 {
                        myvals[i] = max(myvals[i], last[i]);
                    }
                }
            }
            /* next iteration */
            half /= 2;
        }
        /* account for odd thread at end */
        if tid == 0 {
            if numThreads_g % 2 == 1 {
                const ref last = thds[numThreads_g-1].scratch[scratchid].buf;
                for i in 0..nelems-1 {
                    myvals[i] = max(myvals[i], last[i]);
                }
            }
        }
        b.barrier();
    }

    /*****************************
    *
    *   Public Functions
    *
    ******************************/

    /*########################################################################
    #   Description:    Allocate and initialize a number of thd_info classes
    #
    #   Parameters:     nthreads:   Number of threads to allocate for
    #                   nscratch:   Number of scratch arrays to use
    #                   lengths:    Number of 'items' to allocate for each
    #                               scratch array. This is an array itself
    #
    #   Return:         Array of thd_info objects
    ########################################################################*/
    proc thd_init(const nthreads, const nscratch, const lengths)
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

    /*########################################################################
    #   Description:    Performs a reduction 
    #
    #   Parameters:     Stuff
    #
    #   Return:         None
    ########################################################################*/
    proc thd_reduce(const thds, const scratchid, const nelems, const tid, const b, const which)
    {
        if numThreads_g == 1 {
            return;
        }
        b.barrier();
        select which {
            when REDUCE_SUM {
                p_reduce_sum(thds, scratchid, nelems, tid, b);
            }
            when REDUCE_MAX {
                p_reduce_max(thds, scratchid, nelems, tid, b);
            }
        }
    }

}
