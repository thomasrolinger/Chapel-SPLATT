/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/19/2017
    File:   thread_partition.chpl

    Description:    This is a module file for things related to
                    thread partitioning in SPLATT.
*/

module ThreadPartition {
    use Base;
    use Assert;

    /*****************************
    *
    *   Public Functions
    *
    ******************************/
    /*########################################################################
    #   Descriptipn:    Partition weighted items among threads in a load
    #                   balanced manner. This uses chains-on-chains 
    #                   partitioning to optimally load balance
    #
    #   Parameters:     weights (int[]):    Weights of each item
    #                   nitems (int):       Number of items to partition
    #                   nparts (int):       Number of partitions
    #                   bottleneck (int[]): Weight of the heaviest partition
    #
    #   Return:         A balanced partitioning of the data. Thread 't' should
    #                   process items [t, t+1).
    ########################################################################*/
    proc partition_weighted(weights, nitems, nparts, ref bottleneck)
    {
        timers_g.timers["PART1D"].start();
        prefix_sum_inc(weights, nitems);
        var parts : [0..(nparts+1)-1] int;
        nprobes_g = 0;
        var bneck : int = 0;
        
        // actual partitioning
        if nitems > nparts {
            // use recursive bisectioning with 0 tolerance to get exact solution
            bneck = p_eps_rb_partition_1d(weights, nitems, parts, nparts, 0);
            // apply partitioning that we found
            var success : bool = lprobe(weights, nitems, parts, nparts, bneck);
            assert(success == true);
        }
        else {
            // Do a trivial partitioning. Silly, but this happens when tensors
            // have short modes
            for p in 0..nitems-1 {
                parts[p] = p;
                bneck = max(bneck, weights[p]);
            }
            for p in nitems..nparts {
                parts[p] = nitems;
            }
        }
        bottleneck = bneck;
        timers_g.timers["PART1D"].stop();
        return parts;
    }

    /*########################################################################
    #   Descriptipn:    Something
    #
    #   Parameters:     blah
    #
    #   Return:         None
    ########################################################################*/ 
    proc prefix_sum_inc(weights, nitems)
    {
        for x in 1..nitems-1 {
            weights[x] += weights[x-1];
        }
    }

    /*########################################################################
    #   Descriptipn:    Something
    #
    #   Parameters:     blah
    #
    #   Return:         bool
    ########################################################################*/
    proc lprobe(weights, nitems, parts, nparts, bottleneck)
    {
        nprobes_g += 1;
        var wtotal : int = weights[nitems-1];
        // initialize partitioning
        parts[0] = 0;
        for p in 1..nparts {
            parts[p] = nitems;
        }

        var bsum : int = bottleneck;
        var step : int = nitems / nparts;
        for p in 1..nparts-1 {
            // jump to the next bucket      
            while step < nitems && weights[step] < bsum {
                step += nitems / nparts;
            }
            // find the end (exclusive) index of process p
            parts[p] = p_binary_search(weights, step-(nitems/nparts), min(step,nitems), bsum);
            if parts[p] == nitems {
                // check for pathological case when last weight is larger than bottleneck
                var size_last = weights[nitems-1] - weights[parts[p-1]-1];
                return size_last < bottleneck;
            }
            bsum = weights[parts[p]-1] + bottleneck;
        }
        return bsum >= wtotal;
    }

    /*****************************
    *
    *   Prive Functions
    *
    ******************************/

    /*########################################################################
    #   Descriptipn:    Something
    #
    #   Parameters:     blah
    #
    #   Return:         int
    ########################################################################*/
    private proc p_eps_rb_partition_1d(weights, nitems, parts, nparts, eps)
    {
        var tot_weight : int = weights[nitems-1];
        var lower : int = tot_weight/ nparts;
        var upper : int = tot_weight;
        do {
            var mid : int = lower + ((upper-lower)/2);
            if lprobe(weights, nitems, parts, nparts, mid) {
                upper = mid;
            }
            else {
                lower = mid+1;
            }
        } while upper > lower+eps;
    
        return upper;
    }

    /*########################################################################
    #   Description:    Perform a binary search on an array for a value.
    #
    #   Parameters:     weights (int[]):    Array to search
    #                   lo (int):           Lower bound 
    #                   hi (int):           Upper bound
    #                   target (int):       Target value
    #
    #   Return:         The index j, where weights[j] <= target && weights[j+1] > target
    ########################################################################*/
    private proc p_binary_search(weights, lo, hi, target) : int
    {
        var left = lo;
        var right = hi;
        while (right-left) > 8 {
            var mid : int = left + ((right - left) / 2);
            if weights[mid] <= target && weights[mid+1] > target {
                return mid;
            }
            if weights[mid] < target {
                left = mid+1;
            }
            else {
                right = mid;
            }
        }
        return p_linear_search(weights, left, right, target);
    }

    /*########################################################################
    #   Description:    Perform a linear search on an array for a value.
    #
    #   Parameters:     weights (int[]):    Array to search
    #                   left (int):         Lower bound 
    #                   right (int):        Upper bound
    #                   target (int):       Target value
    #
    #   Return:         The index j, where weights[j] <= target && weights[j+1] > target
    ########################################################################*/
    private proc p_linear_search(weights, left, right, target) : int
    {
        for x in left..(right-1)-1 {
            if target < weights[x+1] {
                return x+1;
            }
        }
        return right;
    }

}
