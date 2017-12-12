/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/08/2017
    File:   splatt_sort.chpl

    Description:    This is a module file for sorting related to SPLATT
*/

module splatt_sort {
    use Sort;
    use sptensor;
    use Base;
    use Barriers;

    /**************************************************************************
    *
    *   PUBLIC FUNCTIONS
    *
    ***************************************************************************/

    /*########################################################################
    #   Descriptipn:    Sort a tensor using a permutation of its modes. Sorting
    #                   uses dim_perm to order modes by decreasing priority. If
    #                   dim_perm = {1, 0, 2} then nonzeros will be ordered by
    #                   ind[1] with ties broken by ind[0] and finally deferring
    #                   to ind[2].
    #
    #   Parameters:     tt (sptensor_t):    The sptensor_t to sort
    #                   mode (int):         The primary for sorting
    #                   dim_perm([] int):   A permutation array that
    #                                       defines sorting priority. If nil,
    #                                       a default ordering of {0,1,..,m} 
    #                                       is used.
    #
    #   Return:         None
    ########################################################################*/
    proc tt_sort(tt : sptensor_t, mode : int, dim_perm)
    {
        tt_sort_range(tt, mode, dim_perm, 0, tt.nnz);
    }

    /*########################################################################
    #   Descriptipn:    Sort a tensor using tt_sort on only a range of the
    #                   nonzero elements. Nonzeros in the range [start,end) will
    #                   be sorted.
    #
    #   Parameters:     tt (sptensor_t):    The sptensor_t to sort
    #                   mode (int):         The primary for sorting
    #                   dim_perm([] int):   A permutation array that
    #                                       defines sorting priority. If nil,
    #                                       a default ordering of {0,1,..,m} 
    #                                       is used.
    #                   start (int):        The first nonzero to include in
    #                                       the sorting
    #                   end (int):          The end of the nonzeros to sort (exclusive)
    #
    #   Return:         None
    ########################################################################*/
    proc tt_sort_range(tt : sptensor_t, mode : int, dim_perm : [NUM_MODES_d] int,
                       start : int, end : int)
    {
        var cmplt : [NUM_MODES_d] int;
        // if dim_perm == NULL
        if !c_ptrTo(dim_perm) {
            cmplt[0] = mode;
            for m in 1..tt.nmodes-1 {
                cmplt[m] = (mode+m) % tt.nmodes;
            }
        }
        else {
            cmplt = dim_perm;
        }

        if start == 0 && end == tt.nnz {
            p_counting_sort_hybrid(tt, cmplt);
        }
        else {
            // sort a subtensor
            var cmplt_ptr = c_ptrTo(cmplt);
            select(tt.tensorType) {
                when tt_type.SPLATT_NMODE {
                    p_tt_quicksort(tt, cmplt_ptr, start, end);
                } 
                when tt_type.SPLATT_3MODE {
                    p_tt_quicksort3(tt, cmplt_ptr, start, end);
                }
            }
        }
    }

    /**************************************************************************
    *
    *   PRIVATE FUNCTIONS
    *
    ***************************************************************************/
    
    /*########################################################################
    #   Descriptipn:    Perform a counting sort on the most significant mode
    #                   (cmplt[0]) and then parallel quicksorts on each of the
    #                   slices.
    #
    #   Parameters:     tt (sptensor_t):    The sptensor_t to sort
    #                   cmplt ([] int):     Mode permutation used for defining
    #                                       tie-breaking order
    #
    #   Return:         None
    ########################################################################*/
    private proc p_counting_sort_hybrid(tt : sptensor_t, cmplt : [NUM_MODES_d] int)
    {
        var m : int = cmplt[0];
        var nslices : int = tt.dims[m];
        
        // Matrix with nmodes rows and nnz columns
        var new_ind : [COORD_d] int;

        // Array with nnz elements
        var new_vals : [NNZ_d] real;

        // Histogram array that has nslices*numThreads+1 elements
        var histogram_array : [0..((nslices*numThreads_g)+1)-1] int;
        var arrayPtr = c_ptrTo(histogram_array);

        var b = new Barrier(numThreads_g);

        writeln("Before coforall in counting_sort");
        // For each thread, give it a portion of the mode to sort
        coforall tid in 0..numThreads_g-1 {
            /*
                Get the offset into histogram_array where this TID is.
                Don't know how to do this in native Chapel, so get a
                pointer to histogram_array and do it that way. Note that
                with this way, we can't treat histogram has a Chapel array.
            */
            var histogram = arrayPtr + (nslices * tid);
            // Zero out array
            c_memset(histogram, 0, nslices*8);
            var j_per_thread: int = (tt.nnz + numThreads_g - 1) / numThreads_g;
            var jbegin : int = min(j_per_thread*tid, tt.nnz);
            var jend : int = min(jbegin + j_per_thread, tt.nnz);
            writeln("j_per_thread: ", j_per_thread);
            writeln("jbegin: ", jbegin);
            writeln("jend: ", jend);

            // count
            for j in jbegin..jend-1 {
                var idx: int = tt.ind[m,j];
                writeln("[", j, "] idx = ", idx);
                histogram[idx] += 1;
            }
            
            b.barrier();
            exit(-1);
            // prefix sum
            for j in (tid*nslices)+1..((tid+1)*nslices)-1 {
                var transpose_j : int = p_transpose_idx(j, numThreads_g, nslices);
                var transpose_j_minus_1 : int = p_transpose_idx(j-1, numThreads_g, nslices);
                histogram_array[transpose_j] += histogram_array[transpose_j_minus_1];
            }

            b.barrier();

            // In SPLATT, the following loop is only performed by the master thread.
            // We'll just assume this is tid=0.
            if tid == 0 {
                for t in 1..numThreads_g-1 {
                    var j0 : int = (nslices*t)-1;
                    var j1 : int = nslices * (t+1) -1;
                    var transpose_j0 : int = p_transpose_idx(j0, numThreads_g, nslices);
                    var transpose_j1 : int = p_transpose_idx(j1, numThreads_g, nslices);
                    histogram_array[transpose_j1] += histogram_array[transpose_j0];
                }
            }

            b.barrier();

            // Every thread except for tid=0
            if tid > 0 {
                var transpose_j0 : int = p_transpose_idx(nslices*tid-1, numThreads_g, nslices);
                for j in (tid*nslices)..((tid+1)*nslices-1)-1 {
                    var transpose_j : int = p_transpose_idx(j, numThreads_g, nslices);
                    histogram_array[transpose_j] += histogram_array[transpose_j0];
                }
            }

            b.barrier();

            // now copy values into new structures (but not the mode we're sorting)
            for j_off in 0..(jend-jbegin)-1 {
                // we are actually going backwards
                var j : int = jend-j_off-1; 
                var idx : int = tt.ind[m,j];
                histogram[idx] -= 1;
                var offset : int = histogram[idx];
                new_vals[offset] = tt.vals[idx];
                for mode in 0..tt.nmodes-1 {
                    if mode != m {
                        new_ind[mode,offset] = tt.ind[mode,j];
                    }
                }
            }
        } /* coforall end */
        //writeln("histogram_array: ", histogram_array);
           
        //writeln("End coforall in counting_sort");
        for i in 0..tt.nmodes-1 {
            if i != m {
                tt.ind[i,0..tt.nnz-1] = new_ind[i,0..tt.nnz-1];
            }
        }
        tt.vals = new_vals;

        histogram_array[nslices] = tt.nnz;

        //writeln("Before doing quicksort on left over modes");
        // For 3/4D we can use quicksort on only the leftover modes
        var cmplt_ptr = c_ptrTo(cmplt);
        if tt.nmodes == 3 {
            //writeln("Before forall for quicksort2, nslices = ", nslices);
            for i in 0..nslices-1 {
                //writeln("Calling quicksort with start= ", histogram_array[i], " and end= ", histogram_array[i+1]);
                p_tt_quicksort2(tt, cmplt_ptr+1, histogram_array[i], histogram_array[i+1]);
                //writeln("  finished quicksort");
                for j in histogram_array[i]..histogram_array[i+1]-1 {
                    tt.ind[m,j] = i;
                }
                writeln("Finished with loop iter ", i);
            }
            writeln("End forall for quicksort2");
        }
        else if tt.nmodes == 4 {
            forall i in 0..nslices-1 {
                p_tt_quicksort3(tt, cmplt_ptr+1, histogram_array[i], histogram_array[i+1]);
                for j in histogram_array[i]..histogram_array[i+1]-1 {
                    tt.ind[m,j] = i;
                }
            }
        }
        else {
            // Shift cmplt left one time, then do normal quicksort
            var cmplt_ptr = c_ptrTo(cmplt);
            var saved = cmplt[0];
            c_memmove(cmplt_ptr, cmplt_ptr+1, (tt.nmodes-1)*8);
            cmplt[tt.nmodes-1] = saved;

            forall i in 0..nslices-1 {
                p_tt_quicksort(tt, cmplt, histogram_array[i], histogram_array[i+1]);
                for j in histogram_array[i]..histogram_array[i+1]-1 {
                    tt.ind[m,j] = i;
                }
            }
            // undo cmplt changes
            saved = cmplt[tt.nmodes-1];
            c_memmove(cmplt_ptr+1, cmplt_ptr, (tt.nmodes-1)*8);
            cmplt[0] = saved;
        }
        writeln("End of counting sort");
    }

    /*########################################################################
    #   Descriptipn:    Fill in later
    #
    #   Parameters:     
    #
    #   Return:         None
    ########################################################################*/
    private proc p_transpose_idx(idx: int, dim1: int, dim2: int) : int
    {
        return idx % dim1*dim2 + idx/dim1;
    }

    /*########################################################################
    #   Descriptipn:    Perfrom quicksort on a 2-mode tensor between start
    #                   and end.
    #
    #   Parameters:     tt (sptensor_t):    The tensor to sort
    #                   cmplt (int[]):      Mode permutation used for defining
    #                                       tie-breaking order
    #                   start (int):        First nonzero to sort
    #                   end (int):          Last nonzero to sort
    #
    #   Return:         None
    ########################################################################*/
    private proc p_tt_quicksort2(tt : sptensor_t, cmplt, start : int, end : int)
    {
        writeln("** Entered quicksort **");
        var vmid : real;
        var imid : [0..1] int;

        var ind0 = c_ptrTo(tt.ind[cmplt[0], 0..tt.nnz-1]);
        var ind1 = c_ptrTo(tt.ind[cmplt[1], 0..tt.nnz-1]);
        var vals = c_ptrTo(tt.vals);

        if (end-start) <= MIN_QUICKSORT_SIZE {
            p_tt_insertionsort2(tt, cmplt, start, end);
        }
        else {
            var i = start+1;
            var j = end-1;
            var k = start + ((end-start)/2);

            // grab pivot
            vmid = vals[k];
            vals[k] = vals[start];
            imid[0] = ind0[k];
            imid[1] = ind1[k];
            ind0[k] = ind0[start];
            ind1[k] = ind1[start];

            while i < j {
                // if tt[i] > mid --> tt[i] is on wrong side
                if p_ttqcmp2(ind0, ind1, i, imid) == 1 {
                    // if tt[j] <= mid --> swap tt[i] and tt[j]
                    if p_ttqcmp2(ind0, ind1, j, imid) < 1 {
                        var vtmp = vals[i];
                        vals[i] = vals[j];
                        vals[j] = vtmp;
                        var itmp = ind0[j];
                        ind0[i] = ind0[j];
                        ind0[j] = itmp;
                        itmp = ind1[i];
                        ind1[i] = ind1[j];
                        ind1[j] = itmp;
                        i += 1;
                    }
                    j -= 1;
                }
                else {
                    // if tt[j] > mid --> tt[j] is on right side
                    if p_ttqcmp2(ind0, ind1, j, imid) == 1 {
                        j -= 1;
                    }
                    i += 1;
                }
            }

            // if tt[i] > mid
            if p_ttqcmp2(ind0, ind1, i, imid) == 1 {
                i -= 1;
            }
            vals[start] = vals[i];
            vals[i] = vmid;
            ind0[start] = ind0[i];
            ind1[start] = ind1[i];
            ind0[i] = imid[0];
            ind1[i] = imid[1];
            if i > start+1 {
                p_tt_quicksort2(tt, cmplt, start, i);
            }
            i += 1; // skip pivot element
            if end -i > 1 {
                p_tt_quicksort2(tt, cmplt, i, end);
            }
        }
        writeln("** Exit quicksort **");
    }

    /*########################################################################
    #   Descriptipn:    Perfrom quicksort on a 3-mode tensor between start
    #                   and end.
    #
    #   Parameters:     tt (sptensor_t):    The tensor to sort
    #                   cmplt (int[]):      Mode permutation used for defining
    #                                       tie-breaking order
    #                   start (int):        First nonzero to sort
    #                   end (int):          Last nonzero to sort
    #
    #   Return:         None
    ########################################################################*/
    private proc p_tt_quicksort3(tt : sptensor_t, cmplt, start : int, end : int)
    {
        var vmid : real;
        var imid : [0..2] int;

        var ind0 = c_ptrTo(tt.ind[cmplt[0], 0..tt.nnz-1]);
        var ind1 = c_ptrTo(tt.ind[cmplt[1], 0..tt.nnz-1]);
        var ind2 = c_ptrTo(tt.ind[cmplt[2], 0..tt.nnz-1]);
        var vals = c_ptrTo(tt.vals);

        if (end-start) <= MIN_QUICKSORT_SIZE {
            p_tt_insertionsort3(tt, cmplt, start, end);
        }
        else {
            var i = start+1;
            var j = end-1;
            var k = start + ((end-start)/2);

            // grab pivot
            vmid = vals[k];
            vals[k] = vals[start];
            imid[0] = ind0[k];
            imid[1] = ind1[k];
            imid[2] = ind2[k];
            ind0[k] = ind0[start];
            ind1[k] = ind1[start];
            ind2[k] = ind2[start];
        
            while i < j {
                // if tt[i] > mid --> tt[i] is on wrong side
                if p_ttqcmp3(ind0, ind1, ind2, i, imid) == 1 {
                    // if tt[j] <= mid --> swap tt[i] and tt[j]
                    if p_ttqcmp3(ind0, ind1, ind2, j, imid) < 1 {
                        var vtmp = vals[i];
                        vals[i] = vals[j];
                        vals[j] = vtmp;
                        var itmp = ind0[j];
                        ind0[i] = ind0[j];
                        ind0[j] = itmp;
                        itmp = ind1[i];
                        ind1[i] = ind1[j];
                        ind1[j] = itmp;
                        itmp = ind2[i];
                        ind2[i] = ind2[j];
                        ind2[j] = itmp;
                        i += 1;
                    }
                    j -= 1;
                }
                else {
                    // if tt[j] > mid --> tt[j] is on right side
                    if p_ttqcmp3(ind0, ind1, ind2, j, imid) == 1 {
                        j -= 1;
                    }
                    i += 1;
                }
            }

            // if tt[i] > mid
            if p_ttqcmp3(ind0, ind1, ind2, i, imid) == 1 {
                i -= 1;
            }
            vals[start] = vals[i];
            vals[i] = vmid;
            ind0[start] = ind0[i];
            ind1[start] = ind1[i];
            ind2[start] = ind2[i];
            ind0[i] = imid[0];
            ind1[i] = imid[1];
            ind2[i] = imid[2];
            if i > start+1 {
                p_tt_quicksort3(tt, cmplt, start, i);
            }
            i += 1; // skip pivot element
            if end -i > 1 {
                p_tt_quicksort3(tt, cmplt, i, end);
            }    
        }
    }

    /*########################################################################
    #   Descriptipn:    Perfrom quicksort on a n-mode tensor between start
    #                   and end.
    #
    #   Parameters:     tt (sptensor_t):    The tensor to sort
    #                   cmplt (int[]):      Mode permutation used for defining
    #                                       tie-breaking order
    #                   start (int):        First nonzero to sort
    #                   end (int):          Last nonzero to sort
    #
    #   Return:         None
    ########################################################################*/
    private proc p_tt_quicksort(tt : sptensor_t, cmplt, start : int, end : int)
    {
        var vmid : real;
        var imid : [NUM_MODES_d] int;

        var ind = c_ptrTo(tt.ind[0, 0..tt.nnz-1]);
        var vals = c_ptrTo(tt.vals);
        var nmodes = tt.nmodes;

        if (end-start) <= MIN_QUICKSORT_SIZE {
            p_tt_insertionsort(tt, cmplt, start, end);
        }
        else {
            var i = start+1;
            var j = end-1;
            var k = start + ((end-start)/2);

            // grab pivot
            vmid = vals[k];
            vals[k] = vals[start];
            for m in 0..nmodes-1 {
                ind = c_ptrTo(tt.ind[m, 0..tt.nnz-1]);
                imid[m] = ind[k];
                ind[k] = ind[start];
            }
            while i < j {
                if p_ttqcmp(tt, cmplt, i, imid) == 1 {
                    if p_ttqcmp(tt, cmplt, j, imid) < 1 {
                        p_ttswap(tt, i, j);
                        i += 1;
                    }
                    j -= 1;
                }
                else {
                    if p_ttqcmp(tt, cmplt, j, imid) == 1 {
                        j -= 1;
                    }
                    i += 1;
                }
            }
            if p_ttqcmp(tt, cmplt, i, imid) == 1 {
                i -= 1;
            }
            vals[start] = vals[i];
            vals[i] = vmid;
            for m in 0..nmodes-1 {
                ind = c_ptrTo(tt.ind[m,0..tt.nnz-1]);
                ind[start] = ind[i];
                ind[i] = imid[m];
            }
            if i > start+1 {
                p_tt_quicksort(tt, cmplt, start, i);
            }
            i += 1;
            if end - i > 1 {
                p_tt_quicksort(tt, cmplt, i, end);
            }
        }
    }

    /*########################################################################
    #   Descriptipn:    Perfrom insertion sort on a 2-mode tensor between start
    #                   and end.
    #
    #   Parameters:     tt (sptensor_t):    The tensor to sort
    #                   cmplt (int[]):      Mode permutation used for defining
    #                                       tie-breaking order
    #                   start (int):        First nonzero to sort
    #                   end (int):          Last nonzero to sort
    #
    #   Return:         None
    ########################################################################*/
    private proc p_tt_insertionsort2(tt : sptensor_t, cmplt, start : int, end : int)
    {
        var ind0 = c_ptrTo(tt.ind[cmplt[0], 0..tt.nnz-1]);
        var ind1 = c_ptrTo(tt.ind[cmplt[1], 0..tt.nnz-1]);
        var vals = c_ptrTo(tt.vals);
        var vbuf : real;
        var ibuf : int;

        for i in start+1..end-1 {
            var j = i;
            while j > start && p_ttcmp2(ind0, ind1, i, j-1) < 0 {
                j -= 1;
            }
    
            vbuf = vals[i];
            // Shift all data
            c_memmove(vals+j+1, vals+j, (i-j)*8);
            vals[j] = vbuf;
            ibuf = ind0[i];
            c_memmove(ind0+j+1, ind0+j, (i-j)*8);
            ind0[j] = ibuf;
            ibuf = ind1[i];
            c_memmove(ind1+j+1, ind1+j, (i-j)*8);
            ind1[j] = ibuf;
        }
    }

    /*########################################################################
    #   Descriptipn:    Perfrom insertion sort on a 3-mode tensor between start
    #                   and end.
    #
    #   Parameters:     tt (sptensor_t):    The tensor to sort
    #                   cmplt (int[]):      Mode permutation used for defining
    #                                       tie-breaking order
    #                   start (int):        First nonzero to sort
    #                   end (int):          Last nonzero to sort
    #
    #   Return:         None
    ########################################################################*/
    private proc p_tt_insertionsort3(tt : sptensor_t, cmplt, start : int, end : int)
    {
        var ind0 = c_ptrTo(tt.ind[cmplt[0], 0..tt.nnz-1]);
        var ind1 = c_ptrTo(tt.ind[cmplt[1], 0..tt.nnz-1]);
        var ind2 = c_ptrTo(tt.ind[cmplt[2], 0..tt.nnz-1]);
        var vals = c_ptrTo(tt.vals);
        var vbuf : real;
        var ibuf : int;

        for i in start+1..end-1 {
            var j = i;
            while j > start && p_ttcmp3(ind0, ind1, ind2, i, j-1) < 0 {
                j -= 1;
            }

            vbuf = vals[i];
            // Shift all data
            c_memmove(vals+j+1, vals+j, (i-j)*8);
            vals[j] = vbuf;
            ibuf = ind0[i];
            c_memmove(ind0+j+1, ind0+j, (i-j)*8);
            ind0[j] = ibuf;
            ibuf = ind1[i];
            c_memmove(ind1+j+1, ind1+j, (i-j)*8);
            ind1[j] = ibuf;
            ibuf = ind2[i];
            c_memmove(ind2+j+1, ind2+j, (i-j)*8);
            ind2[j] = ibuf;
        }
    }

    /*########################################################################
    #   Descriptipn:    Perfrom insertion sort on a n-mode tensor between start
    #                   and end.
    #
    #   Parameters:     tt (sptensor_t):    The tensor to sort
    #                   cmplt (int[]):      Mode permutation used for defining
    #                                       tie-breaking order
    #                   start (int):        First nonzero to sort
    #                   end (int):          Last nonzero to sort
    #
    #   Return:         None
    ########################################################################*/
    private proc p_tt_insertionsort(tt : sptensor_t, cmplt, start : int, end : int)
    {
        var ind = c_ptrTo(tt.ind[0, 0..tt.nnz-1]);
        var vals = c_ptrTo(tt.vals);
        var nmodes = tt.nmodes;
        var vbuf : real;
        var ibuf : int;

        for i in start+1..end-1 {
            var j = i;
            while j > start && p_ttcmp(tt, cmplt, i, j-1) < 0 {
                j -= 1;
            }

            vbuf = vals[i];
            // Shift all data
            c_memmove(vals+j+1, vals+j, (i-j)*8);
            vals[j] = vbuf;
            for m in 0..nmodes-1 {
                ind = c_ptrTo(tt.ind[m, 0..tt.nnz-1]);
                ibuf = ind[i];
                c_memmove(ind+j+1, ind+j, (i-j)*8);
                ind[j] = ibuf;
            }
        }
    }

    /*########################################################################
    #   Descriptipn:    Compares ind*[i] and j[*] for 2-mode tensors
    #
    #   Parameters:     ind0 (int[]):   Primary mode to compare
    #                   ind1 (int[]):   Secondary mode to compare
    #                   i (int):        Index into ind*
    #                   j (int):        Second index into ind*
    #
    #   Return:         -1 if ind[i] < ind[j], 1 if ind[i] > ind[j], 0 if equal
    ########################################################################*/
    private proc p_ttcmp2(ind0, ind1, i : int, j: int)
    {
        if ind0[i] < ind0[j] {
            return -1;
        }
        else if ind0[j] < ind0[i] {
            return 1;
        }
        if ind1[i] < ind1[j] {
            return -1;
        }
        else if ind1[j] < ind1[i] {
            return 1;
        }
        return 0;
    }

    /*########################################################################
    #   Descriptipn:    Compares ind*[i] and j[*] for 3-mode tensors
    #
    #   Parameters:     ind0 (int[]):   Primary mode to compare
    #                   ind1 (int[]):   Secondary mode to compare
    #                   ind2 (int[]):   Final tie-breaking mode
    #                   i (int):        Index into ind*
    #                   j (int):        Second index into ind*
    #
    #   Return:         -1 if ind[i] < ind[j], 1 if ind[i] > ind[j], 0 if equal
    ########################################################################*/
    private proc p_ttcmp3(ind0, ind1, ind2, i : int, j: int)
    {
        if ind0[i] < ind0[j] {
            return -1;
        }
        else if ind0[j] < ind0[i] {
            return 1;
        }
        if ind1[i] < ind1[j] {
            return -1;
        }
        else if ind1[j] < ind1[i] {
            return 1;
        }
        if ind2[i] < ind2[j] {
            return -1;
        }
        else if ind2[j] < ind2[i] {
            return 1;
        }
        return 0;
    }

    /*########################################################################
    #   Descriptipn:    Compares ind*[i] and j[*] for n-mode tensors
    #
    #   Parameters:     tt (sptensor_t):    Tensor we are sorting
    #                   cmplt (int[]):      Mode permutation 
    #                   i (int):            Index into ind*
    #                   j (int):            Second index into ind*
    #
    #   Return:         -1 if ind[i] < ind[j], 1 if ind[i] > ind[j], 0 if equal
    ########################################################################*/
    private proc p_ttcmp(tt: sptensor_t, cmplt, i : int, j: int)
    {
        for m in 0..tt.nmodes-1 {
            if tt.ind[cmplt[m], i] < tt.ind[cmplt[m], j] {
                return -1;
            }
            else if tt.ind[cmplt[m], j] < tt.ind[cmplt[m], i] {
                return 1;
            }
        }
        return 0;
    }

    /*########################################################################
    #   Descriptipn:    Compares ind*[i] and j[*] for 2-mode tensors
    #
    #   Parameters:     ind0 (int[]):   Primary mode to compare
    #                   ind1 (int[]):   Secondary mode to compare
    #                   i (int):        Index into ind*[]
    #                   j (int[2]):     Indices we are comparing i against
    #
    #   Return:         -1 if ind[i] < j, 1 if ind[i] > j, 0 if equal
    ########################################################################*/
    private proc p_ttqcmp2(ind0, ind1, i : int, j)
    {
        if ind0[i] < j[0] {
            return -1;
        }
        else if j[0] < ind0[i] {
            return 1;
        }
        if ind1[i] < j[1] {
            return -1;
        }
        else if j[1] < ind1[i] {
            return 1;
        }
        return 0;
    }

    /*########################################################################
    #   Descriptipn:    Compares ind*[i] and j[*] for 3-mode tensors
    #
    #   Parameters:     ind0 (int[]):   Primary mode to compare
    #                   ind1 (int[]):   Secondary mode to compare
    #                   ind2 (int[]):   Final tie-breaking mode
    #                   i (int):        Index into ind*[]
    #                   j (int[3]):     Indices we are comparing i against
    #
    #   Return:         -1 if ind[i] < j, 1 if ind[i] > j, 0 if equal
    ########################################################################*/
    private proc p_ttqcmp3(ind0, ind1, ind2, i : int, j)
    {
        if ind0[i] < j[0] {
            return -1;
        }
        else if j[0] < ind0[i] {
            return 1;
        }
        if ind1[i] < j[1] {
            return -1;
        }
        else if j[1] < ind1[i] {
            return 1;
        }
        if ind2[i] < j[2] {
            return -1;
        }
        else if j[2] < ind2[i] {
            return 1;
        }
        return 0;
    }

    /*########################################################################
    #   Descriptipn:    Compares ind*[i] and j[*] for n-mode tensors
    #
    #   Parameters:     tt (sptensor_t):    Tensor we are sorting
    #                   cmplt (int[]):      Mode permutation 
    #                   i (int):            Index into ind*
    #                   j (int[]):          Coordinate we are comparing against
    #
    #   Return:         -1 if ind[i] < ind[j], 1 if ind[i] > ind[j], 0 if equal
    ########################################################################*/
    private proc p_ttqcmp(tt: sptensor_t, cmplt, i : int, j)
    {
        for m in 0..tt.nmodes-1 {
            if tt.ind[cmplt[m], i] < j[cmplt[m]] {
                return -1;
            }
            else if j[cmplt[m]] < tt.ind[cmplt[m], i] {
                return 1;
            }
        }
        return 0;
    }

    /*########################################################################
    #   Descriptipn:    Swap nonzeros i and j
    #
    #   Parameters:     tt (sptensor_t):    Tensor we are sorting
    #                   i (int):            First nonzero to swap
    #                   j (int[]):          Second nonzero to swap
    #
    #   Return:         -1 if ind[i] < ind[j], 1 if ind[i] > ind[j], 0 if equal
    ########################################################################*/
    private proc p_ttswap(tt: sptensor_t, i : int, j)
    {
        var vtmp = tt.vals[i];
        tt.vals[i] = tt.vals[j];
        tt.vals[j] = vtmp;
        var itmp : int;
        for m in 0..tt.nmodes-1 {
            itmp = tt.ind[m,i];
            tt.ind[m,i] = tt.ind[m,j];
            tt.ind[m,j] = itmp;
        }
    }
}
