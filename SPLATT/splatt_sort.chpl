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
    proc tt_sort(tt : sptensor_t, mode : int, dim_perm : [NUM_MODES_d] int)
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
        if dim_perm == nil {
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
        var new_ind : [COORD_D] int;

        // Array with nnz elements
        var new_vals : [NNZ_d] real;

        /* This will be a bit difficult..... */
        //var histogram_array : 
    }
}
