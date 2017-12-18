/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/18/2017
    File:   mttkrp.chpl

    Description:    This is a module file for things related
                    to the MTTKRP
*/

module MTTKRP {
    use Base;
    use CSF;
    use Args;
    use Assert;
    use Util;

    /*****************************
    *
    *   priv_buffer class
    *
    ******************************/
    class priv_buffer {
        var buffer_d : domain(2) = {0..1, 0..1};
        var buffer : [buffer_d] real;
    }

    /*****************************
    *
    *   tree_partition class
    *
    ******************************/
    class tree_part {
        var tree_part_d : domain(1) = 0..1;
        var buf : [tree_part_d] int;
    }

    /*****************************
    *
    *   MTTKRP workspace class
    *
    ******************************/
    class splatt_mttkrp_ws {
        /* How many CSF representations are available */
        var num_csf : int;
        /* Mapping of modes -> CSF representations (which CSF to use) */  
        var mode_csf_map : [NUM_MODES_d] int;
        /* Number of threads to use */
        var num_threads : int;

        /*
            Partitioning information. If the CSF is tiled, we
            distribute tiles to threads. If the CSF is untiled,
            we distribute slices to threads. Partitioning is performed
            on a per-CSF basis, so we have one for each mode (the
            maximum number of CSF). In all cases, we rely on static
            partitioning via chains-on-chains partitioning.
        */
        var tree_d : domain(1) = 0..1;
        var tile_partition : [0..0] int;
        var tree_partition : [tree_d] tree_part;

        /*
            Privatization information: Privatizing a mode replicates the
            output matrix by each thread in order to avoid lock contention.
            This is useful when the tensor mode is short.
        */
        
        /* Marks if a tensor mode is privatized */
        var is_privatized : [NUM_MODES_d] bool;
        /* 
            The buffer used by each thread for privatization.
            privatize_buffer[thread_id] is large enough to process
            the largest privatized mode.
        */
        var privatize_buffer : [0..numThreads_g] priv_buffer;

        /* Time spent on the latest privatized reduction */
        var reduction_time : real;
    }

    /*****************************
    *
    *   Public Functions
    *
    ******************************/
    
    /*########################################################################
    #   Descriptipn:    Allocates splatt_mttkrp_ws
    #
    #   Parameters:     tensors (splatt_csf[]):     Tensor to factor in CSF
    #                   ncolumns (int):             Number of columns in factor
    #                                               matrices
    #                   args (cpd_cmd_args):        Arguments
    #
    #   Return:         splatt_mttkrp_ws class object
    ########################################################################*/
    proc splatt_mttkrp_alloc_ws(tensors, ncolumns : int, opts : cpd_cmd_args) : splatt_mttkrp_ws
    {
        var ws : splatt_mttkrp_ws = new splatt_mttkrp_ws();
        var num_csf : int = 0;
        ws.num_threads = numThreads_g;

        /* map each MTTKRP mode to a CSF tensor */
        for m in 0..tensors[0].nmodes-1 {
            select(opts.numCSF) {
                when "one" {
                    /* only one tensor, map is easy */
                    ws.mode_csf_map[m] = 0;
                    num_csf = 1;
                }
                when "two" {
                    /* last mode is mapped to 2nd tensor */
                    ws.mode_csf_map[m] = 0;
                    if csf_mode_to_depth(tensors[0], m) == tensors[0].nmodes-1 {
                        ws.mode_csf_map[m] = 1;
                    }
                    num_csf = 2;
                }
                when "all" {
                    /* each mode has its own tensor, map is easy */
                    ws.mode_csf_map[m] = m;
                    num_csf = tensors[0].nmodes;
                }
            }
        }
    
        assert(num_csf > 0);
        ws.num_csf = num_csf;

        /* 
            Setup partition info for each CSF.
            Since we aren't doing tiling yet, this will only
            affect the tree partition.

            Below, we use -1 to stand for nil/null
        */
        //ws.tree_part_d = {0..num_csf-1, 0..(ws.num_threads+1)-1};
        ws.tree_d = 0..num_csf-1;
        for c in 0..num_csf-1 {
            ws.tree_partition[c] = new tree_part();
            ws.tree_partition[c].buf = -1;
        }
        // We're not doing tiling, so we just have a "blank" variable for it
        ws.tile_partition = -1;
        for c in 0..num_csf-1 {
            var csf = tensors[c];
            writeln("csf[0].pt[0].fptr[0]:");
            for e in tensors[0].pt[0].fptr[0].subtree {
                writeln(e);
            }
            exit(-1);
            ws.tree_partition[c] = csf_partition_1d(csf, 0, ws.num_threads);
        }
        
    
        /* Allocate privatization buffer */
        var largest_priv_dim : int = 0;
        for m in 0..tensors[0].nmodes-1 {
            ws.is_privatized[m] = p_is_privatized(tensors, m, opts);
            if ws.is_privatized[m] == true {
                largest_priv_dim = max(largest_priv_dim, tensors[0].dims[m]);
                writeln("# PRIVATIZING-MODE: ", m+1);
            }
        }
        for t in 0..ws.num_threads-1 {
            ws.privatize_buffer[t] = new priv_buffer();
            ws.privatize_buffer[t].buffer_d = {0..largest_priv_dim-1, 0..ncolumns-1};
        }
        if largest_priv_dim > 0 {
            var bytes = ws.num_threads * largest_priv_dim * ncolumns * 8;
            var bstr = bytes_str(bytes);
            writeln("# PRIVATIZATION-BUF: ", bstr);
            writeln("");
        }
        return ws;
    }

    /*########################################################################
    #   Descriptipn:    Matricized Tensor Times Khatri-Rao Product (MTTKRP)
    #                   with a CSF tensor. This is the primary computation
    #                   involved in CPD. Output is written to mats[nmodes].
    #
    #   Parameters:     tensors (splatt_csf[]): Tensor to factor in CSF
    #                   mats (dense_matrix[]):  Input/output matrices
    #                   mode (int):             Which mode we are computing for
    #                   ws (splatt_mttkrp_ws):  MTTKRP workspace
    #                   opts (cpd_cmd_args):    Arguments.
    #
    #   Return:         None
    ########################################################################*/
    proc mttkrp_csf(tensors, mats, mode : int, ws : splatt_mttkrp_ws, opts)
    {
        // Mutex pool???
        
        var nmodes = tensors[0].nmodes;

        /* Clear output matrix */
        var M = mats[nmodes];
        M.I = tensors[0].dims[mode];
        M.vals = 0;

        /* Choose which MTTKRP function to use */
        var which_csf = ws.mode_csf_map[mode];
        var outdepth = csf_mode_to_depth(tensors[which_csf], mode);
        if outdepth == 0 {
            /* root */
        }
    }

    /*****************************
    *
    *   Private Functions
    *
    ******************************/
    
    /*########################################################################
    #   Descriptipn:    Should a certain mode be privatized to avoid locks?
    #
    #   Parameters:     csf (splatt_csf[]): Just used for dimensions
    #                   mode (int):         The mode we're processing
    #                   opts (cpd_cmd_args):Arguments
    #
    #   Return:         bool: True if we should privatize
    ########################################################################*/
    private proc p_is_privatized(csf, mode : int, opts : cpd_cmd_args) : bool
    {
        var length : int = csf[0].dims[mode];
        var nthreads : int = opts.numThreads;
        var thresh : real = DEFAULT_PRIV_THRESH;
        // Don't bother if we're not multi-threaded
        if nthreads == 1 {
            return false;
        }
        var a = (length * nthreads) : real;
        var b = (thresh * csf[0].nnz:real);
        return a <= b;
    }
}
