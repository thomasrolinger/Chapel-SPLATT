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
    use MutexPool;

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
    *   Functions for MTTKRP-root locked
    *
    ******************************/
    class mttkrp_root_locked {
        // Just a container for a function so we can "pass" the
        // function as an argument
        
        // ct is ONE of the CSF tensors; tid is the task ID that called this
        // function.
        proc p_csf_mttkrp_root_locked(ct, tile_id, mats, mode, thds, partition, tid)
        {
            // extract tensor structures
            var nmodes = ct.nmodes;
            var vals = ct.pt[tile_id].vals;

            // empty tile?? I don't think this will ever be the case for us
            
            if nmodes == 3 {
                p_csf_mttkrp_root3_locked(ct, tile_id, mats, mode, thds, partition, tid);
                return;
            }

            var fp = ct.pt[tile_id].fptr;
            var fids = ct.pt[tile_id].fids;
            var nfactors = mats[0].J;

            // mvals is an array of nmodes double pointers.
            var mvals : [NUM_MODES_d] c_ptr(real);
            // buf is an array of nmodes double pointers
            var buf : [NUM_MODES_d] c_ptr(real);
            var idxstack : [NUM_MODES_d] int;
            
            for m in 0..nmodes-1 {
                mvals[m] = c_ptrTo(mats[csf_depth_to_mode(ct, m)].vals);
                // grab next row of buf from thds
                buf[m] = c_ptrTo(thds[tid].scratch[2].buf[nfactors*m]);
                c_memset(buf[m], 0, nfactors*8);
            }
            
            // ovals is a pointer to the 2D matrix vals. If we want to
            // access any elements in vals, we need to use pointer arithmetic
            // and assume the values are in row major order: to access ovals[i][j],
            // we need to say ovals[(i*J)+j] where J is the number of columns in vals
            var ovals = c_ptrTo(mats[nmodes].vals);

            var nfibs = ct.pt[tile_id].nfibs[0];
            assert(nfibs <= mats[nmodes].I);

            var start = partition[tid];
            var stop = partition[tid+1];
            for s in start..stop-1 {
                // checking for -1 is the same as checking for NULL (for us).
                var fid = if fids[0].fiber_ids[0] == -1 then s else fibs[0].fiber_ids[s];
                assert(fid < mats[nmodes].I);
                //TODO: Implement the function below
                //p_propagate_up(buf[0], buf, idxstack, 0, s, fp, fids, vals, mvals, nmodes, nfactors);
                
                // These are 1D arrays.
                var orow = ovals(fid*nfactors);
                var obuf = buf[0];
                mutex_set_lock(pool_g, fid);
                for f in 0..nfactors-1 {
                    orow[f] += obuf[f];
                }
                mutex_unset_lock(pool_g, fid);
            }
            
        }
    }
    
    //***********************************************************************
    private proc p_csf_mttkrp_root3_locked(ct, tile_id, mats, mode, thds, partition, tid)
    {
        assert(ct.nmodes == 3);
        var nmodes = ct.nmodes;
        // pointers to 1D chapel arrays
        var vals = c_ptrTo(ct.pt[tile_id].vals);
        var sptr = c_ptrTo(ct.pt[tile_id].fptr[0].subtree);
        var fptr = c_ptrTo(ct.pt[tile_id].fptr[1].subtree);
        var sids = c_ptrTo(ct.pt[tile_id].fids[0].fiber_ids);
        var fids = c_ptrTo(ct.pt[tile_id].fids[1].fiber_ids);
        var inds = c_ptrTo(ct.pt[tile_id].fids[2].fiber_ids);

        // pointers to 2D chapel matrices
        var avals = c_ptrTo(mats[csf_depth_to_mode(ct,1)].vals);
        var bvals = c_ptrTo(mats[csf_depth_to_mode(ct,2)].vals);
        var ovals = c_ptrTo(mats[nmodes].vals);

        var nfactors = mats[nmodes].J;

        // pointer to 1D chapel array
        var accumF = c_ptrTo(thds[tid].scratch[0].buf);
        
        // write to output
        var writeF = c_ptrTo(thds[tid].scratch[2].buf);
        for r in 0..nfactors-1 {
            writeF[r] = 0.0;
        }

        var nslices = ct.pt[tile_id].nfibs[0];
        var start = partition[tid];
        var stop = partition[tid+1];
        for s in start..stop-1 {
            /* for each fiber in slice */   
            for f in sptr[s]..sptr[s+1]-1 {
                /* first entry of the fiber is used to initialize accumF */
                var jjfirst = fptr[f];
                var vfirst = vals[jjfirst];
                var bv = bvals + (inds[jjfirst] * nfactors);
                for r in 0..nfactors-1 {
                    accumF[r] = vfirst * bv[r];
                }
                /* for each nnz in fiber */
                for jj in fptr[f]+1..fptr[f+1]-1 {
                    var v = vals[jj];
                    var bv = bvals + (inds[jj] * nfactors);
                    for r in 0..nfactors-1 {
                        accumF[r] += v * bv[r];
                    }
                }
                /* scale inner products by row of A and upate to M */
                var av = avals + (fids[f] * nfactors);
                for r in 0..nfactors-1 {
                    writeF[r] += acumF[r] * av[r];
                }
            }
            // checking for sids == NULL
            var fid = if sids[0] == -1 then s else sids[s];
            var mv = ovals + (fid * nfactors);
            
            /* flush to output */
            mutex_set_lock(pool_g, fid);
            for r in 0.. nfactors-1 {
                mv[r] += writeF[r];
                writeF[r] = 0.0;
            }
            mutex_unset_lock(pool_g, fid);
        }
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
    #                   thds (thd_info[]):      Thread scratch space
    #                   ws (splatt_mttkrp_ws):  MTTKRP workspace
    #                   opts (cpd_cmd_args):    Arguments.
    #
    #   Return:         None
    ########################################################################*/
    proc mttkrp_csf(tensors, mats, mode : int, thds, ws : splatt_mttkrp_ws, opts)
    {
        // Set up mutex pool. The pool itself is global
        pool_g = mutex_alloc(DEFAULT_NLOCKS, DEFAULT_LOCK_PAD);
        
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

    /*########################################################################
    #   Descriptipn:    Map MTTKRP functions onto a possibly tiled CSF tensor.
    #                   This function will handle any scheduling required with
    #                   a partitially tiled tensor.
    #
    #   Parameters:     tensors (splatt_csf[]): An array of CSF representations.
    #                                           tensors[csf_id] is processed.
    #                   csf_id (int):           Which tensor are we processing
    #                   atomic_func:            An MTTKRP function which atomically
    #                                           updates the output.
    #                   nosync_func:            An MTTKRP function which does not
    #                                           atomically update.
    #                   mats:                   The matrices, output stored at mat[nmodes].
    #                   mode:                   Which mode of tensors is the output
    #                   ws:                     MTTKRP workspace.
    #
    #   Return:         None
    ########################################################################*/
    //***** NOTE: For now, these function pointers will be objects, which
    // contain the functions.
    private proc p_schedule_tiles(tensors, csf_id, atomic_func, nosync_func, mats,
                                  mode, ws)
    {
        return;
    }


}
