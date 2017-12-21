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
    use Barriers;

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

/*#################################################################################################
#
#   Below are the classes for the various MTTKRP functions, which we need to "pass" into the
#   p_schedule_tiles function. We declare a class for each function, and the only thing within
#   the class is the function, called func. The function within the class will call out to other
#   functions, which are defined further below.
#
####################################################################################################*/


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

        //p_csf_mttkrp_root_locked
        proc func(ct, tile_id, mats, mode, thds, partition, tid)
        {
            // extract tensor structures
            var nmodes = ct.nmodes;
            var vals = c_ptrTo(ct.pt[tile_id].vals);

            // empty tile?? I don't think this will ever be the case for us
            
            if nmodes == 3 {
                p_csf_mttkrp_root3_locked(ct, tile_id, mats, mode, thds, partition, tid);
                return;
            }

            // These get references to fptr and fids, so if
            // we change something in fp and fids, it will
            // change the values within ct as well.
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
                var fid = if fids[0].fiber_ids[0] == -1 then s else fids[0].fiber_ids[s];
                assert(fid < mats[nmodes].I);
                p_propagate_up(buf[0], buf, idxstack, 0, s, fp, fids, vals, mvals, nmodes, nfactors);
                
                // These are 1D arrays.
                var orow = ovals + (fid*nfactors);
                var obuf = buf[0];
                mutex_set_lock(pool_g, fid);
                for f in 0..nfactors-1 {
                    orow[f] += obuf[f];
                }
                mutex_unset_lock(pool_g, fid);
            }
        }/* end func */
    }/* end mttkrp_root_locked class */
 
    /*****************************
    *
    *   Functions for MTTKRP-root nolock
    *
    ******************************/
    class mttkrp_root_nolock {
        //p_csf_mttkrp_root_nolock
        proc func(ct, tile_id, mats, mode, thds, partition, tid)
        {
            /* extract tensor structures */
            var nmodes = ct.nmodes;
            var vals = c_ptrTo(ct.pt[tile_id].vals);

            // empty tile?? I don't think this will ever be the case for us

            if nmodes == 3 {
                p_csf_mttkrp_root3_nolock(ct, tile_id, mats, mode, thds, partition, tid);
                return;
            }
        
            // These get references to fptr and fids, so if
            // we change something in fp and fids, it will
            // change the values within ct as well.
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
                var fid = if fids[0].fiber_ids[0] == -1 then s else fids[0].fiber_ids[s];
                assert(fid < mats[nmodes].I);
                p_propagate_up(buf[0], buf, idxstack, 0, s, fp, fids, vals, mvals, nmodes, nfactors);

                // These are 1D arrays.
                var orow = ovals + (fid*nfactors);
                var obuf = buf[0];
                mutex_set_lock(pool_g, fid);
                for f in 0..nfactors-1 {
                    orow[f] += obuf[f];
                }
                mutex_unset_lock(pool_g, fid);
            }
        }/* end func */
    }/* end mttkrp_root_nolock class */

    /*****************************
    *
    *   Functions for MTTKRP-leaf locked
    *
    ******************************/
    class mttkrp_leaf_locked {
        // p_mttkrp_leaf_locked
        proc func(ct, tile_id, mats, mode, thds, partition, tid)
        {
            /* extract tensor structures */
            var nmodes = ct.nmodes;
            var vals = c_ptrTo(ct.pt[tile_id].vals);

            // empty tile?? I don't think this will ever be the case for us

            if nmodes == 3 {
                p_csf_mttkrp_leaf3_locked(ct, tile_id, mats, mode, thds, partition, tid);
                return;
            }
            
            // These get references to fptr and fids, so if
            // we change something in fp and fids, it will
            // change the values within ct as well.
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
            }

            /* for each outer slice */
            var nslices = ct.pt[tile_id].nfibs[0];
            var start = partition[tid];
            var stop = partition[tid+1];
            for s in start..stop-1 {
                // checking for -1 is the same as checking for NULL (for us).
                var fid = if fids[0].fiber_ids[0] == -1 then s else fids[0].fiber_ids[s];
                idxstack[0] = s;
                /* clear out stale data */
                for m in 1..(nmodes-1)-1 {
                    idxstack[m] = fp[m-1].subtree[idxstack[m-1]];
                }
                /* first buf will always just be a matrix row */
                var rootrow = mvals[0] + (fid*nfactors);
                var rootbuf = buf[0];
                for f in 0..nfactors-1 {
                    rootbuf[f] = rootrow[f];
                }
                var depth = 0;
                var outer_end = fp[0].subtree[s+1];
                while idxstack[1] < outer_end {
                    /* move down to an nnz node */
                    while depth < nmodes-2 {
                        /* propagate buf down */
                        var drow = mvals[depth+1] + (fids[depth+1].fiber_ids[idxstack[depth+1]] * nfactors);
                        p_assign_hada(buf[depth+1], buf[depth], drow, nfactors);
                        depth += 1
                    }
                    /* process all nonzeros [start, end) */
                    var start2 = fp[depth].subtree[idxstack[depth]];
                    var end = fp[depth].subtree[idxstack[depth]+1];
                    p_csf_process_fiber_locked(c_ptrTo(mats[nmodes].vals), buf[depth], nfactors, start2, end, fids[depth+1], vals);

                    /* now move back up to the next unprocessed child */
                    do {
                        idxstack[depth] += 1;
                        depth -= 1;
                    } while depth > 0 && idxstack[depth+1] == fp[depth].subtree[idxstack[depth]+1];
                } /* end DFS */
            } /* end outer slice loop */
        }/* end func */
    } /* end mttkrp_leaf_locked class */

    /*****************************
    *
    *   Functions for MTTKRP-leaf no lock
    *
    ******************************/
    class mttkrp_leaf_nolock {
        // p_csf_mttkrp_leaf_nolock
        proc func(ct, tile_id, mats, mode, thds, partition, tid)
        {
            /* extract tensor structures */
            var nmodes = ct.nmodes;
            var vals = c_ptrTo(ct.pt[tile_id].vals);

            // empty tile?? I don't think this will ever be the case for us

            if nmodes == 3 {
                p_csf_mttkrp_leaf3_nolock(ct, tile_id, mats, mode, thds, partition, tid);
                return;
            }
        } /* end func */
    } /* end mttkrp_leaf_nolock class */
   
/*#################################################################################################
#
#   Below are the various private functions
#
####################################################################################################*/

    /*****************************
    *
    *   Private Functions
    *
    ******************************/
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
                    writeF[r] += accumF[r] * av[r];
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

    //***********************************************************************
    private proc p_csf_mttkrp_root3_nolock(ct, tile_id, mats, mode, thds, partition, tid)
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
            var fid = if sids[0] == -1 then s else sids[s];
            var mv = ovals + (fid * nfactors);
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
                    writeF[r] += accumF[r] * av[r];
                }
            }
            /* flush to output */
            for r in 0.. nfactors-1 {
                mv[r] += writeF[r];
                writeF[r] = 0.0;
            }
        }
    }
    
    //***********************************************************************
    private proc p_csf_mttkrp_leaf3_locked(ct, tile_id, mats, mode, thds, partition, tid)
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
        var avals = c_ptrTo(mats[csf_depth_to_mode(ct,0)].vals);
        var bvals = c_ptrTo(mats[csf_depth_to_mode(ct,1)].vals);
        var ovals = c_ptrTo(mats[nmodes].vals);
        var nfactors = mats[nmodes].J;

        // pointer to 1D chapel array
        var accumF = c_ptrTo(thds[tid].scratch[0].buf);

        var nslices = ct.pt[tile_id].nfibs[0];
        var start = partition[tid];
        var stop = partition[tid+1];

        for s in start..stop-1 {
            var fid = if sids[0] == -1 then s else sids[s];
            /* root row */
            var rv = avals + (fids * nfactors);
            /* for each fiber in slice */
            for f in sptr[s]..sptr[s+1]-1 {
                /* fill fiber with hada */
                var av = bvals + (fids[f] * nfactors);
                for r in 0..nfactors-1 {
                    accumF[r] = rv[r] * av[r];
                }
                /* for each nnz in fiber, scale with hada and write to ovals */
                for jj in fptr[f]..fptr[f+1]-1 {
                    var v = vals[jj];
                    var ov = ovals + (inds[jj] * nfactors);
                    mutex_set_lock(pool_g, inds[jj]);
                    for r in 0..nfactors-1 {
                        ov[r] += v * accumF[r];
                    }
                    mutex_unset_lock(pool_g, inds[jj]);
                }
            }
        }
    }

    //***********************************************************************
    /*
        outBuf :    cptr(real) -> 1D array of reals
        buf:        []cptr(real) -> 1D array of pointers to reals
        idxstack:   int[] -> Chapel array of ints
        init_depth: int
        init_idx:   int
        fp:         []csf_subtree -> array of objects
        fids:       []csf_fiber_ids -> array of objects
        vals:       cptr(real) -> 1D array of reals
        mvals:      []cptr(real) -> 1D array of pointers to reals
    */
    private proc p_propagate_up(outBuf, buf, idxstack, init_depth, init_idx, fp, 
                                fids, vals, mvals, nmodes, nfactors)
    {
        /* push initial idx initialize idxstack */
        idxstack[init_depth] = init_idx;
        for m in init_depth+1..nmodes-1 {
            idxstack[m] = fp[m-1].subtree[idxstack[m-1]];
        }
        assert(init_depth < nmodes-1);

        /* clear out accumulation buffer */
        for f in 0..nfactors-1 {
            buf[init_depth+1][f] = 0;
        }

        while idxstack[init_depth+1] < fp[init_depth].subtree[init_idx+1] {
            /* skip to last internal mode */
            var depth = nmodes - 2;
            /* process all nonzeros [start, end) into buf[depth] */
            var start = fp[depth].subtree[idxstack[depth]];
            var end = fp[depth].subtree[idxstack[depth]+1];
            p_csf_process_fiber(buf[depth+1], nfactors, mvals[depth+1], start, end,
                                fids[depth+1], vals);
            idxstack[depth+1] = end;

            /* exit early if there is no propagation to do */
            if init_depth == nmodes-2 {
                for f in 0..nfactors-1 {
                    outBuf[f] = buf[depth+1][f];
                }
                return;
            }

            /* propagate up until we reach a node with more children to process */
            do {
                /* propagate result up and clear buffer for next sibling */
                // fibrow is a ptr somewhere into mvals. So it's a 1D array
                var fibrow = mvals[depth] + (fids[depth].fiber_ids[idxstack[depth]] * nfactors);
                p_add_hada_clear(buf[depth], buf[depth+1], fibrow, nfactors);
                idxstack[depth] += 1;
                depth -= 1;
            } while depth > init_depth && idxstack[depth+1] == fp[depth].subtree[idxstack[depth]+1];
        }
        /* copy to out */       
        for f in 0..nfactors-1 {
            outBuf[f] = buf[init_depth+1][f];
        }
    }

    //***********************************************************************
    /*
        accumbuf:   real* -> 1D array of reals
        leafmat:    real* -> 1D array of reals
        inds:       csf_fiber_id -> object
        vals:       real* -> 1D array of reals
    */
    private proc p_csf_process_fiber(accumbuf, nfactors, leafmat, start, end,
                                     inds, vals)
    {
        /* for each nnz in fiber */
        for j in start..end-1 {
            var v = vals[j];
            var row = leafmat + (nfactors * inds.fiber_ids[j]);
            for f in 0..nfactors-1 {
                accumbuf[f] += v * row[f];
            }
        }
    }
    private proc p_csf_process_fiber_locked(leafmat, accumbuf, nfactors, start, end
                                            inds, val)
    {
        for jj in start..end-1 {
            var leafrow = leafmat + (inds.fiber_ids[jj] * nfactors);
            var v = vals[jj];
            mutex_set_lock(pool_g, inds.fiber_ids[jj]);
            for f in 0..nfactors-1 {
                leafrow[f] += v * accumbuf[f];
            }
            mutex_unset_lock(pool_g, inds.fiber_ids[jj]);
        }
    }

    //***********************************************************************
    /*
        outBuf: real*
        a:      real*
        b:      real*
    */
    private proc p_add_hada_clear(outBuf, a, b, nfactors)
    {
        for f in 0..nfactors-1 {
            outBuf[f] += a[f] * b[f];
            a[f] = 0;
        }
    }
    private proc p_assign_hada(outBuf, a, b, nfactors)
    {
        for f in 0..nfactors-1 {
            outBuf[f] = a[f] * b[f];
        }
    }

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
    #                   mttkrp_func:            An MTTKRP function
    #                   mats:                   The matrices, output stored at mat[nmodes].
    #                   mode:                   Which mode of tensors is the output
    #                   thds:                   Thread workspace
    #                   ws:                     MTTKRP workspace.
    #
    #   Return:         None
    ########################################################################*/
    private proc p_schedule_tiles(tensors, csf_id, mttkrp_func, mats, mode, thds, ws)
    {
        var csf = tensors[csf_id];
        var nmodes = csf.nmodes;
        var depth = nmodes-1;
        var nrows = mats[mode].I;
        var ncols = mats[mode].J;

        /* Store old pointer */
        var global_output = c_ptrTo(mats[nmodes].vals);
    
        /* barrier used in p_reduce_privatized */
        var b = new Barrier(numThreads_g);

        coforall tid in 0..numThreads_g-1 {
            // this is a buffer/array, not an object
            var tree_partition = ws.tree_partition[csf_id].buf;
            /*
                We may need to edit mats[nmodes].vals so create a
                private copy of the pointers to edit (not the entire
                factor matrix). Updating vals within mats_priv WILL
                change the vals in mats.
            */
            var mats_priv : [0..(nmodes+1)-1] dense_matrix;
            for m in 0..nmodes-1 {
                mats_priv[m] = mats[m];
            }
            /* 
                Each thread gets a separate structure, but do a shallow copy.
                This is where we use the vals_ref field. Need to make sure that
                when we want to modify vals within mats_priv[nmodes] that we
                use vals_ref.
            */
            mats_priv[nmodes] = new dense_matrix();
            mats_priv[nmodes].I = mats[nmodes].I;
            mats_priv[nmodes].J = mats[nmodes].J;
            mats_priv[nmodes].matrix_domain = mats[nmodes].matrix_domain;
            mats_priv[nmodes].vals_ref = mats_priv[nmodes].vals_ref;

            /* 
                Give each thread its own private buffer and overwrite
                atomic function.
            */
            if ws.is_privatized[mode] {
                /* change (thread-private!) output structure */
                var privBufPtr = c_ptrTo(ws.privatize_buffer[tid].buffer);
                c_memset(privBufPtr, 0, nrows*ncols*8);
                mats_priv[nmodes].vals_ref = privBufPtr;
            }

            /* We don't support tiling right now...*/
            if csf.ntiles > 1 {
                writeln("ERROR: Tiling not supported");
            }
            else {
                /* untiled, parallelize within kernel */
                mttkrp_func.func(csf, 0, mats_priv, mode, thds, tree_partition, tid);
            }

            /* If we used privatization, perform a reduction */
            if ws.is_privatized[mode] {
                p_reduce_privatized(ws, global_output, nrows, ncols, tid, b);
            }
        } /* end coforall */

        /* restore pointer */
        mats[nmodes].vals_ref = global_output;
    } 

    /*########################################################################
    #   Descriptipn:    Perform a reduction on thread-local MTTKRP output
    #
    #   Parameters:     ws:             MTTKRP workspace.
    #                   global_output:  The global MTTKRP output we are
    #                                   reducing into
    #                   nrows:          Number of rows in MTTKRP
    #                   ncols:          Number of cols in MTTKRP
    #                   tid:            This task's ID
    #                   b:              Barrier to use
    #
    #   Return:         None
    ########################################################################*/
    private proc p_reduce_privatized(ws, global_output, nrows, ncols, tid, b)
    {
        /* Ensure everyone has finished their local MTTKRP */
        b.barrier();
        var elem_per_thread = (nrows * ncols) / numThreads_g;
        var start = tid * elem_per_thread;
        var stop = if tid == numThreads_g - 1 then (nrows*ncols) else (tid+1)*elem_per_thread;
        
        /* reduction */
        for t in 0..numThreads_g-1 {
            var thread_buf = c_ptrTo(ws.privatize_buffer[t].buffer);
            for x in start..stop-1 {
                global_output[x] += thread_buf[x];
            }
        }
    }

/*#################################################################################################
#
#   Below are the public functions
#
####################################################################################################*/

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
            /*
                In SPLATT, function ptrs are passed into the
                p_schedule_tiles function. We create a class
                which has the function we need and pass that
                class instance in.
    
                If we are privatizing this mode, then 
            */
            var mttkrp_func = nil;
            if ws.is_privatized[mode] {
                /* don't use atomics */
                mttkrp_func = new mttkrp_root_nolock();
            }
            else {
                mttkrp_func = new mttkrp_root_locked();
            }
            p_schedule_tiles(tensors, which_csf, mttkrp_func, mats, mode, thds, ws);
        }
        else if outdepth == nmodes-1 {
            /* leaf */
            var mttkrp_func = nil;
            if ws.is_privatized[mode] {
                /* don't use atomics */
                mttkrp_func = new mttkrp_leaf_nolock();
            }
            else {
                mttkrp_func = new mttkrp_leaf_locked();
            }
            p_schedule_tiles(tensors, which_csf, mttkrp_func, mats, mode, thds, ws);
        }
        else {
            //TODO: Do functions for internal nodes
            writeln("Coming soon");
        }
    }
}
