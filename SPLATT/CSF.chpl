/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   CSF.chpl

    Description:    This is a module file for the compressed sparse fiber
                    (CSF) data structure used to represent sparse tensors.
*/

module CSF {
    use Base;
    use sptensor;
    use Args;

    /*****************************
    *
    *   CSF sparisty patter class
    *
    ******************************/
    class csf_sparsity {
        /** The size of each fptr and fids array */
        var nfibs : [NUM_MODES_d] int;
        /** The pointer structure for each sub-tree. fptr[f]
            marks the start of the children of node 'f'. This
            structure is a generalization of the 'rowptr' array
            used in CSR matrices */
        var fptr : [NUM_MODES_d] int;
        /** The index of each node. These map nodes back to the
            original tensor nonzeros */
        var fids : [NUM_MODES_d] int;
        /** The actual nonzero values. This array is of length
            nfibs[nmodes-1]. Since we do not know this size at
            compile time, we will resize the domain of this
            array later */
        var vals : [CSF_SPARSITY_VALS_d] real;
    }

    /*****************************
    *
    *   CSF class
    *
    ******************************/
    class splatt_csf {
        /** Number of nonzeros */
        var nnz : int;
        /** Number of modes */
        var nmodes : int;
        /** Dimension of each mode */
        var dims : [NUM_MODES_d] int;
        /** This maps levels in the tensor to actual tensor modes.
            dim_perm[0] is the mode stored at the root level and so on.
            NOTE: do not use directly; use 'csf_depth_to_mode' instead */
        var dim_perm : [NUM_MODES_d] int;
        /** Inverse of dim_perm. This maps tensor modes to levels in the
            CSF. NOTE: do not use directly; use 'csf_mode_to_depth' instead */
        var dim_iperm : [NUM_MODES_d] int;
        /** Which tiling scheme this tensor is stored as */
        var which_tile : splatt_tile_type;
        /** How many tiles there are */
        var ntiles : int;
        /** How many modes of the tensor (i.e. CSF levels) are tiled.
            Counted from the leaf (bottom) mode */
        var ntiled_modes : int;
        /** For a dense tiling, how many tiles along each mode */
        var tile_dims : [NUM_MODES_d] int;
        /** Sparsity structures: one for each tile. Since we
            do not know how many tiles there are yet, we'll
            resize this array later */
        var pt : [NUM_TILES_d] csf_sparsity;
    }

    /*****************************
    *
    *   Enums and Constants
    *
    ******************************/

    // The types of mode ordering available
    enum csf_mode_type 
    {  
        CSF_SORTED_SMALLFIRST,  /** sort the modes in non-decreasing order */
        CSF_SORTED_BIGFIRST,    /** sort the modes in non-increasing order */
        CSF_INORDER_MINUSONE,   /** one mode is placed first, rest naturally ordered*/
        CSF_SORTED_MINUSONE,    /** one mode is placed first, rest sorted by size */
        CSF_MODE_CUSTOM         /** custom mode ordering. dim_perm must be set! */
    };

    // CSF tile types
    enum splatt_tile_type
    {
        SPLATT_NOTILE,
        SPLATT_DENSETILE,
        /* DEPRECATED - pending CSF implementations */
        SPLATT_SYNCTILE,
        SPLATT_COOPTILE
    }; 

    /*****************************
    *
    *   Public Functions
    *
    ******************************/
    proc csf_alloc(tt : sptensor_t, args : cpd_cmd_args) : splatt_csf
    {

    }                                                       
}
