/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   CSF.chpl

    Description:    This is a module file for the compressed sparse fiber
                    (CSF) data structure used to represent sparse tensors.
*/

module CSF {
    /*****************************
    *
    *   Enums and Constants
    *
    ******************************/

    // The types of mode ordering available
    enum csf_mode_type = {  CSF_SORTED_SMALLFIRST,  /** sort the modes in non-decreasing order */
                            CSF_SORTED_BIGFIRST,    /** sort the modes in non-increasing order */
                            CSF_INORDER_MINUSONE,   /** one mode is placed first, rest naturally ordered*/
                            CSF_SORTED_MINUSONE,    /** one mode is placed first, rest sorted by size */
                            CSF_MODE_CUSTOM };      /** custom mode ordering. dim_perm must be set! */

    /*****************************
    *
    *   Public Functions
    *
    ******************************/
    
}
