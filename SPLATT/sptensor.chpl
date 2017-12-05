/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   sptensor.chpl

    Description:    This is a module file that defines the spare tensor related
                    functions and structures for SPLATT.
*/

use Base

module sptensor {
    /*****************************
    *
    *   Structures and Enums
    *
    ******************************/

    // Types of tensor supported by SPLATT
    enum tt_type =  
    {
        SPLATT_3MODE, /** Three-mode tensors. */
        SPLATT_NMODE  /** Tensors of with arbitrary numbers of modes.
                          NOTE: support is minimal. */
    };

    // The main data structure for representing sparse tensors in coordinate format
    class sptensor_t {
        /*######################################
        #
        #   Class/struct fields
        #
        ########################################*/
        
        var tensorType: tt_type;        /** Type of tensor represented */
        var nmodes: idx_t;              /** The number of modes in the tensor, denoted as 'm' */
        var nnz: idx_t;                 /** The number of non-zeros in the tensor */
        var dims: [NUM_MODES_d] idx_t;  /** An array that contains the dimension of each mode */
        var ind: [COORD_d] idx_t;       /** An m x nnz matrix that contains the coordinates of each
                                            non-zero. The nth non-zero is accessed via ind[0][n], ind[1][n],
                                            ..., ind[m][n]. */
        var vals: [NNZ_D] val_t;        /** An array containing the values of each non-zero */
        var tiled: int;                 /** Specifies whether sptensor_t has been titled; used by ftensor_t */

        /*######################################
        #
        #   Class functions
        #
        ########################################*/
        // Loads a sparse tensor from the file ifname.
        // Modifies the class/struct fields
        proc tt_read(ifname: string) 
        {

        }
    }
}
