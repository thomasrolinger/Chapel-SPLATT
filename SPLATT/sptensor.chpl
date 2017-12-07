/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   sptensor.chpl

    Description:    This is a module file that defines the spare tensor related
                    functions and structures for SPLATT.
*/


module sptensor {
    use Base;
    use splatt_IO;
    /*****************************
    *
    *   Structures and Enums
    *
    ******************************/

    // Types of tensor supported by SPLATT
    enum tt_type  
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
        var nmodes: int;                /** The number of modes in the tensor, denoted as 'm' */
        var nnz: int;                   /** The number of non-zeros in the tensor */
        var dims: [NUM_MODES_d] int;    /** An array that contains the dimension of each mode */
        var ind: [COORD_d] int;         /** An m x nnz matrix that contains the coordinates of each
                                            non-zero. The nth non-zero is accessed via ind[0][n], ind[1][n],
                                            ..., ind[m][n]. */
        var vals: [NNZ_d] real;         /** An array containing the values of each non-zero */
        var tiled: int;                 /** Specifies whether sptensor_t has been titled; used by ftensor_t */

        /*######################################
        #
        #   Class functions
        #
        ########################################*/

        /*########################################################################
        #   Descriptipn:    Reads tensor from file and creates sptensor_t class
        #
        #   Parameters:     ifname (string):    Name of file to read from
        #
        #   Return:         sptensor_t: Parsed tensor
        ##########################################################################*/
        proc tt_read(ifname: string)
        {
            tt_read_file(ifname, this);
        }

        /*########################################################################
        #   Descriptipn:    Calculates the tensor of the tensor and returns it
        #
        #   Parameters:     None (called on this instance)
        #
        #   Return:         density: real
        ##########################################################################*/
        proc tt_calcDensity() : real
        {
            var root : real = nnz**(1.0/nmodes:real);
            var density : real = 1.0;
            for m in 0..nmodes-1 {
                density *= root / dims[m]:real;
            }
            return density;
        }
    }
}
