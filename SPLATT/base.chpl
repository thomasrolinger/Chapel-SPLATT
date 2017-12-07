/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   base.chpl

    Description:    This is a module file that has various constants and other
                    information used throughout SPLATT
*/

module Base {
    /*****************************
    *
    *   Defaults
    *
    ******************************/
    const MAX_NMODES: int = 8;
    const DEFAULT_TOL: real = 1e-5;
    const DEFAULT_NFACTORS: int  = 10;
    const DEFAULT_ITS: int = 50;
    const DEFAULT_CSF: string = "two";
    const DEFAULT_REG: int = 0;
    const DEFAULT_WRITE: int = 1;
    const DEFAULT_TILE: int = 0;
    const DEFAULT_NNZ: int = 1000;  /** Default number of non-zeros. Used when creating
                                        the domains below. */
    const SPLATT_IDX_TYPEWIDTH: int = 64;
    const SPLATT_VAL_TYPEWIDTH: int = 64;

    /*****************************
    *
    *   Commonly used domains
    *
    ******************************/
    var NNZ_d : domain(1) = 0..DEFAULT_NNZ-1;                   /** Domain for NNZ indices. Gets resized once NNZ is known */
    var NUM_MODES_d: domain(1) = 0..MAX_NMODES-1;               /** Possibly resize this once the
                                                                  number of modes is known */
    var COORD_d: domain(2) = {0..MAX_NMODES-1, 0..DEFAULT_NNZ-1}; /** Domain for 2D array with m rows and
                                                                  NNZ columns. Gets resized once m and NNZ is known */
    var CSF_SPARSITY_VALS_d : domain(1) = 0..DEFAULT_NNZ-1;     /** The csf_sparsity class has a 'vals' array whose length
                                                                    is not known until runtime. We'll resize this then */
    var NUM_TILES_d : domain(1) = 0..MAX_NMODES-1;              /** Number of tiles in CSF; not known until runtime so we'll
                                                                    resize this later. */
}
