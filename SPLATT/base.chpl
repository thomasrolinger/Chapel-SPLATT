/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   CSF.chpl

    Description:    This is a module file that has various constants and other
                    information used throughout SPLATT
*/

module Base {
    /*****************************
    *
    *   Defines
    *
    ******************************/
    const MAX_NMODES: int = 8;
    type idx_t = int(64);       /** 64-bit ints for indices */
    type val_t = real(64);      /** doubles */

    /*****************************
    *
    *   Defaults
    *
    ******************************/
    const DEFAULT_TOL: real(64) = 1e-5;
    const DEFAULT_NFACTORS: idx_t = 10;
    const DEFAULT_ITS: idx_t = 50;
    const DEFAULT_WRITE: int = 1;
    const DEFAULT_TIEL: int = 0;
    const DEFAULT_NNZ: int = 1000;  /** Default number of non-zeros. Used when creating
                                        the domains below. */    

    /*****************************
    *
    *   Commonly used domains
    *
    ******************************/
    var NNZ_d : domain(1) = 0..DEFAULT_NNZ;                   /** Domain for NNZ indices. Gets resized once NNZ is known */
    var NUM_MODES_d: domain(1) = 0..MAX_NMODES;               /** Possibly resize this once the
                                                                  number of modes is known */
    var COORD_d: domain(2) = {0..MAX_NMODES, 0..DEFAULT_NNZ}; /** Domain for 2D array with m rows and
                                                                  NNZ columns. Gets resized once m and NNZ is known */
}
