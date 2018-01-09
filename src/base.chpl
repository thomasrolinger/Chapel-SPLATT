/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   base.chpl

    Description:    This is a module file that has various constants and other
                    information used throughout SPLATT
*/

module Base {
    use Random;
    use Timers;
    use MutexPool;

    /*****************************
    *
    *   Defaults
    *
    ******************************/
    param MAX_NMODES: int = 8;
    param DEFAULT_TOL: real = 1e-5;
    param DEFAULT_NFACTORS: int  = 10;
    param DEFAULT_ITS: int = 50;
    param DEFAULT_CSF: string = "two";
    param DEFAULT_REG: real = 0.0;
    param DEFAULT_WRITE: int = 1;
    param DEFAULT_TILE: int = 0;
    const DEFAULT_NUM_THREADS: int = here.numPUs();
    param MIN_QUICKSORT_SIZE : int = 8;
    param SMALL_SORT_SIZE : int = 1000;
    param DEFAULT_NNZ: int = 1000;  /** Default number of non-zeros. Used when creating
                                        the domains below. */
    param SPLATT_IDX_TYPEWIDTH: int = 64;
    param SPLATT_VAL_TYPEWIDTH: int = 64;
    param SPLATT_SUCCESS: int = 1;
    param RAND_MAX: int = 2147483647;
    param DEFAULT_PRIV_THRESH : real = 0.02; 
    param DEFAULT_NLOCKS : int = 1024;
    param DEFAULT_LOCK_PAD : int = 8; // this is 8*8 bytes = 64 bytes. Used to prevent false sharing
    param MAT_NORM_2 : int = 0;
    param MAT_NORM_MAX : int = 1;
    param REDUCE_SUM : int = 0;
    param REDUCE_MAX : int = 1;

    // Global "numThreads" so we don't have to pass it around all the time.
    // We set this once we parse the args
    var numThreads_g : int;

    // Global random stream to use. We will seed it once we parse command line args
    var randStream_g : RandomStream(int(32));

    // Timers
    var timers_g : splatt_timers;
    
    // Something in thread partition
    var nprobes_g : int = 0;

    // Mutex pool
    var pool_g : mutex_pool;

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
    var LAMBDA_d : domain(1) = 0..DEFAULT_NFACTORS-1;
    var FACTORS_d : domain(1) = 0..DEFAULT_NFACTORS-1; // commonly used domain
}
