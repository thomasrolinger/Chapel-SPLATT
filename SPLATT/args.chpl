/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/07/2017
    File:   args.chpl

    Description:    This is a module file for the arguments passed into SPLATT
*/

module Args {
    use Time;
    use Base;
    /*
        Below are the command-line parameters that can be overriden:

        iters:  Maximum number of iterations to use (default: 50)
        tol:    Minimum change for convergence (default: 1e-5)
        reg:    Regularization parameter (default: 0)
        rank:   Rank of decomposition to find (default: 10)
        csf:    How many CSF to use; options are one, two, all (default: two)
        tile:   Use tiling during SPLATT (default: no (0))
        seed:   Random seed (default: current time)
    */
    config var iters: int = DEFAULT_ITS;
    config var tol: real = DEFAULT_TOL;
    config var reg: real = DEFAULT_REG;
    config var rank: int = DEFAULT_NFACTORS;
    config var csf: string = DEFAULT_CSF;
    config var tile: int = DEFAULT_TILE;
    config var seed: int = getCurrentTime():int;
    config var threads: int = DEFAULT_NUM_THREADS;

    /*
        Utility methods to print help, usage, etc.
    */
    proc printErrorMsg(args: [] string)
    {
        writeln("Usage: ", args[0], " TENSOR-FILE.bin [OPTION...]");
        writeln("Try '", args[0], " --help' or '", args[0], " --usage' for more info");
    }

    proc printUsageMsg(args: [] string)
    {
        writeln("Usage: ", args[0], " [--iters=NITERS] [--tol=TOLERANCE] [--reg=REGULARIZATION]");
        writeln("             [--rank=RANK] [--csf=#CSF] [--tile=TILING] [--seed=SEED] [--threads=NTHREADS]");
    }

    proc printHelpMsg(args: [] string)
    {
        writeln("Usage: ", args[0], " TENSOR-FILE.bin [OPTION...]");
        writeln("SPLATT-CPD: Compute the CPD of a sparse tensor");
        writeln("");
        writeln("  --iters=NITERS       Maximum number of iterations to use (default: 50)");
        writeln("  --tol=TOLERANCE      Minimum change for converage (default: 1e-5)");
        writeln("  --reg=REGULARIZATION Regularization parameter (default: 0)");
        writeln("  --rank=RANK          Rank of decomposition to find (default: 10)");
        writeln("  --csf=#CSF           How many CSF to use {one, two, all} (default: two)");
        writeln("  --tile=TILING        Use tiling during SPLATT {0=no, 1=yes} (default: no)");
        writeln("  --seed=SEED          Random seed (default: current time)");
        writeln("  --threads=NTHREADS   Number of tasks (threads) to yse (default=num cores)");
        writeln("");
        writeln("Report bugs to Thomas Rolinger <tbrolin@cs.umd.edu>");
    }

    /*
        Here, we define a class that holds the arguments so we can
        easily pass them around to functions
    */
    class cpd_cmd_args {
        var numIters : int;
        var tolerance : real;
        var regularization : real;
        var decompRank : int;
        var numCSF : string;
        var tiling : int;
        var rndSeed : int; 
        var tensorFile : string;
        var numThreads : int;

        /*
            Parses the arguments (mostly config arguments) and populates the fields
        */
        proc parseArgs(args: [] string)
        {
            /*
                For a valid execution to happen, the first argument to the
                program needs to be the tensor file. Everything else is optional
                and will be set to default values if not specified.
            */
            // Check for not enough arguments or help/usage
            if args.size <= 1 {
                printErrorMsg(args);
                exit(-1);
            }
            for arg in args {
                if arg == "--help" || arg == "-h" {
                    printHelpMsg(args);
                    exit(-1);
                }
                else if arg == "--usage" {
                    printUsageMsg(args);
                    exit(-1);
                }
            }
            // Get the tensor file (first argument)
            tensorFile = args[1];
            // Set the rest of the fields, which are determined by config variables.
            numIters = iters;
            tolerance = tol;
            regularization = reg;
            decompRank = rank;
            // Resize LAMBDA_d domain
            LAMBDA_d = 0..decompRank-1;
            numCSF = csf;
            tiling = tile;
            rndSeed = seed;
            numThreads = threads;
            // Set global numThreads
            numThreads_g = threads;
        }
    }
}
