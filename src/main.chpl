/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/07/2017
    File:   main.chpl

    Description:    This is the main/driver program.
*/

use Base;
use splatt_IO;
use sptensor;
use Args;
use Stats;
use CSF;
use Kruskal;
use CPD;
use Timers;
use Time;

proc main(args: [] string)
{
    // Initialize timers
    timers_g = new splatt_timers();
    timers_g.initialize_timers();

    // Start ALL timer
    timers_g.timers["TOTAL"].start();

    // Parse command line arguments
    var cpdArgs : cpd_cmd_args = new cpd_cmd_args();
    cpdArgs.parseArgs(args);

    // Seed rng
    randStream_g = new RandomStream(int(32), cpdArgs.rndSeed);
    
    // Print a header
    writeln("****************************************************************");
    writeln("SPLATT + Chapel v1.0");
    writeln("****************************************************************");
    writeln("");

    // Create sptensor and parse the input file
    var tt : sptensor_t;
    tt = new sptensor_t();
    tt.tt_read(cpdArgs.tensorFile);

    // Print basic tensor stats
    stats_tt(tt, cpdArgs.tensorFile);

    // csf is an array. Depending on the csf allocation
    // type we have, each element in the array will be a csf
    // for that correspond mode:
    // one : csf[0] is only CSF we care about
    // two : csf[0] and csf[1]
    // all : all elements in csf
    var csf = csf_alloc(tt, cpdArgs);

    var nmodes : int = tt.nmodes;

    // Now that the CSF is built, we no longer need the sptensor tt
    delete tt;
    
    // Print CPD stats
    cpd_stats(csf, cpdArgs);

    // Create Kruskal tensors (holds info/output of CPD).
    var factored : splatt_kruskal = new splatt_kruskal();

    // Do the factorization
    var ret : int = splatt_cpd_als(csf, cpdArgs, factored);
    if ret != SPLATT_SUCCESS {
        writeln("ERROR: splatt_cpd_als returned ", ret, ". Aborting");
    }

    // End ALL timer and report times
    timers_g.timers["TOTAL"].stop();
    timers_g.report_times();
    writeln("****************************************************************");
}
