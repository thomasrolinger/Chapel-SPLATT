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

proc main(args: [] string)
{
    // Parse command line arguments
    var cpdArgs : cpd_cmd_args = new cpd_cmd_args();
    cpdArgs.parseArgs(args);

    var tt : sptensor_t;
    tt = new sptensor_t();
    tt.tt_read("/home/tbrolin/YELP.bin");
    writeln("num modes = ", tt.nmodes);
    for i in 0..9 {
        write("The #", i+1, " nnz is at (");
        for m in 0..tt.nmodes-1 {
            write(tt.ind[m,i], ",");
        }
        writeln(")");
    }
    for i in 0..9 {
        writeln("The #", i+1, " value is: ", tt.vals[i]);
    }
}
