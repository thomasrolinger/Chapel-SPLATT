use Base;
use splatt_IO;
use sptensor;

// Driver/tester for SPLATT

proc main()
{
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
