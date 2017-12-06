use Base;
use splatt_IO;
proc main()
{
    var header: bin_header;
    header = new bin_header();
    var fin = open_f("/home/tbrolin/YELP.bin", iomode.r);
    var startPos = 0;
    // startPos is updated in this function so we know where to read next
    read_binary_header(fin, header, startPos);
    writeln("Magic: ", header.magic);
    writeln("IDX width: ", header.idx_width);
    writeln("VAL width: ", header.val_width);

    // In C, this would just be a regular int and we would pass in into
    // fill_binary_idx by reference, so it can be treated like an array.
    // But I don't think we can do that in Chapel, so it's a single-element array.
    var temp : [0..0] idx_t;
    // Again, startPos is updated here
    fill_binary_idx(temp, 1, header, fin, startPos);
    var nmodes = temp(0);
    var dims : [0..nmodes-1] int;
    fill_binary_idx(dims, nmodes, header, fin, startPos);
    fill_binary_idx(temp, 1, header, fin, startPos);
    var numNZ = temp(0);
    writeln("NMODES: ", nmodes);
    writeln("Dims: ", dims);
    writeln("NNZ: ", numNZ);    
}
