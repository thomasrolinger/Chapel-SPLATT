use Base;
use splatt_IO;
proc main()
{
    var header: bin_header;
    header = new bin_header();
    var fin = open_f("YELP.bin", iomode.r);
    read_binary_header(fin, header);
    writeln(header.magic);
    writeln(header.idx_width);
    writeln(header.val_width);
}
