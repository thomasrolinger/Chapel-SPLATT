/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   splatt_io.chpl

    Description:    This is a module file for IO related to SPLATT
*/

module splatt_IO {
    use IO;
    // Open a file. Returns file handle/pointer
    proc open_f(fname: string, mode: iomode): file
    {
        var f: file;
        try {
            f = open(fname, mode);
        }
        catch {
            writeln("***ERROR: Cannot open ", fname);
        }
        return f;
    }

    /*****************************
    *
    *   File Types
    *
    ******************************/
    enum splatt_file_type 
    {
        SPLATT_FILE_TEXT_COORD, /** Plain list of tuples + values */
        SPLATT_FILE_BIN_COORD  /** A binary version of the coordinate format */
    };

    // Associative domain/array that maps file extensions (strings) to the enums above
    var file_extensions_d = {".tns", ".coo", ".bin"};
    var file_extensions : [file_extensions_d] int;
    file_extensions[".tns"] = splatt_file_type.SPLATT_FILE_TEXT_COORD;
    file_extensions[".coo"] = splatt_file_type.SPLATT_FILE_TEXT_COORD;
    file_extensions[".bin"] = splatt_file_type.SPLATT_FILE_BIN_COORD;

    // Given a file name, fname, determine its file type and return it
    proc get_file_type(fname: string): splatt_file_type
    {
        // Split fname on ".", and take the last token
        var tokens = fname.split(".");
        var fileExt = ".".join(["", tokens[tokens.size]]);
        // Check for matching extension in file_extensions
        if file_extensions_d.member(fileExt) then
            return file_extensions[fileExt];
        else
            writeln("*** ERROR: extension for ", fname, " not recognized. Defaulting to ASCII coordinate format");
            return splatt_file_type.SPLATT_FILE_TEXT_COORD;
        
    }
    
    /*****************************
    *
    *   Binary Reads
    *
    ******************************/
    enum splatt_magic_type
    {
        SPLATT_BIN_COORD,
        SPLATT_BIN_CSF
    };

    // This class is written to the beginning of any binary tensor file written by SPLATT
    class bin_header {
        var magic: int(32);
        var idx_width: uint(64);
        var val_width: uint(64);
    }

    // Populate a binary header from an input file
    proc read_binary_header(fin: file, header: bin_header)
    {
        try {
            // Create reader to read magic, idx_width and val_width
            var r = fin.reader(kind=ionative, locking=false, start=0, end=32+64+64);
            r.read(header.magic, header.idx_width, header.val_width);
            r.close();
        }
        catch {
            writeln("ERROR: Could not read header from file");
        }
    }
}
