/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   splatt_io.chpl

    Description:    This is a module file for IO related to SPLATT
*/

module splatt_IO {
    use IO;
    use Base;
    // Open a file. Returns file handle/pointer
    proc open_f(fname: string, mode: iomode): file
    {
        var f: file;
        try {
            f = open(fname, mode);
        }
        catch {
            writeln("***ERROR: Cannot open ", fname);
            exit(-1);
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

    // Populate a binary header from an input file. startPos refers to the byte offset
    // in the file to start reading from
    proc read_binary_header(fin: file, header: bin_header, ref startPos: idx_t)
    {
        try {
            // Create reader to read magic, idx_width and val_width
            var r = fin.reader(kind=ionative, locking=false, start=startPos, end=20);
            r.read(header.magic, header.idx_width, header.val_width);
            r.close();
            // updates startPos, since we passed in by reference
            startPos += 20;
        }
        catch {
            writeln("ERROR: Could not read header from file");
            exit(-1);
        }
        // Check for correct widths
        if header.idx_width > SPLATT_IDX_TYPEWIDTH / 8 {
            writeln("ERROR: Input has ", header.idx_width*8, "-bit integers but SPLATT_IDX_TYPEWIDTH is ", SPLATT_IDX_TYPEWIDTH, "-bits");
            exit(-1);
        }
        if header.val_width > SPLATT_VAL_TYPEWIDTH / 8 {
            writeln("ERROR: Input has ", header.val_width*8, "-bit floating point values but SPLATT_VAL_TYPEWIDTH is ", SPLATT_VAL_TYPEWIDTH, "-bits");
            exit(-1);
        }
    }

    // Fill an array of idx_t with values from a binary file. 'header' tells us whether
    // we can just read() the whole array or must read one at a time.
    // If the tensor was large enough that we needed 64-bit ints to represent its indices,
    // then we can read in the entire thing (we hard-code our indices as 64-bit for now).
    // Otherwise, we need to read in 32-bit values.
    proc fill_binary_idx(buffer, count : idx_t, header: bin_header, fin: file, ref startPos: idx_t)
    {
        try {
            var idxSize : idx_t;
            if header.idx_width == SPLATT_IDX_TYPEWIDTH/8 {
                idxSize = SPLATT_IDX_TYPEWIDTH/8;
            }
            else {
                idxSize = 4;
            }
            // Start the reader where we left off from read_binary_header
            var endPos = startPos;
            for i in 0..count-1 {
                startPos = endPos;
                endPos = startPos + idxSize;
                var r = fin.reader(kind=ionative, locking=false, start=startPos, end=endPos);
                var temp : uint(32);
                r.read(temp);
                buffer(i) = temp;
                r.close();
            }   
            // update startPos
            startPos = endPos;
        }
        catch {
            writeln("ERROR: Failed to read indices from file");
            exit(-1);
        }
    }
}
