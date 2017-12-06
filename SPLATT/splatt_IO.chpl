/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/05/2017
    File:   splatt_io.chpl

    Description:    This is a module file for IO related to SPLATT
*/

module splatt_IO {
    use IO;
    use Base;
    use sptensor;
    
    /*****************************
    *
    *   Enums and Classes
    *
    ******************************/
    enum splatt_file_type 
    {
        SPLATT_FILE_BIN_COORD  /** A binary version of the coordinate format */
    };

    // Associative domain/array that maps file extensions (strings) to the enums above
    var file_extensions_d = {".bin"};
    var file_extensions : [file_extensions_d] int;
    file_extensions[".bin"] = splatt_file_type.SPLATT_FILE_BIN_COORD;

    enum splatt_magic_type
    {
        SPLATT_BIN_COORD,
        SPLATT_BIN_CSF
    };

    /*
        This class is written to the beginning of any binary tensor file written by SPLATT
        As we aren't going to bother implementing file writing, we only use this class when
        reading in a file
    */
    class bin_header {
        var magic: int(32);
        var idx_width: uint(64);
        var val_width: uint(64);
    }

    /**************************************************************************
    *
    *   PUBLIC FUNCTIONS
    *
    ***************************************************************************/

    /*########################################################################
    #   Descriptipn:    Reads the file and populates a sptensor_t from it.
    #
    #   Parameters:     fname (string):     Name of file to read
    #                   tt (sptensor_t):    The sptensor_t to populate.
    #
    #   Return:         None
    ########################################################################*/
    proc tt_read_file(fname: string, tt : sptensor_t)
    {
        var fin = p_open_f(fname, iomode.r);
        // Always assume we're reading binary files
        p_tt_read_binary_file(fin, tt);
        try {
            fin.close();
        }
        catch {
            writeln("ERROR: Failed to close file");
        }
    }


    /**************************************************************************
    *
    *   PRIVATE FUNCTIONS
    *
    ***************************************************************************/

    /*########################################################################
    #   Descriptipn:    Open a file. Returns file handle/pointer.
    #
    #   Parameters:     fname (string): Name of file to open
    #                   mode (iomode):  Mode to open the file as
    #
    #   Return:         file: The handle/pointer to the opened file
    ########################################################################*/
    private proc p_open_f(fname: string, mode: iomode): file
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

    /*########################################################################
    #   Descriptipn:    Reads a tensor in coordinate format from a binary
    #                   file.
    #
    #   Parameters:     fin (file):         The opened file to read from
    #                   tt (sptensor_t):    The sptensor_t to populate
    #
    #   Return:         None
    ##########################################################################*/
    private proc p_tt_read_binary_file(fin: file, tt: sptensor_t)
    {
        // Read in header
        var header: bin_header;
        header = new bin_header();
        var startPos = 0;
        // startPos is updated in this function so we know where to read next
        p_read_binary_header(fin, header, startPos);

        // Read in next elements: number of modes, dimensions of each mode and
        // number of non-zeros

        /*
            In C, this would just be a regular int and we would pass in into
            fill_binary_idx by reference, so it can be treated like an array.
            But I don't think we can do that in Chapel, so it's a single-element array.
        */
        var temp : [0..0] idx_t;
        // startPos is also updated by this function
        p_fill_binary_idx(temp, 1, header, fin, startPos);
        // Pull out the nmodes from temp and store it as a single idx_t
        var nmodes = temp(0);
        // dims is an array of ints where each element holds the size of a mode.
        var dims : [0..nmodes-1] int;
        p_fill_binary_idx(dims, nmodes, header, fin, startPos);
        // Do the same thing for nnz that we did for nmodes
        p_fill_binary_idx(temp, 1, header, fin, startPos);
        var nnz = temp(0);
        // Update the domains in the Base module to use the actual sizes
        NNZ_d = 0..nnz-1;
        NUM_MODES_d = 0..nmodes-1;
        COORD_d = {0..nmodes-1, 0..nnz-1};

        if nmodes > MAX_NMODES {
            writeln("ERROR: maximum number of modes supported is ", MAX_NMODES);
            try {
                fin.close();
            }
            catch {
                writeln("ERROR: Failed to close file");
            }
            exit(-1);
        }

        // Read in tensor data into tt
        tt.dims = dims;
        tt.nmodes = nmodes;
        tt.nnz = nnz;
        for m in 0..nmodes-1 {
            p_fill_binary_idx(tt.ind[m,0..nnz-1], nnz, header, fin, startPos);
        }
        p_fill_binary_val(tt.vals, nnz, header, fin, startPos);
    }


    /*########################################################################
    #   Descriptipn:    Populate a binary header from an opened
    #                   file.
    #
    #   Parameters:     fin (file):             The opened file to read from
    #                   header (bin_header):    Binary header to read
    #                                           into. This is created by the
    #                                           caller and then modified by
    #                                           this function.
    #                   startPos (idx_t):       The byte-offset in the file
    #                                           to start reading from. This
    #                                           value is updated by this
    #                                           function to reflect the
    #                                           current position in the file
    #                                           after reading the header.
    #
    #   Return:         None
    ##########################################################################*/
    private proc p_read_binary_header(fin: file, header: bin_header, ref startPos: idx_t)
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
            try {
                fin.close();
            }
            catch {
                writeln("ERROR: Failed to close file");
            }
            exit(-1);
        }
        // Check for correct widths
        if header.idx_width > SPLATT_IDX_TYPEWIDTH / 8 {
            writeln("ERROR: Input has ", header.idx_width*8, "-bit integers but SPLATT_IDX_TYPEWIDTH is ", SPLATT_IDX_TYPEWIDTH, "-bits");
            try {
                fin.close();
            }
            catch {
                writeln("ERROR: Failed to close file");
            }
            exit(-1);
        }
        if header.val_width > SPLATT_VAL_TYPEWIDTH / 8 {
            writeln("ERROR: Input has ", header.val_width*8, "-bit floating point values but SPLATT_VAL_TYPEWIDTH is ", SPLATT_VAL_TYPEWIDTH, "-bits");
            try {
                fin.close();
            }
            catch {
                writeln("ERROR: Failed to close file");
            }
            exit(-1);
        }
    }

    /*########################################################################
    #   Descriptipn:    Reads count idx_t-bit elements from the file associated 
    #                   with fin and stores the elements into buffer. Reading
    #                   starts at the byte-offset represented by startPos.
    #                   This function modifies buffer and startPos. The 
    #                   header is used to determine the actual size of each 
    #                   element in the file to be read.
    #
    #   Parameters:     buffer (idx_t array):   Array of idx_t to store elements
    #                   count (idx_t):          Number of elements to read
    #                   header (bin_header):    Binary header that has been
    #                                           populated previously.
    #                   fin (file):             The opened file to read from
    #                   startPos (idx_t):       The byte-offset in the file
    #                                           to start reading from. This
    #                                           value is updated by this
    #                                           function to reflect the
    #                                           current position in the file
    #                                           after reading the elements.
    #
    #   Return:         None
    ##########################################################################*/
    private proc p_fill_binary_idx(buffer, count : idx_t, header: bin_header, fin: file, ref startPos: idx_t)
    {
        try {
            var idxSize : idx_t;
            if header.idx_width == SPLATT_IDX_TYPEWIDTH/8 {
                idxSize = SPLATT_IDX_TYPEWIDTH/8;
            }
            else {
                idxSize = 4;
            }
            // Read count uint32-bit elements in chunks. Ideally, we'd like to read in all the elements at once,
            // but buffer is an array of idx_t-bit elements, which may not necessarily be 32-bits. However, we
            // can parallelize this.
            var initialStart = startPos;
            const BUF_LEN = 1024*1024;
            forall i in 0..count-1 by BUF_LEN {
                // How many elements to read in
                var read_count = min(BUF_LEN, count-i);
                // Temporary array of uint32's to read into
                var temp : [0..read_count-1] uint(32);
                // Starting and ending pos for this chunk
                var st = (i*idxSize*read_count) + initialStart;
                var endPos = st + read_count*idxSize;
                var r = fin.reader(kind=ionative, locking=false, start=st, end=endPos);
                r.read(temp);
                // Copy over the chunk to the correct position in buffer
                for n in 0..read_count-1 {
                    buffer(n+i) = temp(n);
                }
                r.close();
            }   
            // update startPos
            startPos = (count*idxSize) + startPos;
        }
        catch {
            writeln("ERROR: Failed to read indices from file");
            try {
                fin.close();
            }
            catch {
                writeln("ERROR: Failed to close file");
            }
            exit(-1);
        }
    }

    /*########################################################################
    #   Descriptipn:    Same as fill_binary_idx but for values (floating point) 
    #
    #   Parameters:     buffer (val_t array):   Array of val_t to store elements
    #                   count (idx_t):          Number of elements to read
    #                   header (bin_header):    Binary header that has been
    #                                           populated previously.
    #                   fin (file):             The opened file to read from
    #                   startPos (idx_t):       The byte-offset in the file
    #                                           to start reading from. This
    #                                           value is updated by this
    #                                           function to reflect the
    #                                           current position in the file
    #                                           after reading the elements.
    #
    #   Return:         None
    ##########################################################################*/
    private proc p_fill_binary_val(buffer, count : idx_t, header: bin_header, fin: file, ref startPos: idx_t)
    {
        try {
            var valSize : idx_t;
            if header.val_width == SPLATT_VAL_TYPEWIDTH/8 {
                valSize = SPLATT_VAL_TYPEWIDTH/8;
            }
            else {
                valSize = 4;
            }
            var initialStart = startPos;
            const BUF_LEN = 1024*1024;
            forall i in 0..count-1 by BUF_LEN {
                // How many elements to read in
                var read_count = min(BUF_LEN, count-i);
                // Temporary array of real32's to read into
                var temp : [0..read_count-1] real(32);
                // Starting and ending pos for this chunk
                var st = (i*valSize*read_count) + initialStart;
                var endPos = st + read_count*valSize;
                var r = fin.reader(kind=ionative, locking=false, start=st, end=endPos);
                r.read(temp);
                // Copy over the chunk to the correct position in buffer
                for n in 0..read_count-1 {
                    buffer(n+i) = temp(n);
                }
                r.close();
            }     
            // update startPos
            startPos = (count*valSize) + startPos;
        }
        catch {
            writeln("ERROR: Failed to read indices from file");
            try {
                fin.close();
            }
            catch {
                writeln("ERROR: Failed to close file");
            }
            exit(-1);
        }
    }
}



