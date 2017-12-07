/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/07/2017
    File:   util.chpl

    Description:    This is a module file that contains misc. utility
                    functions
*/

module Util {
    /* 
        Takes a value (number of bytes) and produces a printable
        string that represents the size (i.e. B, KB, MB, etc.)
    */
    proc bytes_str(bytes : int) : string
    {
        var size: real = bytes:real;
        var suff: int = 0;
        var suffix: [0..4] string;
        suffix = ["B", "KB", "MB", "GB", "TB"];
        while (size > 1024 && suff < 5) {
            size /= 1024;
            suff += 1;
        }
        var ret : string = size:string;
        ret += suffix[suff];
        return ret;
    }
    
}
