/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/18/2017
    File:   timers.chpl

    Description:    This is a module file for the timers we'll
                    be using.
*/

module Timers {
    use Time;
    use IO.FormattedIO;
    /*
        Create a class that has a timer for each of
        the various "things" we want to time.
    */
    class splatt_timers {
        // Associative array: maps strings (timer names) to Timer objects
        var timer_d = {"TOTAL", "CPD", "IO", "MTTKRP", "INVERSE", "SORT",
                       "CPD FIT", "MAT MULT", "MAT A^TA", "MAT NORM"};
        var timers : [timer_d] Timer;        

        /* Initialize the timers */
        proc initialize_timers()
        {
            for t in this.timers {
                t.clear();
            }   
        }
    
        /* Print out timers */
        proc report_times()
        {
            writeln("");
            writeln("Timing information ---------------------------------------------");
            for timerName in this.timer_d {
                if this.timers[timerName].elapsed() > 0 {
                    writef("  %-20s%0.3drs\n", timerName, this.timers[timerName].elapsed());
                }
            }
        }
    }
}
