// For a coforall loop. a new thread is created at each iteration through the loop.
// This is useful in cases in which each iteration has a substantial amount of work and the
// number of tasks should be equal to the number of iterations.

use Memory;
config const printLocaleInfo = true;  // permit testing to turn this off

if printLocaleInfo then
  for loc in Locales do
    on loc {
      writeln("locale #", here.id, "...");
      writeln("  ...is named: ", here.name);
      writeln("  ...has ", here.numPUs(), " processor cores");
      writeln("  ...has ", here.physicalMemory(unit=MemUnits.GB, retType=real), 
              " GB of memory");
    }
writeln();

// here.numCores returns the number of cores your processor has
//config const numTasks = here.numCores;
//coforall tid in 1..numTasks {
//    writeln("Hello world from task ", tid, "of ", numTasks);
//}
