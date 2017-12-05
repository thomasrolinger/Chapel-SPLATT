// For a coforall loop. a new thread is created at each iteration through the loop.
// This is useful in cases in which each iteration has a substantial amount of work and the
// number of tasks should be equal to the number of iterations.

// here.numPUs returns the number of cores your processor has
config const numTasks = here.numPUs();
coforall tid in 1..numTasks {
    writeln("Hello world from task ", tid, " of ", numTasks);
}
