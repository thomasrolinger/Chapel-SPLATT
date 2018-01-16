# SPLATT + Chapel

Port of SPLATT to Chapel, with an accompanying benchmark script.

Currently how to compile Chapel code:

```bash
$ make
```

Executable will be in `bin/`. Optionally, you can pass in `PROFILE=1` to the
above command to enable profiling with gprof. You can also pass in `GEN_C_CODE=1`
to generate the C code that the Chapel code is compiled into.

To run the benchmark, go into the `benchmark` directory and then ensure that
the executable you wish to benchmark is in the `exces` directory. To benchmark
Chapel, make sure the executable is named `splatt_CHAPEL` and to benchmark C,
make sure the executable is named `splatt_C`. It is up to you to build and compile
the C-version of SPLATT from its repository: `https://github.com/ShadenSmith/splatt`.
The script `run_benchmark.sh` takes in two arguments: the language to benchmark
(C or CHAPEL) and the number of threads to use for the benchmark. This script is
designed so you can srun it on multiple nodes for different thread counts. The
benchmark will run across a suite of data sets, performing 5 trials per data set.

More details will come later (maybe).
