# SPLATT + Chapel

Port of SPLATT to Chapel, with an accompanying benchmark script.

Currently, this is only the shared-memory implementation of SPLATT.
A multi-locale version will be started soon as a new branch.

Currently how to compile Chapel code:

```bash
$ make
```

This code uses BLAS, so you need to have a BLAS library built on your system
and then set the following environment variables:

```bash
export CBLAS_DIR=/path/to/blas/include/
export BLAS_LIBS=/path/to/blas/lib/
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
(C or CHAPEL) and the data set to use for the benchmark. Right now, the path to 
the actual data set is harcoded for my own system, so if you are actually going to
run this, then you need to edit the script.

More details will come later (maybe).

# Things to work on

1.) Fix/look at CSF allocation and see why it seems to be "broken" for some tensors.
Possibly just reimplement it using SPLATT's newer parallel approach.
2.) Use inline functions where ever SPLATT does
3.) Determine why matrix multiply is slower OR get the more numerically stable
approach for computing the inverse to work.
4.) Start on multi-locale version. This will require investigating the reference
C/MPI code and understanding what is going on.
5.) Implementing mode tiling (perhaps only if it is necessary for the multi-locale
version; otherwise, it is something I rarely use myself in practice).
