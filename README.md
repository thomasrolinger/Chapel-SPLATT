# SPLATT + Chapel

Port of SPLATT to Chapel, with an accompanying benchmark script.

Currently, this is only the shared-memory implementation of SPLATT.
A multi-locale version will be started soon as a new branch.

Currently how to compile Chapel code:

```bash
$ make
```

This code uses BLAS and LAPACK. I have been only using OpenBLAS, which provides
both BLAS and LAPACK routines, so the build process is assuming that you are using OpenBLAS.
You'll need to modify the Makefile if you want to use a different library. You need to have OpenBLAS 
built on your system and then set the following environment variables:

```bash
export CBLAS_DIR=/path/to/blas/include/
export BLAS_LIBS=/path/to/blas/lib/
```

The executable will be in `bin/`. Optionally, you can pass in `PROFILE=1` to the
make command to enable profiling with gprof. You can also pass in `GEN_C_CODE=1`
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

Input tensors must be in SPLATT's binary format. So you'll need to build SPLATT (the C version)
and use its conversion tool to produce input tensors.

More details will come later (maybe).

# Things to work on

1.) Fix/look at CSF allocation and see why it seems to be "broken" for some tensors.
Possibly just reimplement it using SPLATT's newer parallel approach.
2.) Start on multi-locale version. This will require investigating the reference
C/MPI code and understanding what is going on.
3.) Do a pass through of the code and try to make things more "Chapel-like".
4.) Implementing mode tiling (perhaps only if it is necessary for the multi-locale
version; otherwise, it is something I rarely use myself in practice).
