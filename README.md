# SPLATT + Chapel

Port of SPLATT to Chapel, with an accompanying benchmark script.

Currently how to compile Chapel code:

```bash
$ make
```

Executable will be in `bin/`.

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

# Optimizations that I have tried

##1.) Use const variables whenever possible

RESULT: Minimal improvement, like 0.1 seconds on YELP

##2.) Don't use c-pointers in MTTKRP, when possible

RESULT: No difference when doing this for the various CSF arrays but AWFUL performance
when this is done for the factor matrices. I suspect this is because instead of doing
pointer arithmetic to "get" to a certain row, we are doing array slicing.

##3.) When iterating over arrays, use its domain whenever possible

RESULT: For the MTTKRP, this does not help. I suspect this may be because I have the --fast flag on,
which turns on --no-checks.

##4.) Instead of computing the 2norm and maxnorm by hand, use the built-in
    functions. This will involve computing the norms of each matrix column and
    then doing a reduction on those sums.

RESULT: This works for the 2norm but that is only called for the first iteration, so it doesn't
have a big effective on the performance of the decomposition. It doesn't seem that Chapel has a
max norm built in.

##5.) Use BLAS for mat mul
RESULT: Seems to be faster than the dot routine but a hand-written tiled mat mul seems to be the
best option.

##6.) Use 1D arrays for all matrices. The Chapel documentation claims that this will be better
RESULT: No difference, sticking with 2d matrices

