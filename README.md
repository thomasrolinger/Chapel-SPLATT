This branch will be where I try different methods of making the Chapel code faster.

+ It seems that running with less threads (i.e. less than 20) has better performance
than using all 20. Using gprof, a lot of the time is spent doing some coforall stuff.
So maybe the overhead of task creation and what not takes over? For YELP, running with
4 "threads" seems to be the best.


Here are the optimizations I will be trying. I'll try to keep this list up-to-date and
describe the outcomes:

1.) Use const variables whenever possible

RESULT: Minimal improvement, like 0.1 seconds on YELP

2.) Don't use c-pointers in MTTKRP, when possible

RESULT: No difference when doing this for the various CSF arrays but AWFUL performance
when this is done for the factor matrices. I suspect this is because instead of doing
pointer arithmetic to "get" to a certain row, we are doing array slicing.

3.) When iterating over arrays, use its domain whenever possible

RESULT: For the MTTKRP, this does not help. I suspect this may be because I have the --fast flag on,
which turns on --no-checks.

4.) Instead of computing the 2norm and maxnorm by hand, use the built-in
    functions. This will involve computing the norms of each matrix column and
    then doing a reduction on those sums.

RESULT: This works for the 2norm but that is only called for the first iteration, so it doesn't
have a big effective on the performance of the decomposition. It doesn't seem that Chapel has a
max norm built in.

5.) Use BLAS for mat mul
RESULT: Seems to be the best choice, about 1.5x faster than using dot from LinearAlgebra

6.) Implement mat mul from SPLATT
RESULT: Bad performance, wrong results

7.) Use 1D arrays for all matrices. The Chapel documentation claims that this will be better
RESULT: No difference, sticking with 2d matrices
