This branch will be where I try different methods of making the Chapel code faster.

Here are the optimizations I will be trying. I'll try to keep this list up-to-date and
describe the outcomes:

1.) Use const variables whenever possible
RESULT: Minimal improvement, like 0.1 seconds on YELP

2.) Don't use c-pointers in MTTKRP, when possible
RESULT: No difference when doing this for the various CSF arrays but AWFUL performance
when this is done for the factor matrices. I suspect this is because instead of doing
pointer arithmetic to "get" to a certain row, we are doing array slicing.

3.) When iterating over arrays, use its domain whenever possible

4.) Instead of computing the 2norm and maxnorm by hand, use the built-in
    functions. This will involve computing the norms of each matrix column and
    then doing a reduction on those sums.
