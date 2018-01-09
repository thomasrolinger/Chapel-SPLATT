This branch will be where I try different methods of making the Chapel code faster.

Here are the optimizations I will be trying. I'll try to keep this list up-to-date and
describe the outcomes:

1.) Use const variables whenever possible

2.) When iterating over arrays, use its domain whenever possible

3.) Instead of computing the 2norm and maxnorm by hand, use the built-in
    functions. This will involve computing the norms of each matrix column and
    then doing a reduction on those sums.
