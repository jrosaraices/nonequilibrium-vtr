Enclosed are the 'core' (i.e., _application-independent_) modules in the simulation code.
In detail:

- `kinds_m.F90` -- default integer and real 'kinds'; provides basic input/output and error-handling utilities
- `param_m.F90` -- linked list for parameter parsing and assignment
- `iobuf_m.F90` -- buffer for CSV-formatted output
- `prng_m.F90`  -- frontend for [Sebastiano Vigna's implementations](https://prng.di.unimi.it/) of the [_xorshift128+_](https://prng.di.unimi.it/xoshiro128plus.c) and [_xorshift1024*_](https://prng.di.unimi.it/xoroshiro1024star.c) shift-register [PRNGs for uniform random variates](https://en.wikipedia.org/wiki/Xorshift), providing seed generation procedures and interfaces for sampling from other basic distributions
- `prng_xorshift_sm.F90` -- Fortran bindings to Sebastiano Vigna's xorshift PRNGs implementations in `xorshift.c`
- `prng_uniform_sm.F90` -- procedures to sample from uniform distrbutions on integers and reals
- `prng_normal_sm.F90` -- procedures to sample from the univariate standard normal distribution
