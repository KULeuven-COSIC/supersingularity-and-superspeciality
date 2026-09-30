## Supersingularity and Superspeciality Verification of Abelian Surfaces
Magma code accompanying the paper *"Supersingularity and Superspeciality Verification of Abelian Surfaces"* by Maria Corte-Real Santos, Gioella Lorenzon, and Krijn Reijnders.

### Core Functionality
The core algorithms described in sections 4 and 5 are in the `algorithms` folder.

- `probabilistic-algorithms.m` contains the probabilistic verification algorithms from Section 4, such as Algorithm 3, and Algorithm 4 (`IsSupersingular`)
- `conclusive-algorithms.m` contains the conclusive verification algorithms from Section 5, such as those in Corollary 5.5, Corollary 5.10, Lemma 5.13, and Lemma 5.15.

### Benchmarking Results
We run four benchmarks: the probabilistic tests for supersingularity and superspeciality,
and the conclusive tests for supersingularity for generic primes, and primes `p` with a large 2-adic 
factor in `p+1`. The benchmarking scripts are in the `benchmarking` folder.
The results of these benchmarks are in the `benchmark_data` folder, from which
we generate Figures 1 and 2.

When in `bechmarking`, the following scripts perform specific comparisons.
- `test_supersingular.m` runs benchmarking of probabilistic supersingularity tests, results in `benchmark_data_ssing_prob.dat`
- `test_superspecial.m` runs benchmarking of probabilistic superspeciality tests, results in `benchmark_data_sspecial_prob.dat`
- `test_conclusive_big_p.m` and `test_conclusive_big_p_count.m` run benchmarking of conclusive supersingularity tests, results in `benchmark_data_conc_big_p.dat` and `benchmark_data_conc_big_p_count.dat` respectively, combined together in `benchmark_data_conc_big_p_together.dat`
- `test_conclusive_crypto.m` runs benchmarking of conclusive supersingularity testing for primes of the form `p = c*2^k - 1`, results in `benchmark_data_conc_crypto.dat`

### Helper functions
The folder `helpers` contains several helper functions to do arithmetic on Jacobians, sample random Jacobians, and so on.

### LMFDB Data
The LMFDB [1] has been instrumental for use to compute data on supersingular genus-2 hyperelliptic curves

- `helpers/ss_curves_29_211.m` contains data on supersingular genus-2 hyperelliptic curves for primes `29 <= p <= 211`


### Citations
[1] The LMFDB Collaboration, The L-functions and modular forms database, https://www.lmfdb.org, 2026.
