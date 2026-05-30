# Parallelisation

*This vignette covers parallel vs sequential execution in the estimation
function.*

*dineR* supports parallel execution of the ADMM optimisation across the
regularisation path via the `cores` argument of the `estimation`
function. By default, the function runs sequentially (`cores = 1`).
Setting `cores` to a value greater than 1 distributes each lambda value
across the specified number of worker processes, which can yield
substantial reductions in wall-clock time for medium-to-large problems.

## Switching Between Sequential and Parallel Mode

The `cores` argument is the only change required to move between modes.
The value must be a positive integer no greater than
`parallel::detectCores() - 1`.

``` r

# Check how many cores are available on your machine
parallel::detectCores()
```

``` r

# Generate data for the examples below
data <- data_generator(n_X = 150, p = 100, seed = 42)
X <- data$X
Y <- data$Y
```

**Sequential** (default, `cores = 1`):

``` r

result_seq <- estimation(X, Y, nlambda = 15, cores = 1)
result_seq$elapse  # elapsed time in seconds
```

**Parallel** (`cores` set to the desired number of workers):

``` r

result_par <- estimation(X, Y, nlambda = 15, cores = 4)
result_par$elapse  # elapsed time in seconds
```

Both calls return identical results — `cores` affects only computation
time, not the estimates or any other output field.

## Performance Benchmarks

The table below reports wall-clock times and speed-up factors measured
on an Apple M-series machine (8 logical cores) using 4 parallel workers.
Each scenario uses the default LASSO loss and no tuning.

| Scenario | Observations (n) | Dimensions (p) | Lambda values (nlambda) | Sequential | Parallel (4 cores) | Speed-up |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|
| Small | 100 | 50 | 10 | 0.68 s | 0.70 s | 1.0× |
| Medium | 150 | 100 | 15 | 3.89 s | 1.76 s | 2.2× |
| Large | 200 | 150 | 20 | 22.4 s | 6.98 s | 3.2× |
| High-dim | 100 | 200 | 20 | 17.9 s | 5.71 s | 3.1× |
| Many-λ | 150 | 100 | 30 | 7.72 s | 2.77 s | 2.8× |

The code used to produce these results is provided below for
reproducibility:

``` r

library(dineR)

run_bench <- function(label, n_X, p, nlambda, cores_par) {
  data <- data_generator(n_X = n_X, p = p, seed = 42)
  X <- data$X; Y <- data$Y

  r_seq <- estimation(X, Y, nlambda = nlambda, cores = 1)
  r_par <- estimation(X, Y, nlambda = nlambda, cores = cores_par)

  cat(sprintf(
    "[%s] seq: %.3fs | par(%d cores): %.3fs | speed-up: %.2fx\n",
    label, r_seq$elapse, cores_par, r_par$elapse,
    r_seq$elapse / r_par$elapse
  ))
}

run_bench("Small",    n_X = 100, p = 50,  nlambda = 10, cores_par = 4)
run_bench("Medium",   n_X = 150, p = 100, nlambda = 15, cores_par = 4)
run_bench("Large",    n_X = 200, p = 150, nlambda = 20, cores_par = 4)
run_bench("High-dim", n_X = 100, p = 200, nlambda = 20, cores_par = 4)
run_bench("Many-lam", n_X = 150, p = 100, nlambda = 30, cores_par = 4)
```

## Interpreting the Results

Several patterns emerge from the benchmarks:

- **Small problems** (`p = 50`, `nlambda = 10`): parallelisation
  provides no benefit. The overhead of spawning worker processes and
  distributing work exceeds the per-lambda computation time. Sequential
  execution is preferred here.

- **Medium-to-large problems** (`p ≥ 100` or `nlambda ≥ 15`): consistent
  2–3× speed-ups are observed with 4 cores. Gains are driven by both the
  dimensionality (`p` increases the per-lambda solve time) and the
  number of lambda values (`nlambda` increases the number of independent
  tasks that can be distributed).

- **Speed-up scales with workload**: the largest absolute savings occur
  for high-dimensional or fine-grained regularisation paths, where each
  lambda solve is computationally expensive.

## Choosing the Number of Cores

A practical rule of thumb:

- Use `cores = 1` when `p < 100` and `nlambda < 15`.
- Use `cores = min(nlambda, parallel::detectCores() - 1)` for larger
  problems, reserving one core for the main R process.

``` r

# Recommended cores selection for larger problems
n_cores <- min(nlambda, parallel::detectCores() - 1)
result <- estimation(X, Y, nlambda = nlambda, cores = n_cores)
```

Note that `estimation` will automatically reduce `cores` to `nlambda` if
more cores are requested than there are lambda values, and will warn
accordingly.
