# NPN - Non paranormal Transformation

This functions allows us to transform non-normal multivariate data to
that of non paranormal data.

## Usage

``` r
npn(x, npn_func = "shrinkage", npn_thresh = NULL, verbose = TRUE)
```

## Arguments

- x:

  The multivariate non-normal data to be transformed.

- npn_func:

  Optional parameter - The method of transformation to be applied. Can
  either be "shrinkage" or "truncation" but defaults to "shrinkage".

- npn_thresh:

  Optional parameter - The truncation threshold that is used when making
  use of truncation.

- verbose:

  Optional parameter - Prints additional output of the selected
  approach. Can either be "TRUE" or "FALSE" and defaults to "TRUE".

## Value

Returns the transformed data matrix.

## References

Liu, H., Han, F., Yuan, M., Lafferty, J. and Wasserman, L., 2012. The
nonparanormal skeptic. arXiv preprint arXiv:1206.6488.

Liu, H., Lafferty, J. and Wasserman, L., 2009. The nonparanormal:
Semiparametric estimation of high dimensional undirected graphs. Journal
of Machine Learning Research, 10(10).

Xue, L. and Zou, H., 2012. Regularized rank-based estimation of
high-dimensional nonparanormal graphical models. The Annals of
Statistics, 40(5), pp.2541-2571.

## Examples

``` r
data <- data_generator(n_X = 100, p = 50, seed = 123)
X <- data$X
X_transformed <- npn(X, npn_func = "truncation")
#> Nonparanomral transformation via truncated ECDF.
```
