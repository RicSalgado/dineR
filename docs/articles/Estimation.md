# Estimation

*This Vignette covers the estimation function in depth.*

The main component of *dineR* is the *estimation* function that enables
any user the ability to perform differential network estimation with
great ease. The need for such a function stems for the fact that in high
dimensional settings, the inverse of any matrix cannot be determined
analytical. One such avenue in which this is overcome is through the use
of estimation procedures.

*dineR* includes the ability to perform this estimation under a number
of different loss functions available within literature as well as the
option to fine tune each of the estimation parameters to improve the
speed and accuracy of the process given specific problem
characteristics.

The focus of *dineR* lies primarily on differential networks, and the
subsequent estimation thereof. However, the structure of data suitable
for differential network estimation may not be immediately apparent. As
such, this function provides users with an easy, scalable manner to
generate appropriate data as they familiarize themselves with the
estimation options available within *dineR*

## Getting Started

The *data_generator* function enables users to generate two multivariate
normal samples of varying sizes, whose number of features/dimensions are
equal. These two samples are representative of the two experimental
conditions under which one would typically determine a differential
network.

Due to the technicalities of differential network estimation, it is not
possible to validate the accuracy of estimates obtained through the use
of traditional statistical methodologies. *data_generator* aims to
address this by generating data whose covariance matrices, and precision
matrices are known a-priori and as such so too is the differential
network. This allows users to validate the estimation procedure for
their specific experimental criteria.

## Arguments

The most important arguments available within *data_generator* are as
follows:

``` r

# Number of observations in sample 1
n_X <- 100 

# Number of observations in sample 2
n_Y <- 150

# Number of features in each of the samples
p <- 50

# The form of the precision covariance matrices
case <- "sparse"

# The seed of the simulation process to ensure reproducibility of results
seed <- 123
```

## Data Generation

Having defined each of the key arguments above, the *data_generator*
function can now be called:

``` r

data <- data_generator(n_X = n_X, n_Y = n_Y, p = p, case = case, seed = seed)
```

## Outputs

Given the data has now been generated, the data can be “wrangled” for
use within the estimation function.

``` r

# Extract the first sample
X <- data$X
cat("The number of observations in the first sample is:", nrow(X))
#> The number of observations in the first sample is: 100
cat("The number of features/dimensions in the first sample is:", ncol(X))
#> The number of features/dimensions in the first sample is: 50

# Extract the second sample
Y <- data$Y
cat("The number of observations in the second sample is:", nrow(Y))
#> The number of observations in the second sample is: 150
cat("The number of features/dimensions in the second sample is:", ncol(Y))
#> The number of features/dimensions in the second sample is: 50
```

It is also possible to extract the sample covariance matrices for each
of the samples:

``` r

# Extract the first sample's covariance matrix
Sigma_X <- data$Sigma_X

# Extract the second sample's covariance matrix
Sigma_Y <- data$Sigma_Y
```

In addition to the above matrices, the precision matrices can also be
extracted:

``` r

# Extract the first sample's precision matrix
Omega_X <- data$Omega_X

# Extract the second sample's precision matrix
Omega_Y <- data$Omega_Y
```

Lastly, it is possible to extract the differential network for
comparison to the estimate:

``` r

# Extract the differential network
Delta <- data$Delta
```

Having covered each of the key arguments and outputs available within
*data_generator*, please kindly note that further details are available
within the function documentation as well as the listed reference.
