## R CMD check results
There were no ERRORs or WARNINGs. There were 2 NOTEs:

* checking for future file timestamps: "unable to verify current time" — this is
  a transient network issue on the check machine and is not related to the
  package.

* Uses the superseded package 'doSNOW' — doSNOW is used intentionally as it is
  the only foreach parallel backend that supports passing a progress callback to
  workers via .options.snow. The alternative (doParallel) does not provide an
  equivalent mechanism, and switching would require removing progress reporting
  from the parallel execution path entirely.

## Updates
I have updated the package to include the ability to perform the optimization through the use of parallel computing.

I have also added in several of the core references whose theory was consulted when developing the package.

Lastly, I have updated the vignettes to be more comprehensive and cover the usage of the parallel parameter to perform parallel optimization.

