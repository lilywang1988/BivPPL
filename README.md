Bivariate Survival Outcome Analysis Using Penalized Partial Likelihood (BivPPL)
================
[Lili Wang](mailto:lilywang@umich.edu)
2019-06-14

Purpose
-------

This R package is to analyze bivariate survival outcomes. Here specifically, we implement this method to analyze alternating recurrent events using a bivariate correlated frailty model. Both the regression parameters and the variance-covariance matrix will be estimated. This R package allow data to have many clusters. The estimation procedure is similar to the traditional PPL method as established in [coxme](https://cran.r-project.org/web/packages/coxme/index.html), while we do not require knowing the correlation direction between the correlated two events.

This R package will be improved and upgraded in the near future.

Install the package
-------------------

To install the R package from Github, you will need to install another R package "devtools". Please uncomment the codes to install them.

``` r
# install.packages("devtools")
# library(devtools)
# install_github("lilywang1988/BivPPL")
library(BivPPL)
#> Loading required package: survival
#> Loading required package: Matrix
#> Loading required package: mvtnorm
```

Vignettes
---------

Generate alternating recurrent event data using `gen.data`, `beta1` and `beta2` are regression parameters for the two events, `theta` is defining the 3 different entries of the variance-covariance matrix of the two correlated frailties. `N` is the number of clusters or subjects in the framework of alternating recurrent events.

``` r
set.seed(100)
N<-100 # number of clusters
beta1 <- c(0.5,-0.3,0.5) # regression parameters for event type 1
beta2 <- c(0.8,-0.2,0.3) # regression parameters for event type 2
beta  <- c(beta1,beta2)
theta <- c(0.25,0.25,-0.125) # variance-covariance matrix for the bivariate frailty (denoted as D), it is a vector (D[1,1],D[2,2], D[1,2])
lambda01 <- 1
lambda02 <- 1
cen <- 10 # maximum censoring time
centype <- TRUE # fixed censoring time at cen

data <- gen.data(N,beta1,beta2,theta,lambda01,lambda02)
ptm<-proc.time()
res <- BivPPL(data)
proc.time() - ptm
#>    user  system elapsed 
#>  20.554   0.940  21.875
res$beta_hat
#> [1]  0.52126635 -0.26535133  0.45676406  0.74664981 -0.06705498  0.24962866
res$beta_ASE
#> [1] 0.05652019 0.05266173 0.05729070 0.06239355 0.05456852 0.05988814
res$D_hat
#>            r1_hat     r2_hat
#> r1_hat  0.3791955 -0.2053474
#> r2_hat -0.2053474  0.2104612


# fit the model assuming the independence between the two frailties
ptm<-proc.time()
res2 <- BivPPL(data,independence=T)
proc.time() - ptm
#>    user  system elapsed 
#>   9.451   0.425   9.961
res2$beta_hat
#> [1]  0.51974445 -0.26578661  0.44356750  0.73519612 -0.06245907  0.23199444
res2$beta_ASE
#> [1] 0.05700419 0.05297851 0.05777039 0.06270781 0.05549611 0.06043905
res2$D_hat
#>           [,1]      [,2]
#> [1,] 0.3742242 0.0000000
#> [2,] 0.0000000 0.1866389

# A likelihood ratio test
LRT<-2*(res$LogMargProb-res2$LogMargProb)
print(round(pchisq(abs(LRT),df=1,lower.tail = F),3))
#>       [,1]
#> [1,] 0.002


# sparsen the Hessian matrix
ptm<-proc.time()
res3 <- BivPPL(data,huge=TRUE)
proc.time() - ptm
#>    user  system elapsed 
#>  12.116   0.037  12.185
res3$beta_hat
#>         z11         z12         z13         z21         z22         z23 
#>  0.52042943 -0.26522250  0.45584999  0.74594886 -0.06716247  0.24938898
res3$beta_ASE
#> [1] 0.05649124 0.05267874 0.05726610 0.06227382 0.05454575 0.05985640
res3$D_hat
#>            r1_hat     r2_hat
#> r1_hat  0.3728382 -0.2018795
#> r2_hat -0.2018795  0.2062161
```
