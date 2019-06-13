\name{BivPPL-package}
\alias{BivPPL-package}
\alias{BivPPL}
\docType{package}
\title{
\packageTitle{BivPPL}
}
\description{
\packageDescription{BivPPL}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BivPPL}
\packageIndices{BivPPL}
The R package
}
\author{
\packageAuthor{BivPPL}

Maintainer: \packageMaintainer{BivPPL}
}
\references{
%~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
%\code{\link[coxme]{coxme}} 
}
\examples{
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
res$beta_hat
res$beta_ASE
res$D_hat

# fit the model assuming the independence between the two frailties
ptm<-proc.time()
res2 <- BivPPL(data,independence=T)
proc.time() - ptm
res2$beta_hat
res2$beta_ASE
res2$D_hat

# A likelihood ratio test
LRT<-2*(res$LogMargProb-res2$LogMargProb)
print(round(pchisq(abs(LRT),df=1,lower.tail = F),3))


# sparsen the Hessian matrix
ptm<-proc.time()
res3 <- BivPPL(data,huge=TRUE)
proc.time() - ptm
res3$beta_hat
res3$beta_ASE
res3$D_hat
}
