\name{gen.data}
\alias{gen.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate Data for BivPPL
}
\description{
  Generate input data for \code{BivPPL}
}
\usage{
gen.data(N, beta1, beta2, theta,lambda01,lambda02, c = 10, Ctype = TRUE, const_cov = F, same_cov = F)
}
\arguments{
  \item{N}{
  Sample size for alternating recurrent events, or number of clusters for other bivariate clustered events
}
  \item{beta1}{
  Regression parameters for event type 1. 
}
  \item{beta2}{
  Regression parameters for event type 2. 
}
  \item{theta}{
  A vector of 3 entries in the variance-covariance matrix of the bivariate frailty vector. The first two are variance values for the first and second event types, and the third is the covariance. 
}
  \item{c}{
  Maximum censoring time.
}
  \item{Ctype}{
  logical; if TRUE, fixed censoring time at c; or if FALSE, randomly draw from a uniform(0,c).
}
  \item{const_cov}{
  logical; if TRUE, the covariates are not updated across the event times. 
}
  \item{same_cov}{
  logical; if TRUE, the covariates are shared between two event types.
}
}

\value{
a list object of data that recorded all the covariates, event times, indicators, and the true frailty values. 
\item{x}{a list of follow-up time for the first event type}
\item{y}{a list of follow-up time for the second event type}
\item{z1}{a list of covariats for event type 1}
\item{z2}{a list of covariats for event type 2}
\item{delta1}{a list of type 1 event status}
\item{delta1}{a list of type 2 event status}
\item{r}{true values of the bivariate frailties}
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

data <- gen.data(N,beta1,beta2,theta,lambda01,lambda02,c=cen,Ctype=centype)
ptm<-proc.time()
res <- BivPPL(data)
proc.time() - ptm
res$beta_hat
res$beta_ASE
res$D_hat

}

\keyword{ Bivariate frailty }
\keyword{ Clustered events }
\keyword{ Simulation }
