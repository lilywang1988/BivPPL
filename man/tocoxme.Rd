\name{tocoxme}
\alias{tocoxme}
\alias{tocoxme_n}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Transform data to fit coxme for comparison
}
\description{
  Transform the output dataset from \code{gen.data} to fit in \code{coxme} for comparison.
}
\usage{
tocoxme(data)
tocoxme_n(data)
}
\details{
  \code{tocoxme} is designed for positively correlated two events, and \code{tocoxme_n} is designed for negatively correlated two events.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
  a list object from \code{gen.data}
}
}

\value{
\item{time}{follow-up times}
\item{delta}{event status}
\item{joint}{indicator for the two event types}
\item{b0,b1,b2}{indicators of the three frailties, b0 is the shared one which accounts for the correlations, b1 and b2 are for the first and second event types}
\item{Z1,Z2}{covariates for the two event types}
\item{Z0}{Only in tocoxme_n, }
\item{mi}{number of events within each cluster}
\item{N}{number of clusters}

}
\references{
%% ~put references to the literature/web site here ~
}
\examples{
# install coxme by:
#   install.packages("coxme")

library(coxme) 

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
res <- BivPPL(data,check.log.lik=TRUE) # check the convergence of the log-likelihood
proc.time() - ptm
res$beta_hat
res$beta_ASE
res$D_hat
data_coxme <- tocoxme(data) # assume that we do not know the two events are negatively associated and we transform the data to be positively associated
ptm<-proc.time()
res_coxme <- coxme(Surv(data_coxme$time,data_coxme$delta)~data_coxme$Z1+data_coxme$Z2+(1|data_coxme$b0)+(1|data_coxme$b1)+(1|data_coxme$b2)+strata(data_coxme$joint),eps=1e-6,list(reltol = 1e-5),sparse=c(1001,0.001),refine.n=0) # disable the Hessian matrix sparsening and likelihood refining
proc.time() - ptm
res_coxme$coefficients 
sqrt(diag(vcov(res_coxme)))
assemble(as.vector(unlist(res_coxme$vcoef)))

data_coxme_n <- tocoxme_n(data) # assume that we know the two events are negatively associated 
ptm<-proc.time()
res_coxme_n<- coxme(Surv(data_coxme_n$time,data_coxme_n$delta)~data_coxme_n$Z1+data_coxme_n$Z2+(data_coxme_n$Z0|data_coxme_n$b0)+(1|data_coxme_n$b1)+(1|data_coxme_n$b2)+strata(data_coxme_n$joint),eps=1e-6,list(reltol = 1e-5),sparse=c(1001,0.001),refine.n=0)
proc.time() - ptm
res_coxme_n$coefficients 
sqrt(diag(vcov(res_coxme_n)))
assemble_n(as.vector(unlist(res_coxme_n$vcoef)))
}

\seealso{
\code{\link[coxme:coxme]{coxme}} 
}

\keyword{Bivariate frailty}
\keyword{Clustered events}
