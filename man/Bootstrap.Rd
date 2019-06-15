\name{Bootstrap}
\alias{Bootstrap}
\title{
Boostrap for BivPPL
}
\description{
Bootstrap for the standard error of the variance-covariance matrix estimators.
}
\usage{
Bootstrap(data.complete, B = 50, size.limit = 200, raw = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.complete}{
The complete dataset for bootstrap
}
  \item{B}{
The number of Boostrap replicates.
}
  \item{size.limit}{
To control the computation time, the Boostrap sample size is controlled to be less or equal to this value. If size.limit is larger or equal to the number of subjects in data.complete, this is a regular bootstrap; but if it is less than the number of subjects in data.complete, this becomes a m-out-of-n Boostrap. 
}
  \item{raw}{
Indicator of whether the input data is the output from \code{gen.data}. 
}
  \item{\dots}{
All the other arguments for \code{BivPPL}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{beta_B_mean}{Mean of the estimates from the Bootstrap replicates}
\item{beta_B_SE}{Estimated standard errors for the regression parameters via Boostrap}
\item{D_B_mean}{mean of variance-covariance matrix estimates from the Bootstrap replicates}
\item{D_B_SE}{Estimated standard errors for the variance-covariance matrix via Boostrap}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
lambda01 <- 2
lambda02 <- 2
cen <- 10 # maximum censoring time
centype <- TRUE # fixed censoring time at cen

data <- gen.data(N,beta1,beta2,theta,lambda01,lambda02)
ptm<-proc.time()
res <- BivPPL(data)
proc.time() - ptm
res$beta_hat
res$beta_ASE
res$D_hat
Bootstrap(data)$D_B_SE

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
