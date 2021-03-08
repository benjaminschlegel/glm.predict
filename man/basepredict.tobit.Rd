\name{basepredict.tobit}
\alias{basepredict.tobit}

\title{
predicted value
}
\description{
The function calculates the predicted value with the confidence interval for a tobit model.
}
\usage{
\method{basepredict}{tobit}(model, values, sim.count=1000, conf.int=0.95, sigma=NULL, set.seed=NULL,
  type = c("any", "simulation", "bootstrap"), summary = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the model Object generated with tobit()
}
  \item{values}{
the values of the case as vector in the order how they appear in the summary(model) Estimate
}
  \item{sim.count}{
OPTIONAL numbers of simulations to be done by the function. default: 1000
}
  \item{conf.int}{
OPTIONAL the confidence interval used by the function. default: 0.95
}
  \item{sigma}{
OPTIONAL the variance-covariance matrix, can be changed when having for example robust or clustered vcov. default: vcov(model)
}
  \item{set.seed}{
OPTIONAL set a seed for the random number generator
}
  \item{type}{
OPTIONAL only simulation is supported for tobit()
}
  \item{summary}{
OPTIONAL if mean/quantiles should be return or all simulated values (default: TRUE)	
}
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a 3x3 matrix having in the first column the predicted value, in the second
  column the lower value of the confidence interval and in the third column the upper value
  of the confidence interval.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
library(AER)
model1 = tobit(Age ~ Height + Pulse, right = 65, data=MASS::survey)
summary(model1)
# Person with a height of 160 and a pulse of 80
basepredict(model1, values = c(1,160,80))
}

\keyword{ models }
