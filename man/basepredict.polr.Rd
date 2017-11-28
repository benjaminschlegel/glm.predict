\name{basepredict.polr}
\alias{basepredict.polr}

\title{
predicted value
}
\description{
The function calculates the predicted value with the confidence interval. It can be used for any
polr model.
}
\usage{
\method{basepredict}{polr}(model, values, sim.count=1000, conf.int=0.95, sigma=NULL, set.seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the model Object generated with polr()
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
OPTIONAL the variance-covairance matrix, can be changed when having for exaple robust or clustered vcov. default: vcov(model)
}
  \item{set.seed}{
OPTIONAL set a seed for the random number generator
}
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a matrix have in the first colone the predicted value, in the second
  colone the lower value of the confidence interval and in the third colone the upper value
  of the confidence interval.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
## For examples see ?basepredict
}

\keyword{ models }
