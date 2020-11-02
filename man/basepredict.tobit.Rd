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
  type = c("any", "simulation", "bootstrap"))
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
OPTIONAL the variance-covairance matrix, can be changed when having for exaple robust or clustered vcov. default: vcov(model)
}
  \item{set.seed}{
OPTIONAL set a seed for the random number generator
}
  \item{type}{
OPTIONAL choose between simulation and boostrap, "any" chooses between those two
according to the number of cases (bootstrap if n < 1000)
}
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a 3x3 matrix haveing in the first colone the predicted value, in the second
  colone the lower value of the confidence interval and in the third colone the upper value
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
