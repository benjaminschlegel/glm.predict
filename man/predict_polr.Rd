\name{polr.predict}
\alias{polr.predict}

\title{
predicted value for polr
}
\description{
The function calculates the predicted value with the confidence interval. It can be used for a polr model.
}
\usage{
polr.predict(model, values, 
  sim.count=1000, conf.int=0.95, sigma=NULL, sed.seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the polr-Object generated with polr()
}
  \item{values}{
the values of the case as vector in the order how they appear in the summary(model) estimate
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
  The output is a nx3 matrix with the levels as rows and the different results as columns.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
data = MASS::survey
data$Smoke = ordered(MASS::survey$Smoke,levels=c("Never","Occas","Regul","Heavy"))
model1 = MASS::polr(Smoke ~ Sex + Height, data=data, Hess=TRUE)
summary(model1)
# comparing a man to a woman
polr.predict(model1, c(1,170),sim.count=10000)
}

\keyword{ models }
