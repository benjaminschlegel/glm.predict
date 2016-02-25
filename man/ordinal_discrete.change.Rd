\name{ordinal_discrete.change}
\alias{ordinal_discrete.change}

\title{
predicted values and discrete change for polr
}
\description{
The function calcultes the predicted values and the difference of
two cases with the confidence interval. It can be used for a polr model.
}
\usage{
ordinal_discrete.change(model, values1, values2, sim.count=1000, conf.int=0.95, sigma=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the polr-Object generated with polr()
}
  \item{values1}{
the values of case 1 as vector in the order how they appear in the summary(model) Estimate
}
  \item{values2}{
the values of case 2 as vector in the order how they appear in the summary(model) Estimate
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
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a nx9 matrix with the levels as rows and the different results as columns.
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
ordinal_discrete.change(model1, c(1,170),
  c(0,170),sim.count=10000)
}

\keyword{ models }
