\name{discrete.change}
\alias{discrete.change}

\title{
predicted values and discrete change
}
\description{
The function calculates the predicted values and the difference of
two cases with the confidence interval. It can be used for any
glm model.
}
\usage{
discrete.change(model, values1, values2, sim.count=1000, conf.int=0.95, 
  sigma=NULL, sed.seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the glm-Object generated with glm() or glm.nb()
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
  \item{set.seed}{
OPTIONAL set a seed for the random number generator
}
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a 3x3 matrix have in the first colone the predicted values, in the second
  colone the lower value of the confidence interval and in the third colone the upper value
  of the confidence interval.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
model1 = glm(Sex ~ Height + Smoke + Pulse, data=MASS::survey, family=binomial(link=logit))
summary(model1)
# comparing a person with the height 150cm to 151cm
discrete.change(model1, c(1,150,1,0,0,mean(MASS::survey$Pulse,na.rm=TRUE)),
  c(1,151,1,0,0,mean(MASS::survey$Pulse,na.rm=TRUE)))
# the higher person has a greater probability to be a man
# the difference is significant, because the confidence interval
# does not include the 0
}

\keyword{ models }
