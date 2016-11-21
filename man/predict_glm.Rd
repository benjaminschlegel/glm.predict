\name{glm.predict_}
\alias{glm.predict_}

\title{
predicted value
}
\description{
The function calculates the predicted value with the confidence interval. It can be used for any
glm model.
}
\usage{
glm.predict_(model, values, sim.count=1000, conf.int=0.95, sigma=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the glm-Object generated with glm() or glm.nb()
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
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a 1x3 matrix have in the first colone the predicted value, in the second
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
glm.predict_(model1, c(1,150,1,0,0,mean(MASS::survey$Pulse,na.rm=TRUE)))
}

\keyword{ models }
