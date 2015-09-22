\name{discrete.changes}
\alias{discrete.changes}

\title{
predicted values and discrete change
}
\description{
The function calcultes the predicted values and the difference of
a range of cases with the confidence interval. It can be used for any
glm model. - NOTICE: THE FUNCTION DOES NOT WORK WITH INTERAKTION THERM!
}
\usage{
discrete.changes(model, values, sim.count=1000, conf.int=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the glm-Object generated with glm()
}
  \item{values}{
the values of cases as character in the order how they appear in the summary(model) Estimate
The values must be in the following way: "value1;value2;value3;...". Each one of the value can be obe of the following:
"mean": takes the mean of that variable (can only be used when the variable is numeric) - 
"median": takes the median of that variable (can only be used when the variable is numeric) - 
"mode": takes the mode of that variable - 
from-to,by: takes all values from "from" to "to" with the distance "by" (for example: "160-180,5" --> 160,165,170,175,180) - 
from-to: same as from-to,by with by=1 (for example: "2-8" --> 2,3,4,5,6,7,8) - 
value1,value2,value3,...: takes the given values (for example: "160,180" --> 160,180) - 
value1: takes the given value (for example: "5.34" --> 5.34)
}
  \item{sim.count}{
OPTIONAL numbers of simulations to be done by the function. default: 1000
}
  \item{conf.int}{
OPTIONAL the confidence interval used by the function. default: 0.95
}
}
\details{
  The function makes a simulation for the all combination of cases and compares them to each other.
}
\value{
  The output is a data.frame with the predicted values and discrete changes.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
model1 = glm(Sex ~ Height + Smoke + Pulse, data=MASS::survey, family=binomial(link=logit))
summary(model1)
# comparing person with hight 150 to 160, 160 to 170, 170 to 180, 180 to 190 
# with all combination of(non-)smokers and a median of pulse 
discrete.changes(model1, "150-190,10;F4;median")
}

\keyword{ models }
