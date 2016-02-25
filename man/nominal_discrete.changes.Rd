\name{nominal_discrete.changes}
\alias{nominal_discrete.changes}

\title{
predicted values and discrete changes for multinom models
}
\description{
The function calcultes the predicted values and the difference of
a range of cases with the confidence interval. It can be used for a multinom
model. - NOTICE: THE FUNCTION DOES NOT WORK WITH : INTERAKTIONS, ONLY WITH *!
}
\usage{
nominal_discrete.changes(model, values, data, position=1, sim.count=1000, conf.int=0.95, sigma=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the multinom-Object generated with multinom() from package nnet
}
  \item{values}{
the values of cases as character in the order how they appear in the summary(model) Estimate
The values must be in the following way: "value1;value2;value3;...". Each one of the value can be obe of the following:
"mean": takes the mean of that variable (can only be used when the variable is numeric) - 
"median": takes the median of that variable (can only be used when the variable is numeric) - 
"mode": takes the mode of that variable - 
"Q4": takes the quartiles (0,0.25,0.5,0.75,1) of that variable (other number for other quantiles) -
"min": takes the mimimum of that variable - 
"max": takes the maximum of that variable - 
from-to,by: takes all values from "from" to "to" with the distance "by" (for example: "160-180,5" --> 160,165,170,175,180) - 
from-to: same as from-to,by with by=1 (for example: "2-8" --> 2,3,4,5,6,7,8) - 
value1,value2,value3,...: takes the given values (for example: "160,180" --> 160,180) - 
value1: takes the given value (for example: "5.34" --> 5.34) -
"F4": takes all values of a factor with four levels (put as number the number of levels the factor has) -
"F5(2)": takes the second level of a factor with five levels (put as first number the number of levels the factor has) -
}

  \item{data}{
data.frame used to estimate the model
}

  \item{position}{
OPTIONAL which variable should be taken for the discrete change, the variable must have at least two values. default: 1
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
  The function makes a simulation for the all combination of cases and compares them to each other.
}
\value{
  The output is a data.frame with the predicted values and discrete changes.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
data = MASS::survey
data$Clap = relevel(data$Clap,ref="Neither")
model1 = nnet::multinom(Clap ~ Sex + Height, data=data)
summary(model1)
# comparing a woman to a man
# with the values 160, 170, 180 and 190 for the height
nominal_discrete.changes(model1, "0-1;160-190,10", data)
}

\keyword{ models }
