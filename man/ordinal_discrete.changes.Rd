\name{ordinal_discrete.changes}
\alias{ordinal_discrete.changes}

\title{
predicted values and discrete changes for polr models
}
\description{
The function calculates the predicted values and the difference of
a range of cases with the confidence interval. It can be used for a polr
model. - NOTICE: THE FUNCTION DOES NOT WORK WITH : INTERAKTIONS, ONLY WITH *!
}
\usage{
ordinal_discrete.changes(model, values, position=1, 
  sim.count=1000, conf.int=0.95, sigma=NULL, set.seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the polr-Object generated with polr() from package MASS
}
  \item{values}{
The values of cases as character in the order how they appear in the summary(model) Estimate.
The values must be in the following way: \emph{"value1;value2;value3;..."}. Each one of the values can be one of the following:

\itemize{
	\item \strong{"mean"}: takes the mean of that variable (can only be used when the variable is numeric)
	\item \strong{"median"}: takes the median of that variable (can only be used when the variable is numeric)
	\item \strong{"mode"}: takes the mode of that variable
	\item \strong{"Q4"}: takes the quartiles (0,0.25,0.5,0.75,1) of that variable (other number for other quantiles)
	\item \strong{"min"}: takes the mimimum of that variable
	\item \strong{"max"}: takes the maximum of that variable
	\item \strong{from-to,by}: takes all values from "from" to "to" with the distance "by" (for example: "160-180,5" --> 160,165,170,175,180)
	\item \strong{from-to}: same as from-to,by with by=1 (for example: "2-8" --> 2,3,4,5,6,7,8)
	\item \strong{value1,value2,value3,...}: takes the given values (for example: "160,180" --> 160,180)
	\item \strong{value1}: takes the given value (for example: "5.34" --> 5.34)
	\item \strong{"F"}: takes all values of a factor/character
	\item \strong{"F(2)"}: takes the second level of a factor/character
}
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
  \item{set.seed}{
OPTIONAL set a seed for the random number generator
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
data$Smoke = ordered(MASS::survey$Smoke,levels=c("Never","Occas","Regul","Heavy"))
model1 = MASS::polr(Smoke ~ Sex + Height, data=data, Hess=TRUE)
summary(model1)
# comparing a woman to a man
# with the values 160, 170, 180 and 190 for the height
ordinal_discrete.changes(model1, "0-1;160-190,10")
}

\keyword{ models }
