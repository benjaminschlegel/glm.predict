\name{predicts}
\alias{predicts}

\title{
predicted values and discrete change
}
\description{
The function calculates the predicted values and the difference of
a range of cases with the confidence interval. It can be used for any
glm, polr or multinom model.
}
\usage{
predicts(model, values, position = NULL, sim.count = 1000, conf.int = 0.95, 
  sigma = NULL, set.seed = NULL, doPar = TRUE)
}
\arguments{
  \item{model}{
the model-Object generated with glm(), glm.nb(), polr() or multinom()
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
		\item \strong{log(from-to,by)}: takes the log of all values from "from" to "to" with the distance "by" (for example: "160-180,5" --> 160,165,170,175,180)
	\item \strong{log(from-to)}: same as log(from-to,by) with by=1 (for example: "2-8" --> 2,3,4,5,6,7,8)
	\item \strong{log(value1,value2,value3,...)}: takes the  log of the given values (for example: "160,180" --> 160,180)
	\item \strong{log(value1)}: takes the log of the given value (for example: "5.34" --> 5.34)
	\item \strong{"F"}: takes all values of a factor/character
	\item \strong{"F(1,4,7)"}: takes the first, fourth and seventh level of a factor/character
	\item \strong{"F(2)"}: takes the second level of a factor/character
}
}

  \item{position}{
OPTIONAL which variable should be taken for the discrete change, the variable must have at least two values. default: only predicted probabilities
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
  \item{doPar}{
OPTIONAL if the code should run parallel if more than 2 cores are detected
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
predicts(model1, "150-190,10;F;median", position = 1, doPar = FALSE)
}

\keyword{ models }
