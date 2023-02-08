\name{dc.vglm}
\alias{dc.vglm}

\title{
predicted values and discrete change
}
\description{
The function calculates the predicted values and the difference of
two cases with the confidence interval. It can be used for a vglm model.
}
\usage{
\method{dc}{vglm}(model, values = NULL, sim.count = 1000, conf.int = 0.95, 
  sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
  type = c("any", "simulation", "bootstrap"), summary = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the model-Object generated with vglm()
}
  \item{values}{
the values of case 1 and 2 as vector in the order how they appear in the summary(model) Estimate. Values is if values1 and values2 are specified after each other in the same vector. Either values or values1 and values2 have to be specified.
}
  \item{sim.count}{
OPTIONAL numbers of simulations to be done by the function. default: 1000
}
  \item{conf.int}{
OPTIONAL the confidence interval used by the function. default: 0.95
}
  \item{sigma}{
OPTIONAL the variance-covariance matrix, can be changed when having for example robust or clustered vcov. default: vcov(model)
}
  \item{set.seed}{
OPTIONAL set a seed for the random number generator
}
  \item{values1}{
the values of case 1 as vector in the order how they appear in the summary(model) Estimate. Has to be defined if values is not defined.
}
  \item{values2}{
the values of case 2 as vector in the order how they appear in the summary(model) Estimate. Has to be defined if values is not defined.
}
  \item{type}{
OPTIONAL choose between simulation and bootstrap, "any" chooses between those two
according to the number of cases (bootstrap if n < 500) (boostrap is very slow for vglm models)
}
  \item{summary}{
OPTIONAL if mean/quantiles should be return or all simulated values (default: TRUE)	
}
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a matrix have in the first column the predicted values, in the second
  column the lower value of the confidence interval and in the third column the upper value
  of the confidence interval.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
\dontrun{
data = MASS::survey
data$Smoke = ordered(data$Smoke, levels = c("Never", "Occas", "Regul", "Heavy"))
model1 = VGAM::vglm(Smoke ~ Height + Pulse, data=data, 
  family = cumulative(parallel=FALSE ~ Pulse), maxit=1000)
summary(model1)
dc(model1, values1 = c(150,mean(MASS::survey$Pulse,na.rm=TRUE)),
  values2 = c(151,mean(MASS::survey$Pulse,na.rm=TRUE)), type = "simulation")
# all differences are significant as the confidence intervals do not include 0
}
}

\keyword{ models }
