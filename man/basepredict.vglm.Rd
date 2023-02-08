\name{basepredict.vglm}
\alias{basepredict.vglm}

\title{
predicted value
}
\description{
The function calculates the predicted value with the confidence interval. It can be used for any
vglm model.
}
\usage{
\method{basepredict}{vglm}(model, values, sim.count=1000, conf.int=0.95, sigma=NULL, set.seed=NULL,
  type = c("any", "simulation", "bootstrap"), summary = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the model Object generated with vglm()
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
OPTIONAL the variance-covariance matrix, can be changed when having for example robust or clustered vcov. default: vcov(model)
}
  \item{set.seed}{
OPTIONAL set a seed for the random number generator
}
  \item{type}{
OPTIONAL choose between simulation and bootstrap, "any" chooses between those two
according to the number of cases (bootstrap if n < 1000). Note: bootstrap is very slow for vglm() models.
}
  \item{summary}{
OPTIONAL if mean/quantiles should be return or all simulated values (default: TRUE)	
}
}
\details{
  The function makes a simulation for the two cases and compares them to each other. 
}
\value{
  The output is a matrix have in the first column the predicted value, in the second
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
# predicted probability of smoking of a person with height 170 and an average pulse
basepredict(model1, c(170,mean(MASS::survey$Pulse,na.rm=TRUE)))
}
}

\keyword{ models }
