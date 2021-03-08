\name{basepredict.mlogit}
\alias{basepredict.mlogit}

\title{
predicted value
}
\description{
The function calculates the predicted value with the confidence interval. It can be used for any
mlogit model.
}
\usage{
\method{basepredict}{mlogit}(model, values, sim.count=1000, conf.int=0.95, sigma=NULL, set.seed=NULL,
  type = c("any", "simulation", "bootstrap"), summary = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the model Object generated with mlogit()
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
type is ignored as only simulation is implemented
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
# devtools::install_github("benjaminschlegel/schlegel")
df_selects = schlegel::selects2015
df_selects_withoutNA = df_selects %>% 
  filter(vote_choice != "other") %>% 
  mutate(vote_choice = factor(vote_choice)) %>% 
  select(age, gender, vote_choice, starts_with("lr_")) %>% 
  na.omit()
mlogit_data = dfidx::dfidx(df_selects_withoutNA, varying = 5:11, 
                    sep = "_", shape = "wide", 
                    choice = "vote_choice")

mlogit_data$distance = abs(mlogit_data$lr - mlogit_data$lr_self)

model1 = mlogit::mlogit(vote_choice ~ distance | lr_self + 
                          gender, data = mlogit_data)
summary(model1)
# predicted probability of a left male person with a distance of 2
basepredict(model1, c(1, 2, 0, 0))
}
}

\keyword{ models }
