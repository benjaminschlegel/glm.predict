\name{dc.mlogit}
\alias{dc.mlogit}

\title{
predicted values and discrete change
}
\description{
The function calculates the predicted values and the difference of
two cases with the confidence interval. It can be used for a
mlogit model.
}
\usage{
\method{dc}{mlogit}(model, values = NULL, sim.count = 1000, conf.int = 0.95, 
  sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
  type = c("any", "simulation", "bootstrap"), summary = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the model-Object generated with mlogit()
}
  \item{values}{
the values of case 1 and 2 as list in the order how they appear in the summary(model) Estimate. Values is if values1 and values2 are specified after each other in the same vector. Either values or values1 and values2 have to be specified.
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
the values of case 1 as list in the order how they appear in the summary(model) Estimate. Has to be defined if values is not defined. Can be single values or a vector of values with the length of choices - 1 for condition variables.
}
  \item{values2}{
the values of case 2 as list in the order how they appear in the summary(model) Estimate. Has to be defined if values is not defined. Can be single values or a vector of values with the length of choices - 1 for condition variables.
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
  The output is a matrix have in the first column the predicted values, in the second
  column the lower value of the confidence interval and in the third column the upper value
  of the confidence interval.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
\dontrun{
df_selects_withoutNA = selects2015 |>
  filter(vote_choice != "other") |> 
  mutate(vote_choice = factor(vote_choice)) |>
  select(age, gender, vote_choice, starts_with("lr_")) |>
  na.omit()
mlogit_data = dfidx::dfidx(df_selects_withoutNA, varying = 5:11, 
                    sep = "_", shape = "wide", 
                    choice = "vote_choice")

mlogit_data$distance = abs(mlogit_data$lr - mlogit_data$lr_self)

model1 = mlogit::mlogit(vote_choice ~ distance | lr_self + 
                          gender, data = mlogit_data)
summary(model1)
# predicted probability of a left male person with a distance of 2
dc(model1, list(1, c(2, 5, 5, 5, 5, 5), 0, 0, 1, c(2, 5, 5, 5, 5, 5), 10, 0))
}
}

\keyword{ models }
