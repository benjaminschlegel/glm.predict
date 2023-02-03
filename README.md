# glm.predict
R GNU Package to simulate/bootstrap predicted values/probabilies and discrete changes for lm(), glm(), glm.nb(), polr(), multinom(), tobit() and lmer() models. 

## Example

Here an example how you can use the function `predicts()`. First we estimate an logistic regression to predict the gender of a person using hight, smoking and pulse as independent variables.

`model1 = glm(Sex ~ Height + Smoke + Pulse, data=MASS::survey, family=binomial(link=logit))`

Now we can estimate predicted probabilies. For Height we take 150, 170 and 190, for Smoke we take all levels of the factor and for Pulse we take the average value of the data. We set a seed to always get the same results.

```
library(glm.predict)
predicts(model1, "150-190,20;all;mean", set.seed = 1848)
```

By setting position to the position of a variable, we can get discrete changes for that specific variable. If we want to get the difference between 150 and 170 body size, we can do that in the following way:

`predicts(model1, "150,170;all;mean", set.seed = 1848, position = 1)`

To get all possible values and other parameters you can read the documentation:

`?predicts()`

## Installation

If you want to install the stable version from CRAN you can do that with the following command:

`install.packages("glm.predict")`
  
If you want to install the newest developer version from GitHub, you can do that with the package remotes:

```
install.packages("remotes")
remotes::install_github("benjaminschlegel/glm.predict")
```

You might need to install Rtools (for Windows) or XCode (for Mac OS) to install the package from GitHub.
