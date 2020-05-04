# EXAMEN FINAL

library(inspectdf)
library(dplyr)
library(leaps)
library(MASS)
library(Rcmdr)
library(Hmisc)
library(car)
library(lmtest)
library(olsrr)
library(sandwich)

# Data load
# setwd('/Users/daniel/Documents/MCD/Quantitative Analysis/bicycle-rent')
orData = read.csv('./day.csv')
data = orData
str(data)

# Data transformation
# Date
data$dteday = as.Date(as.character(data$dteday), format="%Y-%m-%d")
summary(data$dteday)

# Season
data$season = factor(data$season, levels=c(1,2,3,4), labels=c("springer","summer", "fall","winter"))
summary(data$season)

# Year (yr)
data$yr = factor(data$yr, levels=c(0,1), labels=c("2011","2012"))
summary(data$yr)

# Month (mnth)
data$mnth = factor(data$mnth, levels=c(1:12), labels=c("Jan","Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
summary(data$mnth)

# Holiday (holiday)
data$holiday = factor(data$holiday, levels=c(0,1), labels=c("Not Holiday","Holiday"))
summary(data$holiday)

# weekday
data$weekday = factor(data$weekday, levels=c(0:6), labels=c("Mon","Tue", "Wed", "Thr", "Fri", "Sat", "Sun"))
summary(data$weekday)

# workingday
data$workingday = factor(data$workingday, levels=c(0,1), labels=c("Holiday or Weekend","Workday"))
summary(data$workingday)

# weathersit
data$weathersit = factor(
  data$weathersit,
  levels=c(1:4),
  labels=c(
    "Clear-FewClouds-PartlyCloudy",
    "Mist+Cloudy-Mist+BrokenClouds-Mist+FewClouds-Mist",
    "LightSnow-LightRain+Thunderstorm+ScatteredClouds-LightRain+ScatteredClouds",
    "HeavyRain+IcePallets+Thunderstorm+Mist-Snow+Fog"
  ))
summary(data$weathersit)


# Modelo inicial - Regresion multiple, Inferencia, Variables dummy
# firstModelData = data %>% select(-dteday, -instant)
firstModelData = data[,c(-1,-2,-5,-6,-11,-14,-15)]
firstModel = lm(cnt ~ ., data = firstModelData)
summary(firstModel)

# Multicolinealidad
# Prueba VIF
vif(firstModel)

# Prueba de Belskey, Kuh y Welsh
firstModelXTX <- model.matrix(firstModel)
eFirstModel <- eigen(t(firstModelXTX) %*% firstModelXTX)

lambda.1 <- max(eFirstModel$val)
lambda.k <- min(eFirstModel$val)
kappa <- sqrt(lambda.1/lambda.k)
kappa

# Heterocedasticidad
bptest(firstModel, studentize = FALSE)

bptest(firstModel, studentize = TRUE)

ols_test_breusch_pagan(firstModel, rhs = TRUE, multiple = TRUE, p.adj = 'bonferroni')
ols_test_normality(firstModel)

coeftest(res1, vcov = (vcovHC(res1)))


# Modelo seleccionado automaticamente
# fwd.model <- regsubsets(x = firstModelData, y = data$cnt, nvmax=500, method = "forward")
fwd.model = stepwise(
  firstModel, 
  direction = c("forward")
)
length(fwd.model$coefficients)
summary(fwd.model)

bwd.model = stepwise(
  firstModel, 
  direction = c("backward")
)
length(bwd.model$coefficients)
summary(bwd.model)

bwd.fwd.model = stepwise(
  firstModel, 
  direction = c("backward/forward")
)
length(bwd.fwd.model$coefficients)
summary(bwd.fwd.model)

fwd.bwd.model = stepwise(
  firstModel, 
  direction = c("forward/backward")
)
length(fwd.bwd.model$coefficients)
summary(fwd.bwd.model)

# Multicolinealidad

# Heterocedasticidad

# Autocorrelacion

# Coeficientes estandarizados
