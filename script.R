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
library(AER)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(relaimpo)
library(lm.beta)
library(tseries)
library(car)
# ----------------- HELPER FUNCTIONS
getKappa = function(model){
  modelXTX <- model.matrix(model)
  eModel <- eigen(t(modelXTX) %*% modelXTX)
  
  lambda.1 <- max(eModel$val)
  lambda.k <- min(eModel$val)
  kappa <- sqrt(lambda.1/lambda.k)
  kappa
}

removeVIF <- function(modelo, u){
  require(car)
  data <- modelo$model
  all_vifs <- car::vif(modelo)
  names_all <- rownames(all_vifs)
  dep_var <- all.vars(formula(modelo))[1]
  while(any(all_vifs > u)){
    var_max_vif <- names(which(all_vifs == max(all_vifs)))
    names_all <- names_all[!(names_all) %in% var_max_vif]
    myForm <- as.formula(paste(paste(dep_var, "~ "),
                               paste (names_all, collapse=" + "), sep=""))
    modelo.prueba <- lm(myForm, data= data)
    all_vifs <- car::vif(modelo.prueba) 
  }
  modelo.limpio <- modelo.prueba
  return(modelo.limpio)
}


# ----------------- DATA LOADING
# setwd('/Users/daniel/Documents/MCD/Quantitative Analysis/bicycle-rent')
orData = read.csv('./day.csv')
data = orData
str(data)

# ----------------- Data transformation
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
    "ClearFewCloudsPartlyCloudy",
    "MistCloudyMistBrokenCloudsMistFewCloudsMist",
    "LightSnowLightRainThunderstormScatteredCloudsLightRainScatteredClouds",
    "HeavyRainIcePalletsThunderstormMisSnowFog"
  ))
summary(data$weathersit)


# ----------------- Modelo inicial
firstModelData = data[,c(-1,-2,-3,-8,-11,-14,-15)]
# firstModelData = data[,c(-1,-2,-3,-8,-11,-14,-15)]
firstModel = lm(cnt ~ ., data = firstModelData)
length(firstModel$coefficients)
summary(firstModel)
removeVIF(firstModel,8)
# Multicolinealidad
# Prueba VIF
vifs = vif(firstModel)
rownames(vifs)

# Prueba de Belskey, Kuh y Welsh
getKappa(firstModel)


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

jtest(bwd.model, fwd.model)
anova(bwd.model, fwd.model)

tab_model(fwd.model, bwd.model, fwd.bwd.model,
          show.ci = FALSE, 
          dv.labels = c("Reduced Forward Model", "Backward Model"),
          transform = NULL
)

getKappa(bwd.model)

# Heterocedasticidad
# Plot the errors against the data
residuals <-resid(firstModel)
attach(firstModelData)
# par(mfrow=c(2,2))
# plot(season, residuals)
# plot(yr, residuals)
# plot(hum, residuals)
# plot(temp, residuals)
# plot(weathersit, residuals)
# plot(weekday, residuals)
# plot(windspeed, residuals)
# plot(workingday, residuals)
# plot(yr, residuals)


bptest(firstModel, studentize = FALSE)

ols_test_breusch_pagan(firstModel, rhs = TRUE, multiple = TRUE, p.adj = 'bonferroni')
ols_test_normality(firstModel)

bptest(fwd.model, studentize = TRUE)

calc.relimp(fwd.model, type = "betasq")


modelo2.s <- lm.beta(fwd.model)
coeftest(modelo2.s, vcov = (vcovHC(fwd.model))) 

# Autocorrelacion

errors = residuals(fwd.model)
sig_errors = factor(errors > 0)
summary(sig_errors)
runs.test(sig_errors)

dwtest(fwd.model, alternative = "two.sided")


coeftest(fwd.model, vcov = NeweyWest(fwd.model))
coeftest(fwd.model, vcov = kernHAC(fwd.model))
# Coeficientes estandarizados
