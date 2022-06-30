library(ggplot2)
library(tidyverse)
library(dplyr)
library(optmatch)
library(sm)
library(MatchIt)
library(survey)
library(lmtest) 
library(sandwich) 
library(kableExtra)
library(rpart)
library(boot)
library(e1071)
library(rlang)
library(dplyr)
library(RItools)

#read in file
df_ch =  read_csv("df_ch.csv")

#Boxplot
df_ch %>%
  ggplot(aes(x=  as.character(Chol_Healthy), y =Diastolic_BP)) +
  geom_boxplot(notch=TRUE) +
  labs(x="Cholesterol Ratio Healthy", y="Diastolic Blood Pressure") +
  ggtitle("Chol. Ratio Healthy vs. Diastolic B.P.")

# create value labels
attach(df_ch)
ch.f <- factor(Chol_Healthy, levels= c(0,1),
               labels = c("No","Yes"))

# plot densities
sm.density.compare(df_ch$Diastolic_BP,df_ch$Chol_Healthy, 
                   xlab="Diastolic Blood Pressure", lwd=2)
title(main="Density Function of Diastolic Blood Pressure")
# add legend via mouse click
colfill<-c(2:(2+length(levels(ch.f))))
legend(locator(1), levels(ch.f), fill=colfill, cex=1, 
       title = "Healthy Cholesterol Ratio?", y.intersp=0.6, bty="n")


###############################
## Propensity score matching ##
###############################
# Logistic regression to fit a propensity score model
prop.score <- glm(Chol_Healthy ~ Age + Gender + Race_Eth + Education +
                    Married + Citizen + Phys_Activity +  HH_Income 
                  + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                    Cholesterol_Daily + Fat_Daily + 
                    Vig_Activity + BMI, data = df_ch)

#full matching on propensity scores
ps_full_match <- fullmatch(prop.score, data = df_ch, max.controls = 10,
                           min.controls = 0.1)
summary(ps_full_match)

#plot balancing on covariates before and after matching based on propensity scores
ps_full_match_balance <- xBalance(Chol_Healthy ~ Age + Gender + Race_Eth + Education +
                                    Married + Citizen + Phys_Activity + HH_Income +
                                    + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                                    Cholesterol_Daily + Fat_Daily +
                                    Vig_Activity + BMI +
                                    strata(ps_full_match), data = df_ch,
                                  report = c("adj.means", "std.diffs", 
                                             "z.scores", "chisquare.test"))
plot(ps_full_match_balance, main ="Balance Before/After 
     Full Matching (Cholest. Ratio)")

#chi square before and after matching 
print(ps_full_match_balance,digits = 1)

#######################################################################
#######Using MatchIt to full match on propensity scores (another method)

result_full_ps <- matchit(Chol_Healthy ~ Age + Gender + Race_Eth + Education +
                            Married + Citizen + Phys_Activity + HH_Income +
                            + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                            Cholesterol_Daily + Fat_Daily +
                            Vig_Activity + BMI, 
                          data = df_ch, method = "full", 
                          distance = "glm", link = "logit",max.controls=10,
                          min.controls=0.1)
summary(result_full_ps)
plot(summary(result_full_ps))
m.data <- match.data(result_full_ps)

match.result_full_ps <- as.factor(result_full_ps$subclass)
xBalance(Chol_Healthy ~Age + Gender + Race_Eth + Education +
           Married + Citizen + Phys_Activity + HH_Income +
           + Carb_Daily + Alcohol_Daily + Protein_Daily + 
           Cholesterol_Daily + Fat_Daily +
           Vig_Activity + BMI, 
         data = df_ch, 
         strata = match.result_full_ps, 
         report = "chisquare.test")

#distribution of propensity scores across treatment/control units
plot(result_full_ps, type = "jitter", interactive = FALSE)


#########################
#Doubly Robust Estimators
#########################

#produce propensity score model and extract fitted probabilities
ps_mod <- glm(Chol_Healthy ~ Age + Gender + Race_Eth + Education +
                Married + Citizen + Phys_Activity + HH_Income +
                + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                Cholesterol_Daily + Fat_Daily +
                Vig_Activity + BMI, data=df_ch, family = "binomial")
e_hat <- ps_mod$fitted.values

#create alternate copies of data set 
#within which all subjects are treatment
#or all subjects are control (all other values held constant)
dat_1 <- transform(df_ch, Chol_Healthy  = 1)
dat_0 <- transform(df_ch, Chol_Healthy  = 0)

#generate linear regression model to predict DBP
out_mod <- lm(Diastolic_BP ~ Chol_Healthy + Age + Gender + Race_Eth + Education +
                Married + Citizen + Phys_Activity + HH_Income +
                + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                Cholesterol_Daily + Fat_Daily +
                Vig_Activity + BMI, data=df_ch)

#apply linear regression model to dat_0 and dat_1
y_hat_1 <- predict(out_mod, newdata = dat_1) 
y_hat_0 <- predict(out_mod, newdata = dat_0)
boxplot(y_hat_1, y_hat_0)

#create doubly robust estimator: 
DR_1 <- (e_hat^-1)*df_ch$Diastolic_BP*df_ch$Chol_Healthy - (e_hat^-1)*(y_hat_1*(Chol_Healthy  - e_hat))

DR_0 <- (1 - e_hat)^-1*(df_ch$Diastolic_BP*(1 - df_ch$Chol_Healthy ))+ (1 - e_hat)^-1*(y_hat_0*(df_ch$Chol_Healthy  - e_hat))

# estimated treatment effect (ATE) of having healthy (rather than unhealthy)
#cholesterol ratio
mean(DR_1) - mean(DR_0) 


###define functions to use three different regression methods to generate DRE

#### linear regression function 
function_lm <- function(data, i){
  d2 <- data[i,]
  ps_mod <- glm(Chol_Healthy ~ Age + Gender + Race_Eth + Education +
                  Married + Citizen + Phys_Activity + HH_Income +
                  + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                  Cholesterol_Daily + Fat_Daily +
                  Vig_Activity + BMI, data = d2, family = "binomial")
  e_hat <- ps_mod$fitted.values
  dat_1 <- transform(d2[d2$Chol_Healthy==1,], Chol_Healthy = 1)
  dat_0 <- transform(d2[d2$Chol_Healthy==1,], Chol_Healthy = 0)
  out_mod <- lm(Diastolic_BP ~ Chol_Healthy + Age + Gender + Race_Eth + Education +
                  Married + Citizen + Phys_Activity + HH_Income +
                  + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                  Cholesterol_Daily + Fat_Daily +
                  Vig_Activity + BMI, data = d2)
  y_hat_1 <- predict(out_mod, newdata = dat_1) 
  y_hat_0 <- predict(out_mod, newdata = dat_0)
  DR_1 <- (e_hat^-1)*d2$Diastolic_BP*d2$Chol_Healthy - (e_hat^-1)*(y_hat_1*(d2$Chol_Healthy - e_hat))
  DR_0 <- (1 - e_hat)^-1*(d2$Diastolic_BP*(1 - d2$Chol_Healthy)) + (1 - e_hat)^-1*(y_hat_0*(d2$Chol_Healthy - e_hat))
  treateffect <- mean(DR_1) - mean(DR_0)
  return(treateffect)
}

#Support Vector Machine function
function_svm <- function(data, i){
  d2 <- data[i,]
  ps_mod <- glm(Chol_Healthy ~ Age + Gender + Race_Eth + Education +
                  Married + Citizen + Phys_Activity + HH_Income +
                  + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                  Cholesterol_Daily + Fat_Daily +
                  Vig_Activity + BMI, data = d2, family = "binomial")
  e_hat <- ps_mod$fitted.values
  dat_1 <- transform(d2[d2$Chol_Healthy==1,], Chol_Healthy = 1)
  dat_0 <- transform(d2[d2$Chol_Healthy==1,], Chol_Healthy = 0)
  modelsvm = svm(Diastolic_BP ~ Chol_Healthy+ Age + Gender + Race_Eth + Education +
                   Married + Citizen + Phys_Activity + HH_Income +
                   + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                   Cholesterol_Daily + Fat_Daily +
                   Vig_Activity + BMI, data=d2)
  y_hat_1 <- predict(modelsvm, newdata = dat_1) 
  y_hat_0 <- predict(modelsvm, newdata = dat_0)
  DR_1 <- (e_hat^-1)*d2$Diastolic_BP*d2$Chol_Healthy - (e_hat^-1)*(y_hat_1*(d2$Chol_Healthy - e_hat))
  DR_0 <- (1 - e_hat)^-1*(d2$Diastolic_BP*(1 - d2$Chol_Healthy)) + (1 - e_hat)^-1*(y_hat_0*(d2$Chol_Healthy - e_hat))
  treateffect <- mean(DR_1) - mean(DR_0)
  return(treateffect)
}

#Regression Tree function
function_regtree <- function(data, i){
  d2 <- data[i,]
  ps_mod<- glm(Chol_Healthy ~ Age + Gender + Race_Eth + Education +
                 Married + Citizen + Phys_Activity + HH_Income +
                 + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                 Cholesterol_Daily + Fat_Daily +
                 Vig_Activity + BMI, data = d2, family = "binomial")
  e_hat <- ps_mod$fitted.values
  dat_1 <- transform(d2[d2$Chol_Healthy==1,], Chol_Healthy = 1)
  dat_0 <- transform(d2[d2$Chol_Healthy==1,], Chol_Healthy = 0)
  regtree= rpart(formula = Diastolic_BP ~ Chol_Healthy + Age + Gender + Race_Eth + Education +
                   Married + Citizen + Phys_Activity + HH_Income +
                   + Carb_Daily + Alcohol_Daily + Protein_Daily + 
                   Cholesterol_Daily + Fat_Daily +
                   Vig_Activity + BMI, data=d2, method = "anova")
  y_hat_1 <- predict(regtree, newdata = dat_1) 
  y_hat_0 <- predict(regtree, newdata = dat_0)
  DR_1 <- (e_hat^-1)*d2$Diastolic_BP*d2$Chol_Healthy  - (e_hat^-1)*(y_hat_1*(d2$Chol_Healthy - e_hat))
  DR_0 <- (1 - e_hat)^-1*(d2$Diastolic_BP*(1 - d2$Chol_Healthy)) + (1 - e_hat)^-1*(y_hat_0*(d2$Chol_Healthy - e_hat))
  treateffect <- mean(DR_1) - mean(DR_0)
  return(treateffect)
}


#bootstrapping (1000 resampled samples)
bootstrap_rt_te <- boot(df_ch,function_regtree,R=1000, sim = "balanced")
bootstrap_svm_te <- boot(df_ch,function_svm,R=1000, sim = "balanced")
bootstrap_lm_te <- boot(df_ch,function_lm,R=1000, sim="balanced")

#bootstrapped treatment effects and std. errors
bootstrap_svm_te
bootstrap_lm_te
bootstrap_rt_te
