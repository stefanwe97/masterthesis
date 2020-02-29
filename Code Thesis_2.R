## ---------------------------------------------------------------------------##
##                                                                            ##
## Script name: Code Master Thesis 2                                          ##
##                                                                            ##  
## Author: Raphael Bierschenk, Stefan Wennemar                                ##
##                                                                            ##
## Date Created: 2018-12-09                                                   ##
##                                                                            ##
## ---------------------------------------------------------------------------##

rm(list = ls())
options(warn=-1)

start_time <- Sys.time()

library(readr)
library(ggplot2)
library(scales)
library(dplyr)
library(lubridate)
library(tidyverse)
library(zoo)
library(xts)
library(forecast)
library(tseries)
library(urca)
library(sweep)
library(broom)
library(ggthemes)
library(moments)
library(optimx)
library(stargazer)
library(rugarch)
library(sandwich)
library(lmtest)
library(extrafont)
library(RColorBrewer)
library(gridExtra)
#font_import()

# ***** Import Data *****
FF_daily <- read_csv("F-F_Research_Data_Factors_daily.CSV", col_names = TRUE, skip = 3)
FF_monthly <- read_csv("F-F_Research_Data_Factors.CSV", col_names = TRUE, skip = 3)

FF_daily <- FF_daily %>% rename(Date = X1)
FF_monthly <- FF_monthly %>% rename(Date = X1)

first_day <- 19260701
first_month <- 192607
last_day <- 20190830
last_month <- 201908

FF_daily <- FF_daily %>% subset(subset = Date <= last_day & Date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = Date <= last_month & Date >= first_month)

n_days <- as.integer(count(FF_daily))
n_months <- as.integer(count(FF_monthly))


# ***** Manipulate Data *****
FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")
FF_monthly$Date <- as.Date(FF_monthly$Date)

# Format of u: percentages as decimal numbers 
FF_daily$u <- log(1+FF_daily$`Mkt-RF`/100)
FF_daily$u_sq <- FF_daily$u^2

# Calculate EWMA Variances
# Estimate Parameters
ewma_function_daily <- function(lambda)
{
  ewma_variances <- c(1:nrow(FF_daily))
  ewma_variances[1] <- FF_daily$u_sq[1]
  for(i in 2:nrow(FF_daily)) {
    ewma_variances[i] <- lambda*ewma_variances[i-1] + (1-lambda)*FF_daily$u_sq[i]
  }
  ewma_likelihood <- c(1:(nrow(FF_daily)-1))
  for(i in 1:(nrow(FF_daily)-1)) {
    ewma_likelihood[i] <- -log(ewma_variances[i])-FF_daily$u_sq[i+1]/ewma_variances[i]
  }
  return (sum(ewma_likelihood))
}
print(ewma_max_daily <- optimize(ewma_function_daily, interval = c(0, 1), 
                                 maximum = TRUE, tol = 0.000000000000001))

# Calculate Variance
lambda = ewma_max_daily$maximum
FF_daily$EWMA_vars <- c(1:nrow(FF_daily))
FF_daily$EWMA_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$EWMA_vars[i] <- lambda*FF_daily$EWMA_vars[i-1]+(1-lambda)*FF_daily$u_sq[i]
}

# Calculate GARCH Variances
# Estimate Parameters
long_term_variance <- mean(FF_daily$u_sq)
garch_function_daily <- function(alpha, beta)
{
  omega <- max(0,long_term_variance*(1-alpha-beta))
  garch_variances <- c(1:nrow(FF_daily))
  garch_variances[1] <- FF_daily$u_sq[1]
  for(i in 2:nrow(FF_daily)) {
    garch_variances[i] <- omega + beta*garch_variances[i-1] + alpha*FF_daily$u_sq[i]
  }
  garch_likelihood <- c(1:(nrow(FF_daily)-1))
  for(i in 1:(nrow(FF_daily)-1)) {
    garch_likelihood[i] <- -log(garch_variances[i])-FF_daily$u_sq[i+1]/garch_variances[i]
  }
  return (sum(garch_likelihood))
}
print(garch_max_daily <- optimx(c(0.1, 0.9), function(x) garch_function_daily(x[1], x[2]), 
                                method = "Nelder-Mead", control = list(maximize = TRUE)))

# Calculate Variance
omega = max(0,long_term_variance*(1-alpha_daily-beta_daily))
alpha = garch_max_daily$p1
beta = garch_max_daily$p2
FF_daily$GARCH_vars <- c(1:nrow(FF_daily))
FF_daily$GARCH_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$GARCH_vars[i] <- omega+alpha*FF_daily$u_sq[i]+beta*FF_daily$GARCH_vars[i-1]
}

################################################################################
#                    Strategies with flexible time periods                     #
################################################################################

trading_days <- 264

# Generate output lists
returns_flex_list <- list()
factor_flex_list <- c(1:50)

reg_models_flex_mkt <- list()
reg_models_flex_mkt_1bps <- list()
reg_models_flex_mkt_10bps <- list()
reg_models_flex_mkt_14bps <- list()

reg_models_flex_ff3 <- list()

# Determine intervals for variances
min_days_for_var = 5
intervals <- c(min_days_for_var, 11, 22, 44, 66, 132, 264, 518)

# Create data frama that contains daily variances for all strategies
daily_vars <- data.frame(0)
daily_vars <- data.frame(FF_daily$Date[-c(1:(min_days_for_var-1))],
                         FF_daily$Mkt[-c(1:(min_days_for_var-1))],
                         FF_daily$RF[-c(1:(min_days_for_var-1))],
                         FF_daily$SMB[-c(1:(min_days_for_var-1))],
                         FF_daily$HML[-c(1:(min_days_for_var-1))])
colnames(daily_vars) <- c("Date", "Mkt", "RF", "SMB", "HML")

# For example, on the fifth day, vars are calculated for days 1 to 5
for (i in 1:nrow(daily_vars))
{
  for (interval in intervals)
  {
    col_name <- paste("Var_", interval, sep = "")
    # Calculate vars for period smaller than interval by using returns from
    # day one to day i until length of interval is reached
    if (i <= interval - min_days_for_var + 1) {
      daily_vars[[col_name]][i] <- var(FF_daily$Mkt[1:(i+min_days_for_var-1)])*(i+4-1)/(i+4)
    }
    else {
      mkt_temp <- FF_daily$Mkt[(i-interval+min_days_for_var):(i+min_days_for_var-1)]
      daily_vars[[col_name]][i] <- var(mkt_temp)*(length(mkt_temp)-1)/length(mkt_temp)
        
    }
  }
  # Set EWMA and GARCH variance
  daily_vars$EWMA[i] <- FF_daily$EWMA_vars[i+min_days_for_var-1]*10000
  daily_vars$GARCH[i] <- FF_daily$GARCH_vars[i+min_days_for_var-1]*10000
}

# ggplot(daily_vars, aes(Date)) + 
#   geom_line(aes(y=sqrt(264*EWMA), colour="EWMA")) + 
#   geom_line(aes(y=sqrt(264*GARCH), colour="GARCH")) + 
#   geom_line(aes(y=sqrt(264*Var_5), colour="Var_5")) + 
#   geom_line(aes(y=sqrt(264*Var_11), colour="Var_11")) + 
#   geom_line(aes(y=sqrt(264*Var_22), colour="Var_22")) + 
#   geom_line(aes(y=sqrt(264*Var_44), colour="Var_44")) + 
#   geom_line(aes(y=sqrt(264*Var_66), colour="Var_66")) + 
#   geom_line(aes(y=sqrt(264*Var_126), colour="Var_126")) +
#   geom_line(aes(y=sqrt(264*Var_264), colour="Var_264")) +
#   geom_line(aes(y=sqrt(264*Var_504), colour="Var_504"))

# Calculate percentage deviations of all variances
for (i in 6:15) {
  daily_vars[[paste(colnames(daily_vars[i]), "_perc_dev", sep = "")]][1] <- 0
}
for (i in 2:nrow(daily_vars)) {
  for (j in 16:25) {
    daily_vars[i,j] <- daily_vars[i,j-10]/daily_vars[i-1,j-10] - 1
  }
}

# Create a list containing quantiles of variances
quantiles <- list()
quantiles[[1]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(1/44, 1-1/44))
quantiles[[2]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.05, 0.95))
quantiles[[3]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.10, 0.90))
quantiles[[4]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.25, 0.75))
quantiles[[5]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.50, 0.50))
quantiles

strategies <- colnames(daily_vars[16:25])

# Warning: Loop runs for ~ 10 mins and requires 400 MB of RAM
for (quantile in 1:length(quantiles))
{
  for (strategy in strategies) # c("Var_264_perc_dev"))
  {
    returns_flex <- data.frame(ymd("1900/01/01"), 0, 0, 0, 0, 0)
    colnames(returns_flex) <- c("Date", "Variance", "Mkt", "RF", "SMB", "HML")
    
    last_i = 0
    for (i in 1:nrow(daily_vars)) {
      if (daily_vars[[strategy]][i] > quantiles[[quantile]][2,strategy] || 
          daily_vars[[strategy]][i] < quantiles[[quantile]][1,strategy]) {
        mkt_temp = 0
        rf_temp = 0
        smb_temp = 0
        hml_temp = 0
        for (j in c((last_i+1):i)) {
          mkt_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + mkt_temp/100) - 1)*100
          rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
          smb_temp <- ((1 + daily_vars$SMB[j]/100)*(1 + smb_temp/100) - 1)*100
          hml_temp <- ((1 + daily_vars$HML[j]/100)*(1 + hml_temp/100) - 1)*100
        }
        df_temp <- data.frame(daily_vars$Date[i],
                              daily_vars[[substr(strategy,1,nchar(strategy)-9)]][i],
                              mkt_temp, rf_temp, smb_temp, hml_temp)
        returns_flex[nrow(returns_flex) + 1,] <- df_temp
        last_i = i
      }
    }
    # if last day's variance was not outside interval, add last period manually
    if (last_i != nrow(daily_vars)) {
      mkt_temp = 0
      rf_temp = 0
      smb_temp = 0
      hml_temp = 0
      for (j in c((last_i+1):i)) {
        mkt_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + mkt_temp/100) - 1)*100
        rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
        smb_temp <- ((1 + daily_vars$SMB[j]/100)*(1 + smb_temp/100) - 1)*100
        hml_temp <- ((1 + daily_vars$HML[j]/100)*(1 + hml_temp/100) - 1)*100
      }
      df_temp <- data.frame(daily_vars$Date[i],
                            daily_vars[[substr(strategy,1,nchar(strategy)-9)]][i],
                            mkt_temp, rf_temp, smb_temp, hml_temp)
      returns_flex[nrow(returns_flex) + 1,] <- df_temp
    }
    
    returns_flex$`Mkt-RF` <- returns_flex$Mkt - returns_flex$RF
    
    returns_flex <- returns_flex[-1,]
    
    
    # Set shortcuts
    n <- nrow(returns_flex)
    denom <- returns_flex$Variance[-n]
    mktrf <- returns_flex$`Mkt-RF`[-1]
    rf <- returns_flex$RF[-1]
    
    # Calculate c
    # We do not have to convert var to sample variance because the prefactor cancels out
    a_qe_flex <- var(mktrf/denom)
    b_qe_flex <- 2*cov(mktrf/denom,rf)
    c_qe_flex <- var(rf)-var(mktrf+rf)
    
    c_flex <- 1/(2*a_qe_flex)*(-b_qe_flex+sqrt((b_qe_flex^2-4*a_qe_flex*c_qe_flex)))
    
    # Calculate weights
    weights_flex <- c_flex/denom
    quantile(weights_flex, probs = c(0.5, 0.75, 0.9, 0.99))
    
    returns_flex$Weight <- c(1,weights_flex)
    
    # Calculate volatility managed returns (incl. various transaction costs)
    returns_flex <- returns_flex %>%
      mutate(AbsDeltaW = abs(diff(c(0,Weight))),
             VMR = Weight*`Mkt-RF` + RF,
             VMR_1bps = Weight*(`Mkt-RF`) + RF - AbsDeltaW*0.01,
             VMR_10bps = Weight*(`Mkt-RF`) + RF - AbsDeltaW*0.10,
             VMR_14bps = Weight*(`Mkt-RF`) + RF - AbsDeltaW*0.14)
    
    # Check whether VMR and Mkt have the same variance
    apply(returns_flex[-1,c(3,10)], 2, var)
    
    # Calculate total returns (incl. various transaction costs)
    for (i in 1:n) {
      returns_flex$tot_ret_VM[i] <- c(1,returns_flex$tot_ret_VM)[i]*
        (1 + returns_flex$VMR[i]/100)
      returns_flex$tot_ret_VM_1bps[i] <- c(1,returns_flex$tot_ret_VM_1bps)[i]*
        (1 + returns_flex$VMR_1bps[i]/100)
      returns_flex$tot_ret_VM_10bps[i] <- c(1,returns_flex$tot_ret_VM_10bps)[i]*
        (1 + returns_flex$VMR_10bps[i]/100)
      returns_flex$tot_ret_VM_14bps[i] <- c(1,returns_flex$tot_ret_VM_14bps)[i]*
        (1 + returns_flex$VMR_14bps[i]/100)
    }
    
    # Annualize returns
    trading_days <- 264
    factor <- (nrow(returns_flex) - 1) / nrow(daily_vars) * trading_days
    returns_flex <- returns_flex %>%
      mutate(Mkt = factor*Mkt,
             RF = factor*RF,
             SMB = factor*SMB,
             HML = factor*HML,
             `Mkt-RF` = factor*`Mkt-RF`,
             VMR = factor*VMR,
             VMR_1bps = factor*VMR_1bps,
             VMR_10bps = factor*VMR_10bps,
             VMR_14bps = factor*VMR_14bps)
    
    # Regressions on Mkt factor
    reg_flex_mkt <- lm(VMR[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    reg_flex_mkt_1bps <- lm(VMR_1bps[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    reg_flex_mkt_10bps <- lm(VMR_10bps[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    reg_flex_mkt_14bps <- lm(VMR_14bps[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    
    # Regression on FF3 factors
    reg_flex_ff3 <- lm(VMR[-1] - RF[-1] ~ `Mkt-RF`[-1] + SMB[-1] + HML[-1], returns_flex)
    
    # Write return data frame and regression models to output lists
    index <- match(strategy, strategies) + 10*(quantile - 1)
    
    returns_flex_list[[index]] <- returns_flex
    factor_flex_list[index] <- factor
    
    reg_models_flex_mkt[[index]] <- reg_flex_mkt
    reg_models_flex_mkt_1bps[[index]] <- reg_flex_mkt_1bps
    reg_models_flex_mkt_10bps[[index]] <- reg_flex_mkt_10bps
    reg_models_flex_mkt_14bps[[index]] <- reg_flex_mkt_14bps
    
    reg_models_flex_ff3[[index]] <- reg_flex_ff3
    
    print(paste("Progress: ",index*2,"%", sep = ""))
  }
}
end_time <- Sys.time()
end_time - start_time

# function for showing significance stars of FF3 alpha in output table
ff3_alpha_stars <- function(model)
{
  p_val <- coeftest(model, vcovHC(model, type = "HC"))[1,4]
  stars <- ""
  if (p_val < 0.01) {
    stars <- "<sup>***</sup>"
  } else if (p_val < 0.05) {
    stars <- "<sup>**</sup>"
  } else if (p_val < 0.1) {
    stars <- "<sup>*</sup>"
  }
  return (stars)
}

################################################################################
#                                   Tables                                     #
################################################################################

# Regression table - realized variance, one year interval strategy
indices_one_year <- c(7, 17, 27, 37, 47)
returns_one_year <- list()
factors_one_year <- c(1:length(quantiles))
models_mkt_one_year <- list()
models_ff3_one_year <- list()
robust_se_mkt_one_year <- list()
robust_se_ff3_one_year <- list()
SR_one_year <- c(1:length(quantiles))
appr_ratio_one_year <- c(1:length(quantiles))
alphas_ff3_one_year <- c(1:length(quantiles))
se_ff3_one_year <- c(1:length(quantiles))

for (i in 1:length(quantiles)) {
  returns_one_year[[i]] <- returns_flex_list[[indices_one_year[i]]]
  factors_one_year[[i]] <- factor_flex_list[[indices_one_year[i]]]
  models_mkt_one_year[[i]] <- reg_models_flex_mkt[[indices_one_year[i]]]
  models_ff3_one_year[[i]] <- reg_models_flex_ff3[[indices_one_year[i]]]
  robust_se_mkt_one_year[[i]] <- sqrt(diag(vcovHC(models_mkt_one_year[[i]], type = "HC")))
  robust_se_ff3_one_year[[i]] <- sqrt(diag(vcovHC(models_ff3_one_year[[i]], type = "HC")))
  SR_one_year[i] <- round(mean(returns_one_year[[i]]$VMR - returns_one_year[[i]]$RF) /
    sd(returns_one_year[[i]]$VMR - returns_one_year[[i]]$RF) * sqrt(factors_one_year[i]),2)
  appr_ratio_one_year[i] <- round(models_mkt_one_year[[i]]$coefficients[1] /
    sqrt(mean(residuals(models_mkt_one_year[[i]])^2)) * sqrt(factors_one_year[i]),2)
  alphas_ff3_one_year[i] <- paste(round(models_ff3_one_year[[i]]$coefficients[1],2),
                               ff3_alpha_stars(models_ff3_one_year[[i]]),sep = "")
  se_ff3_one_year[i] <- paste("(",round(robust_se_ff3_one_year[[i]][1],2),")",sep = "")
}

stargazer(models_mkt_one_year, 
          se = robust_se_mkt_one_year,
          type = "text", out = "table_flex_one_year.htm",
          title = "Panel A: Interval of one year",
          column.labels = c("95.5%", "90%", "80%", "50%", "0%"),
          column.separate = c(1,1,1,1,1),
          dep.var.labels = "<i>Volatility-Managed Return",
          dep.var.caption = "Quantile",
          covariate.labels = c("Mkt-RF", "Alpha (&#945;)"),
          omit.stat = c("f", "adj.rsq"), df = FALSE, no.space = TRUE,
          table.layout = "-dl#c-t-s-a-",
          digits = 2, digits.extra = 2,
          add.lines = list(c("Vol-Managed Sharpe", SR_one_year),
                           c("Appraisal Ratio", appr_ratio_one_year),
                           c("Alpha (&#945;) FF3", alphas_ff3_one_year),
                           c("", se_ff3_one_year)))

# Regression table - realized variance, one year interval strategy
indices_one_month <- c(3, 13, 23, 33, 43)
returns_one_month <- list()
factors_one_month <- c(1:length(quantiles))
models_mkt_one_month <- list()
models_ff3_one_month <- list()
robust_se_mkt_one_month <- list()
robust_se_ff3_one_month <- list()
SR_one_month <- c(1:length(quantiles))
appr_ratio_one_month <- c(1:length(quantiles))
alphas_ff3_one_month <- c(1:length(quantiles))
se_ff3_one_month <- c(1:length(quantiles))

for (i in 1:length(quantiles)) {
  returns_one_month[[i]] <- returns_flex_list[[indices_one_month[i]]]
  factors_one_month[[i]] <- factor_flex_list[[indices_one_month[i]]]
  models_mkt_one_month[[i]] <- reg_models_flex_mkt[[indices_one_month[i]]]
  models_ff3_one_month[[i]] <- reg_models_flex_ff3[[indices_one_month[i]]]
  robust_se_mkt_one_month[[i]] <- sqrt(diag(vcovHC(models_mkt_one_month[[i]], type = "HC")))
  robust_se_ff3_one_month[[i]] <- sqrt(diag(vcovHC(models_ff3_one_month[[i]], type = "HC")))
  SR_one_month[i] <- round(mean(returns_one_month[[i]]$VMR - returns_one_month[[i]]$RF) /
                            sd(returns_one_month[[i]]$VMR - returns_one_month[[i]]$RF) * sqrt(factors_one_month[i]),2)
  appr_ratio_one_month[i] <- round(models_mkt_one_month[[i]]$coefficients[1] /
                                    sqrt(mean(residuals(models_mkt_one_month[[i]])^2)) * sqrt(factors_one_month[i]),2)
  alphas_ff3_one_month[i] <- paste(round(models_ff3_one_month[[i]]$coefficients[1],2),
                                  ff3_alpha_stars(models_ff3_one_month[[i]]),sep = "")
  se_ff3_one_month[i] <- paste("(",round(robust_se_ff3_one_month[[i]][1],2),")",sep = "")
}

stargazer(models_mkt_one_month, 
          se = robust_se_mkt_one_month,
          type = "text", out = "table_flex_one_month.htm",
          title = "Panel B: Interval of one month",
          column.labels = c("95.5%", "90%", "80%", "50%", "0%"),
          column.separate = c(1,1,1,1,1),
          dep.var.labels = "<i>Volatility-Managed Return",
          dep.var.caption = "Quantile",
          covariate.labels = c("Mkt-RF", "Alpha (&#945;)"),
          omit.stat = c("f", "adj.rsq"), df = FALSE, no.space = TRUE,
          table.layout = "-dl#c-t-s-a-",
          digits = 2, digits.extra = 2,
          add.lines = list(c("Vol-Managed Sharpe", SR_one_month),
                           c("Appraisal Ratio", appr_ratio_one_month),
                           c("Alpha (&#945;) FF3", alphas_ff3_one_month),
                           c("", se_ff3_one_month)))

# Regression table - EWMA strategy
indices_EWMA <- c(9, 19, 29, 39, 49)
returns_EWMA <- list()
factors_EWMA <- c(1:length(quantiles))
models_mkt_EWMA <- list()
models_ff3_EWMA <- list()
robust_se_mkt_EWMA <- list()
robust_se_ff3_EWMA <- list()
SR_EWMA <- c(1:length(quantiles))
appr_ratio_EWMA <- c(1:length(quantiles))
alphas_ff3_EWMA <- c(1:length(quantiles))
se_ff3_EWMA <- c(1:length(quantiles))

for (i in 1:length(quantiles)) {
  returns_EWMA[[i]] <- returns_flex_list[[indices_EWMA[i]]]
  factors_EWMA[[i]] <- factor_flex_list[[indices_EWMA[i]]]
  models_mkt_EWMA[[i]] <- reg_models_flex_mkt[[indices_EWMA[i]]]
  models_ff3_EWMA[[i]] <- reg_models_flex_ff3[[indices_EWMA[i]]]
  robust_se_mkt_EWMA[[i]] <- sqrt(diag(vcovHC(models_mkt_EWMA[[i]], type = "HC")))
  robust_se_ff3_EWMA[[i]] <- sqrt(diag(vcovHC(models_ff3_EWMA[[i]], type = "HC")))
  SR_EWMA[i] <- round(mean(returns_EWMA[[i]]$VMR - returns_EWMA[[i]]$RF) /
                             sd(returns_EWMA[[i]]$VMR - returns_EWMA[[i]]$RF) * sqrt(factors_EWMA[i]),2)
  appr_ratio_EWMA[i] <- round(models_mkt_EWMA[[i]]$coefficients[1] /
                                     sqrt(mean(residuals(models_mkt_EWMA[[i]])^2)) * sqrt(factors_EWMA[i]),2)
  alphas_ff3_EWMA[i] <- paste(round(models_ff3_EWMA[[i]]$coefficients[1],2),
                                   ff3_alpha_stars(models_ff3_EWMA[[i]]),sep = "")
  se_ff3_EWMA[i] <- paste("(",round(robust_se_ff3_EWMA[[i]][1],2),")",sep = "")
}

stargazer(models_mkt_EWMA, 
          se = robust_se_mkt_EWMA,
          type = "text", out = "table_flex_EWMA.htm",
          title = "Panel A: EWMA",
          column.labels = c("95.5%", "90%", "80%", "50%", "0%"),
          column.separate = c(1,1,1,1,1),
          dep.var.labels = "<i>Volatility-Managed Return",
          dep.var.caption = "Quantile",
          covariate.labels = c("Mkt-RF", "Alpha (&#945;)"),
          omit.stat = c("f", "adj.rsq"), df = FALSE, no.space = TRUE,
          table.layout = "-dl#c-t-s-a-",
          digits = 2, digits.extra = 2,
          add.lines = list(c("Vol-Managed Sharpe", SR_EWMA),
                           c("Appraisal Ratio", appr_ratio_EWMA),
                           c("Alpha (&#945;) FF3", alphas_ff3_EWMA),
                           c("", se_ff3_EWMA)))

# Regression table - GARCH strategy
indices_GARCH <- c(10, 20, 30, 40, 50)
returns_GARCH <- list()
factors_GARCH <- c(1:length(quantiles))
models_mkt_GARCH <- list()
models_ff3_GARCH <- list()
robust_se_mkt_GARCH <- list()
robust_se_ff3_GARCH <- list()
SR_GARCH <- c(1:length(quantiles))
appr_ratio_GARCH <- c(1:length(quantiles))
alphas_ff3_GARCH <- c(1:length(quantiles))
se_ff3_GARCH <- c(1:length(quantiles))

for (i in 1:length(quantiles)) {
  returns_GARCH[[i]] <- returns_flex_list[[indices_GARCH[i]]]
  factors_GARCH[[i]] <- factor_flex_list[[indices_GARCH[i]]]
  models_mkt_GARCH[[i]] <- reg_models_flex_mkt[[indices_GARCH[i]]]
  models_ff3_GARCH[[i]] <- reg_models_flex_ff3[[indices_GARCH[i]]]
  robust_se_mkt_GARCH[[i]] <- sqrt(diag(vcovHC(models_mkt_GARCH[[i]], type = "HC")))
  robust_se_ff3_GARCH[[i]] <- sqrt(diag(vcovHC(models_ff3_GARCH[[i]], type = "HC")))
  SR_GARCH[i] <- round(mean(returns_GARCH[[i]]$VMR - returns_GARCH[[i]]$RF) /
                        sd(returns_GARCH[[i]]$VMR - returns_GARCH[[i]]$RF) * sqrt(factors_GARCH[i]),2)
  appr_ratio_GARCH[i] <- round(models_mkt_GARCH[[i]]$coefficients[1] /
                                sqrt(mean(residuals(models_mkt_GARCH[[i]])^2)) * sqrt(factors_GARCH[i]),2)
  alphas_ff3_GARCH[i] <- paste(round(models_ff3_GARCH[[i]]$coefficients[1],2),
                              ff3_alpha_stars(models_ff3_GARCH[[i]]),sep = "")
  se_ff3_GARCH[i] <- paste("(",round(robust_se_ff3_GARCH[[i]][1],2),")",sep = "")
}

stargazer(models_mkt_GARCH, 
          se = robust_se_mkt_GARCH,
          type = "text", out = "table_flex_GARCH.htm",
          title = "Panel B: GARCH",
          column.labels = c("95.5%", "90%", "80%", "50%", "0%"),
          column.separate = c(1,1,1,1,1),
          dep.var.labels = "<i>Volatility-Managed Return",
          dep.var.caption = "Quantile",
          covariate.labels = c("Mkt-RF", "Alpha (&#945;)"),
          omit.stat = c("f", "adj.rsq"), df = FALSE, no.space = TRUE,
          table.layout = "-dl#c-t-s-a-",
          digits = 2, digits.extra = 2,
          add.lines = list(c("Vol-Managed Sharpe", SR_GARCH),
                           c("Appraisal Ratio", appr_ratio_GARCH),
                           c("Alpha (&#945;) FF3", alphas_ff3_GARCH),
                           c("", se_ff3_GARCH)))

# Alphas, Appraisal Ratios & RMSE (incl transaction costs)
output_flex_alpha_list <- list()
output_flex_appr_ratio_list <- list()
output_flex_rmse_list <- list()
for (quantile in 1:length(quantiles)) {
  names_cost <- c("| dw |", "Default", "1 bps", "10 bps", "14 bps")
  output_flex_alpha <- data.frame(matrix(ncol = length(intervals)+2, nrow = length(names_cost)))
  rownames(output_flex_alpha) <- names_cost
  colnames(output_flex_alpha) <- c("1 week", "2 weeks", "1 months", "2 months",
                                   "3 months", "6 months", "1 year", "2 years",
                                   "EWMA", "GARCH")
  output_flex_appr_ratio <- data.frame(matrix(ncol = length(intervals)+2, nrow = length(names_cost)-1))
  rownames(output_flex_appr_ratio) <- names_cost[-1]
  colnames(output_flex_appr_ratio) <- c("1 week", "2 weeks", "1 months",
                                        "2 months", "3 months", "6 months",
                                        "1 year", "2 years", "EWMA", "GARCH")
  output_flex_rmse <- data.frame(matrix(ncol = length(intervals)+2, nrow = length(names_cost)-1))
  rownames(output_flex_rmse) <- names_cost[-1]
  colnames(output_flex_rmse) <- c("1 week", "2 weeks", "1 months",
                                        "2 months", "3 months", "6 months",
                                        "1 year", "2 years", "EWMA", "GARCH")
  for (i in 1:(length(intervals)+2)) {
    index <- i + (quantile - 1)*10
    output_flex_alpha[1,i] <- mean(returns_flex_list[[index]]$AbsDeltaW[-1])
    output_flex_alpha[2,i] <- reg_models_flex_mkt[[index]]$coefficients[1]
    output_flex_alpha[3,i] <- reg_models_flex_mkt_1bps[[index]]$coefficients[1]
    output_flex_alpha[4,i] <- reg_models_flex_mkt_10bps[[index]]$coefficients[1]
    output_flex_alpha[5,i] <- reg_models_flex_mkt_14bps[[index]]$coefficients[1]
    output_flex_appr_ratio[1,i] <- reg_models_flex_mkt[[index]]$coefficients[1] /
      sqrt(mean(residuals(reg_models_flex_mkt[[index]])^2)) * sqrt(factor_flex_list[index])
    output_flex_appr_ratio[2,i] <- reg_models_flex_mkt_1bps[[index]]$coefficients[1] /
      sqrt(mean(residuals(reg_models_flex_mkt_1bps[[index]])^2)) * sqrt(factor_flex_list[index])
    output_flex_appr_ratio[3,i] <- reg_models_flex_mkt_10bps[[index]]$coefficients[1] /
      sqrt(mean(residuals(reg_models_flex_mkt_10bps[[index]])^2)) * sqrt(factor_flex_list[index])
    output_flex_appr_ratio[4,i] <- reg_models_flex_mkt_14bps[[index]]$coefficients[1] /
      sqrt(mean(residuals(reg_models_flex_mkt_14bps[[index]])^2)) * sqrt(factor_flex_list[index])
    output_flex_rmse[1,i] <- sqrt(mean(residuals(reg_models_flex_mkt[[index]])^2))
    output_flex_rmse[2,i] <- sqrt(mean(residuals(reg_models_flex_mkt_1bps[[index]])^2))
    output_flex_rmse[3,i] <- sqrt(mean(residuals(reg_models_flex_mkt_10bps[[index]])^2))
    output_flex_rmse[4,i] <- sqrt(mean(residuals(reg_models_flex_mkt_14bps[[index]])^2))
  }
  output_flex_alpha_list[[quantile]] <- output_flex_alpha
  output_flex_appr_ratio_list[[quantile]] <- output_flex_appr_ratio
  output_flex_rmse_list[[quantile]] <- output_flex_rmse
}

stargazer(output_flex_alpha_list, type = "text",
          summary = c(FALSE,FALSE,FALSE,FALSE,FALSE),
          out = "alpha_flex.htm", digits = 2)
stargazer(output_flex_appr_ratio_list, type = "text",
          summary = c(FALSE,FALSE,FALSE,FALSE,FALSE),
          out = "appr_ratio_flex.htm", digits = 2)
stargazer(output_flex_rmse_list, type = "text",
          summary = c(FALSE,FALSE,FALSE,FALSE,FALSE),
          out = "rmse.htm", digits = 2)


################################################################################
#                                   Graphs                                     #
################################################################################

# General settings
display.brewer.pal(9, name = "Greys")
display.brewer.pal(9, name = "Blues")
grey_col <- brewer.pal(9, name = "Greys")
blue_col <- brewer.pal(9, name = "Blues")
line_size <- 0.7
loadfonts(device = "win")
windowsFonts(Times = windowsFont("TT Times New Roman"))

# Replicate base method to show in graphs
# Calculate Monthly Variances
monthly_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(`Mkt-RF`)*(length(`Mkt-RF`)-1)/length(`Mkt-RF`))
plot(sqrt(monthly_vars$variance*trading_days), type = "l")

denom_base <- monthly_vars$variance[-n_months]

# Calculate c
# We do not have to convert var to sample variance because the prefactor cancels out
a_qe_base <- var(FF_monthly$`Mkt-RF`[-1]/denom_base)
b_qe_base <- cov(FF_monthly$`Mkt-RF`[-1]/denom_base, FF_monthly$RF[-1])
c_qe_base <- var(FF_monthly$RF[-1])-var(FF_monthly$Mkt[-1])
c_base <- 1/(a_qe_base)*(-b_qe_base+sqrt((b_qe_base)^2-a_qe_base*c_qe_base))

# Calculate weights and volatility managed returns
weights <- c_base / denom_base

returns <- data.frame(matrix(ncol = 4, nrow = n_months - 1))
colnames(returns) <- c("Date", "Mkt", "rf", "Base")
returns <- returns %>% mutate(Date = FF_monthly$Date[-1], Mkt = FF_monthly$Mkt[-1], 
                              rf = FF_monthly$RF[-1])
for (month in 1:(n_months-1)) {
  returns$Base[month] <- weights[month] * (returns$Mkt[month] - returns$rf[month]) +
    returns$rf[month]
}

apply(returns[,c(2,4)], 2, var)

# Calculate total returns
tot_ret <- data.frame(matrix(ncol = 3, nrow = n_months))
colnames(tot_ret) <- c("Date", "Mkt", "Base")
tot_ret <- tot_ret %>% mutate(Date = FF_monthly$Date)

tot_ret[1,-1] <- 1

for (month in 2:n_months) {
  tot_ret$Mkt[month] <- tot_ret$Mkt[month-1] * 
    (1 + (returns$Mkt[month-1]/100))
  tot_ret$Base[month] <- tot_ret$Base[month-1] * 
    (1 + (returns$Base[month-1]/100))
}

# Figure 1: Distribution of GARCH percentage deviations
daily_vars %>%
  mutate(GARCH_perc_dev_temp = ifelse(GARCH_perc_dev > 1, 1, GARCH_perc_dev)) %>%
  ggplot(aes(GARCH_perc_dev_temp)) +
  geom_histogram(aes(log(1+GARCH_perc_dev_temp)), binwidth = 0.005) +
  geom_vline(xintercept = quantiles[[1]][1,10], color = "black", linetype = "dashed") +
  geom_vline(xintercept = quantiles[[1]][2,10], color = "black", linetype = "dashed") +
  ggtitle("Distribution of percentage deviations of variances") + xlab("") + ylab("") +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Figure 2: Total returns Market, Moreira and Muir, Var 1 year, EWMA, GARCH, 50% quantile
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot() +
  geom_line(aes(y=tot_ret$Mkt, x=tot_ret$Date, colour="Buy and Hold"), size = line_size) +
  geom_line(aes(y=tot_ret$Base, x=tot_ret$Date, colour="Base Strategy"), size = line_size) +
  geom_line(aes(y=returns_flex_list[[37]]$tot_ret_VM,
                x=returns_flex_list[[37]]$Date, colour="Var 1 year"), size = line_size) +
  geom_line(aes(y=returns_flex_list[[39]]$tot_ret_VM,
                x=returns_flex_list[[39]]$Date, colour="EWMA"), size = line_size) +
  geom_line(aes(y=returns_flex_list[[40]]$tot_ret_VM,
                x=returns_flex_list[[40]]$Date, colour="GARCH"), size = line_size) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks('log10', function(x) 10^x),
                     minor_breaks = scale,
                     labels = trans_format('log10', math_format(10^.x)),
                     limits = c(0.1,100000),
                     expand = c(0,0)) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Cumulative Performance") + xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = grey_col[9],
                                "Base Strategy" = grey_col[8],
                                "Var 1 year" = grey_col[7],
                                "EWMA" = grey_col[6],
                                "GARCH" = grey_col[5]),
                     breaks = c("Buy and Hold", "Base Strategy",
                                "Var 1 year", "EWMA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# Figure 3: Alphas for one month, one year interval, EWMA and GARCH strategy by quantile
alphas_figure_3 <- list(0,0,0,0)
for (i in 1:4) {
  j <- c(3,7,9,10)[i]
  alphas_figure_3[[i]][1] <- output_flex_alpha_list[[1]][2,j]
  alphas_figure_3[[i]][2] <- output_flex_alpha_list[[2]][2,j]
  alphas_figure_3[[i]][3] <- output_flex_alpha_list[[3]][2,j]
  alphas_figure_3[[i]][4] <- output_flex_alpha_list[[4]][2,j]
  alphas_figure_3[[i]][5] <- output_flex_alpha_list[[5]][2,j]
}
ggplot() +
  geom_line(aes(y=alphas_figure_3[[1]], x=c(22,10,5,2,1), colour="one month"), size = line_size) +
  geom_point(aes(y=alphas_figure_3[[1]], x=c(22,10,5,2,1), colour="one month")) +
  geom_line(aes(y=alphas_figure_3[[2]], x=c(22,10,5,2,1), colour="one year"), size = line_size) +
  geom_point(aes(y=alphas_figure_3[[2]], x=c(22,10,5,2,1), colour="one year")) +
  geom_line(aes(y=alphas_figure_3[[3]], x=c(22,10,5,2,1), colour="EWMA"), size = line_size) +
  geom_point(aes(y=alphas_figure_3[[3]], x=c(22,10,5,2,1), colour="EWMA")) +
  geom_line(aes(y=alphas_figure_3[[4]], x=c(22,10,5,2,1), colour="GARCH"), size = line_size) +
  geom_point(aes(y=alphas_figure_3[[4]], x=c(22,10,5,2,1), colour="GARCH")) +
  scale_x_continuous(limits = c(0,25), breaks = c(22,10,5,2,1),  minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,6), breaks = c(0:6), minor_breaks = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Alphas for selected strategies") + xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("one month" = blue_col[9],
                                "one year" = blue_col[8],
                                "EWMA" = blue_col[7],
                                "GARCH" = blue_col[6]),
                     breaks = c("one month", "one year",
                                "EWMA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))


# Figure 4: Average absolute weights by quantile over interval
ggplot() +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[1]][1,c(1:8)]), x=intervals,
                colour="95.5% quantile"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[1]][1,c(1:8)]), x=intervals,
                colour="95.5% quantile")) +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[2]][1,c(1:8)]), x=intervals,
                colour="90% quantile"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[2]][1,c(1:8)]), x=intervals,
                colour="90% quantile")) +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[3]][1,c(1:8)]), x=intervals,
                colour="80% quantile"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[3]][1,c(1:8)]), x=intervals,
                colour="80% quantile")) +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[4]][1,c(1:8)]), x=intervals,
                colour="50% quantile"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[4]][1,c(1:8)]), x=intervals,
                colour="50% quantile")) +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[5]][1,c(1:8)]), x=intervals,
                colour="0% quantile"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[5]][1,c(1:8)]), x=intervals,
                colour="0% quantile")) +
  scale_x_log10(breaks = intervals, minor_breaks = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Average absolute change in weights by quantile over interval") + xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("95.5% quantile" = blue_col[9],
                                "90% quantile" = blue_col[8],
                                "80% quantile" = blue_col[7],
                                "50% quantile" = blue_col[6],
                                "0% quantile" = blue_col[5]),
                     breaks = c("95.5% quantile", "90% quantile",
                                "80% quantile", "50% quantile", "0% quantile")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# Figure 5: Alphas for one month, one year interval, EWMA and GARCH strategy for
# the 0% quantile considering transaction costs
trading_costs <- c(0.00, 0.01, 0.10, 0.14)
ggplot() +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,3]), x=trading_costs,
                colour="one month"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,3]), x=trading_costs,
                colour="one month")) +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,7]), x=trading_costs,
                colour="one year"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,7]), x=trading_costs,
                colour="one year")) +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,9]), x=trading_costs,
                colour="EWMA"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,9]), x=trading_costs,
                colour="EWMA")) +
  geom_line(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,10]), x=trading_costs,
                colour="GARCH"), size = line_size) +
  geom_point(aes(y=as.numeric(output_flex_alpha_list[[4]][-1,10]), x=trading_costs,
                colour="GARCH")) +
  scale_x_continuous(breaks = trading_costs, minor_breaks = NULL) +
  scale_y_continuous(limits = c(0,6), breaks = c(0:6), minor_breaks = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Alphas for selected strategies considering transaction costs") + xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("one month" = blue_col[9],
                                "one year" = blue_col[8],
                                "EWMA" = blue_col[7],
                                "GARCH" = blue_col[6]),
                     breaks = c("one month", "one year",
                                "EWMA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))


# Figure 6: Example for when better to exclude last returns: April 1927
# Monthly variance of the whole month
var(FF_daily$Mkt[225:249])*24/25*25
# Monthly variance excluding the last six days
var(FF_daily$Mkt[225:243])*24/25*25
# Return of May 1927: 5.74%
fig12 <- ggplot() +
  geom_line(aes(y=as.numeric(FF_daily$Mkt[225:249]), x=as.Date(FF_daily$Date[225:249]),
                colour="Mkt"), size = line_size) +
  scale_x_date() +
  scale_y_continuous(limits = c(-2,2), breaks = c(-2:2), minor_breaks = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Market returns in April 1927") + xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Mkt" = blue_col[9]),
                     breaks = c("Mkt"))

# Figure 7: Example for when better to exclude first returns: February 1940
# Monthly variance of the whole month
var(FF_daily$Mkt[4052:4074])*24/25*25
# Monthly variance excluding the first seven days
var(FF_daily$Mkt[4059:4074])*24/25*25
fig13 <- ggplot() +
  geom_line(aes(y=as.numeric(FF_daily$Mkt[4052:4074]), x=as.Date(FF_daily$Date[4052:4074]),
                colour="Mkt"), size = line_size) +
  scale_x_date() +
  scale_y_continuous(limits = c(-2,2), breaks = c(-2:2), minor_breaks = NULL) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Market returns in February 1940") + xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Mkt" = blue_col[9]),
                     breaks = c("Mkt"))

# Plot both next to each other
grid.arrange(fig12, fig13, ncol = 2)


# Figure 8: Example for when variance keeps the same for period longer than a month
trailing_vars <- c(1:153)
for (i in trailing_vars) {
  trailing_vars[i] <- var(FF_daily$Mkt[(10351+i-1):(10351+21+i-1)])*21/22*22
}
ggplot() +
  geom_line(aes(y=as.numeric(FF_daily$Mkt[10351:10525]), x=as.Date(FF_daily$Date[10351:10525]),
                colour="Mkt"), size = line_size) +
  geom_line(aes(y=as.numeric(trailing_vars)/8-4, x=as.Date(FF_daily$Date[10373:10525]), 
                colour="22-day trailing Variance"), size = 0.8) +
  scale_x_date() +
  scale_y_continuous(limits = c(-4,4), breaks = c(-4:4), minor_breaks = NULL, name = "Return in %",
                     sec.axis = sec_axis(~.*8+32, name = "Monthly Variance")) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Market Returns and their Variance in 1963") + xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Mkt" = blue_col[6],
                                "22-day trailing Variance" = blue_col[9]),
                     breaks = c("Mkt", "22-day trailing Variance"))

