## ---------------------------------------------------------------------------##
##                                                                            ##
## Script name: Code Master Thesis 1                                          ##
##                                                                            ##  
## Author: Raphael Bierschenk, Stefan Wennemar                                ##
##                                                                            ##
## Date Created: 2018-12-09                                                   ##
##                                                                            ##
## ---------------------------------------------------------------------------##

# If first time: set working directory and run command:
# font_import()

# Set Up: Initiate Packages
rm(list = ls())

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
library(cowplot)
library(gridExtra)

# Import Data
FF_daily <- read_csv("F-F_Research_Data_Factors_daily.CSV", col_names = TRUE, skip = 3)
FF_monthly <- read_csv("F-F_Research_Data_Factors.CSV", col_names = TRUE, skip = 3)
recession <- read_csv("USREC.csv", col_names = TRUE)

# Rename and Crop Data
FF_daily <- FF_daily %>% rename(date = X1, "mkt-rf" = "Mkt-RF", rf = RF)
FF_monthly <- FF_monthly %>% rename(date = X1, "mkt-rf" = "Mkt-RF", rf = RF)
recession <- recession %>% rename(date = DATE, indicator = USREC)

first_day <- 19260701
first_month <- 192607
last_day <- 20190831
last_month <- 201908

FF_daily <- FF_daily %>% subset(subset = date <= last_day & date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = date <= last_month & date >= first_month)
recession <- recession %>% subset(subset = date <= ymd(last_day) & date >= ymd(first_day))

# Mutate Data
n_days <- FF_daily %>% count() %>% as.numeric()
n_months <- FF_monthly %>% count() %>% as.numeric()

FF_monthly <- FF_monthly %>% mutate(mkt = `mkt-rf` + rf)
FF_daily <- FF_daily %>% mutate(mkt = `mkt-rf` + rf)

FF_daily$date <- ymd(FF_daily$date)
FF_monthly$date <- as.character(FF_monthly$date)
FF_monthly$date <- parse_date_time(FF_monthly$date, "ym")
FF_monthly$date <- ymd(FF_monthly$date)

FF_daily$u_sq <- log(1+FF_daily$`mkt-rf` / 100)^2
FF_monthly$u_sq <- log(1+FF_monthly$`mkt-rf` / 100)^2

################################################################################
#******************************* Monthly Level *******************************#

# Define Used Strategies
names <- c("var_managed","ARIMA_var_managed", "ESA_var_managed", "GARCH_var_managed")
names_clean <- c("Var", "ARIMA", "ESA", "GARCH")

# Calculate Monthly Variance and Volatility
trading_days <- 22
trading_months <- 12
trading_year <- trading_days * trading_months

var_m <- FF_daily %>%
  mutate(month = month(date), year = year(date)) %>%
  group_by(year, month) %>%
  summarise(var = var(`mkt-rf`) * trading_days * 
              (length(`mkt-rf`)-1) / length(`mkt-rf`))

var_m <- var_m %>% 
  as.data.frame() %>% 
  mutate(date = FF_monthly$date, vol = sqrt(var)) %>% 
  select(date, var, vol)

# Conduct Tests for Stationarity and Calculate ARIMA Variance 
variance_ts_m <- xts(var_m$var, order.by = var_m$date)

kpss.test(var_m$var)
kpss.test(diff(var_m$var))

# Only Change Display of Graphs if Window is Large Enough, Otherwise Margin Error
# par(mfrow=c(2,1))
Acf(diff(variance_ts_m), lag = 24, main = "Autocorrelation Function")
Pacf(diff(variance_ts_m), lag = 24, main = "Partial Autocorrelation Function")

ARIMA_model_m <- auto.arima(variance_ts_m,stepwise = FALSE, approximation = FALSE,
                            max.p = 10, max.q = 10, ic = "aic", test = "kpss")

var_m$ARIMA_var <- c(fitted(ARIMA_model_m)[-1], forecast(ARIMA_model_m, h = 1)$mean)
var_m <- var_m %>% mutate(ARIMA_vol = sqrt(ARIMA_var))

# Optimize Lambda for EWMA Model
EWMA_function_m <- function(lambda) {
  EWMA_var <- c(1:n_months)
  EWMA_var[1] <- FF_monthly$u_sq[1]
  for(i in 2:n_months) {
    EWMA_var[i] <- lambda * EWMA_var[i-1] + (1 - lambda) * FF_monthly$u_sq[i]
  }
  EWMA_likelihood <- c(1:(n_months-1))
  for(i in 1:(n_months-1)) {
    EWMA_likelihood[i] <- -log(EWMA_var[i]) - FF_monthly$u_sq[i+1] / EWMA_var[i]
  }
  return (sum(EWMA_likelihood))
}

EWMA_max_m <- optimize(EWMA_function_m, interval = c(0, 1), maximum = TRUE, 
                       tol = 0.000000000000001)
lambda_m <- EWMA_max_m$maximum

# Calculate ESA Variance and Volatility
var_m$ESA_var <- c(1:n_months)
var_m$ESA_var[1] <- var_m$var[1]

for (i in 2:n_months) {
  var_m$ESA_var[i] <- lambda_m * var_m$ESA_var[i - 1] + (1 - lambda_m) * var_m$var[i]
}

var_m <- var_m %>% mutate(ESA_vol = sqrt(var_m$ESA_var))

# Optimize GARCH Parameter
GARCH_function_m <- function(alpha, beta) {
  omega <- max(0,mean(FF_monthly$u_sq)*(1-alpha-beta))
  GARCH_var <- c(1:n_months)
  GARCH_var[1] <- FF_monthly$u_sq[1]
  for(i in 2:n_months) {
    GARCH_var[i] <- omega + beta*GARCH_var[i-1] + alpha*FF_monthly$u_sq[i]
  }
  GARCH_likelihood <- c(1:(n_months - 1))
  for(i in 1:(n_months-1)) {
    GARCH_likelihood[i] <- -log(GARCH_var[i])-FF_monthly$u_sq[i+1]/GARCH_var[i]
  }
  return (sum(GARCH_likelihood))
}

GARCH_max_m <- optimx(c(0.1, 0.9), function(x) GARCH_function_m(x[1], x[2]), 
                      method = "Nelder-Mead", control = list(maximize = TRUE))

alpha_m <- GARCH_max_m$p1
beta_m <- GARCH_max_m$p2
omega_m <- max(0, mean(FF_monthly$u_sq) * (1 - alpha_m-beta_m))

# Calculate GARCH Variance and Volatility
var_m$GARCH_var <- c(1:n_months)
var_m$GARCH_var[1] <- FF_monthly$u_sq[1]
for (i in 2:n_months) {
  var_m$GARCH_var[i] <- (omega_m +
    alpha_m*FF_monthly$u_sq[i] + beta_m * var_m$GARCH_var[i-1])
}
var_m <- var_m %>% mutate(GARCH_var = GARCH_var * 10000, GARCH_vol = sqrt(GARCH_var))

# Set Up Parameter for Following GGPlots
black_col <- brewer.pal(9, name = "Greys")[4:8]
blue_col <- brewer.pal(9, name = "Blues")
green_col <- brewer.pal(9, name = "Greens")[5]
grey_col <- brewer.pal(9, name = "Greys")[5]

line_size <- 0.7
loadfonts(device = "win")
windowsFonts(Times = windowsFont("TT Times New Roman"))

# Plot Variance Estimates by 4 Measures
ggplot(var_m, aes(x = date)) +
  geom_line(aes(y = vol * sqrt(trading_months), color = "Realized Variance"), size = line_size) +
  geom_line(aes(y = ARIMA_vol * sqrt(trading_months), color = "ARIMA"), size = line_size) +
  geom_line(aes(y = ESA_vol * sqrt(trading_months), color = "ESA"), size = line_size) +
  geom_line(aes(y = GARCH_vol * sqrt(trading_months), color = "GARCH"), size = line_size) +
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
  ggtitle("Annualized Volatility") +
  xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# Analyze Predictive Power of Variance Estimates with Regressions
var_names <- c("var", "ARIMA_var", "ESA_var", "GARCH_var")
var_var_reg_m <- vector(mode = "list", length = length(var_names))
ret_var_reg_m <- vector(mode = "list", length = length(var_names))
retvar_var_reg_m <- vector(mode = "list", length = length(var_names))

var_var_reg_se_m <- vector(mode = "list", length = length(var_names))
ret_var_reg_se_m <- vector(mode = "list", length = length(var_names))
retvar_var_reg_se_m <- vector(mode = "list", length = length(var_names))

for (i in 1:length(var_names)) {
  var_var_reg_m[[i]] <- lm(var_m$var[-1] ~ var_m[-n_months, var_names[i]])
  ret_var_reg_m[[i]] <- lm(FF_monthly$`mkt-rf`[-1] ~ var_m[-n_months, var_names[i]])
  retvar_var_reg_m[[i]] <- lm((FF_monthly$`mkt-rf`[-1] / var_m$var[-1]) ~
                                var_m[-n_months, var_names[i]])
  names(var_var_reg_m[[i]]$coefficients) <- 
    c(names(var_var_reg_m[[i]]$coefficients[1]), var_names[i])
  names(ret_var_reg_m[[i]]$coefficients) <- 
    c(names(var_var_reg_m[[i]]$coefficients[1]), var_names[i])
  names(retvar_var_reg_m[[i]]$coefficients) <- 
    c(names(var_var_reg_m[[i]]$coefficients[1]), var_names[i])
}

for (i in 1:length(names)) {
  var_var_reg_se_m[[i]] <- 
    var_var_reg_m[[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
  ret_var_reg_se_m[[i]] <- 
    ret_var_reg_m[[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
  retvar_var_reg_se_m[[i]] <- 
    retvar_var_reg_m[[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
}

# Set Up Function to Calculate Significance Stars
ff3_alpha_stars <- function(model) {
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

# Output Relevant Statistics
stargazer(var_var_reg_m[[1]], var_var_reg_m[[2]], 
          var_var_reg_m[[3]], var_var_reg_m[[4]],
          ret_var_reg_m[[1]], ret_var_reg_m[[2]], 
          ret_var_reg_m[[3]], ret_var_reg_m[[4]],
          retvar_var_reg_m[[1]], retvar_var_reg_m[[2]], 
          retvar_var_reg_m[[3]], retvar_var_reg_m[[4]],
          se = list(var_var_reg_se_m[[1]], var_var_reg_se_m[[2]], 
                    var_var_reg_se_m[[3]], var_var_reg_se_m[[4]], 
                    ret_var_reg_se_m[[1]], ret_var_reg_se_m[[2]], 
                    ret_var_reg_se_m[[3]], ret_var_reg_se_m[[4]], 
                    retvar_var_reg_se_m[[1]], retvar_var_reg_se_m[[2]], 
                    retvar_var_reg_se_m[[3]], retvar_var_reg_se_m[[4]]),
          type = "text", 
          table.layout = "-d#-t-s-",
          dep.var.labels = c("Variance Next Month","Return Next Month",
                             "Return per Unit of Variance"), 
          covariate.labels = c("Realized Variance", "ARIMA", "ESA", "GARCH"),
          keep.stat = c("n", "rsq"),
          out = "table.htm", digits = 3, digits.extra = 0)

# Calculate Parameter c with Midnight Formula
c_m <- data.frame(matrix(ncol = length(names)))
colnames(c_m) <- names

a_qe_m <- c(1:length(names))
b_qe_m <- c(1:length(names))
c_qe_m <- c(1:length(names))

# Usage of Population Variance, since the Prefactor Cancels out on Aggregate
for (i in 1:length(names)) {
  a_qe_m[i] <- var(FF_monthly$`mkt-rf`[-1] / var_m[-n_months,var_names[i]])
  b_qe_m[i] <- 2 * cov(FF_monthly$`mkt-rf`[-1] / var_m[-n_months,var_names[i]], 
                       FF_monthly$rf[-1])
  c_qe_m[i] <- var(FF_monthly$rf[-1]) - var(FF_monthly$mkt[-1])
  
  c_m[i] <- 1 / (2 * a_qe_m[i]) * 
    (-b_qe_m[i] + sqrt((b_qe_m[i])^2 - 4 * a_qe_m[i] * c_qe_m[i]))
}

# Calculate Weights, Absolute Weight Deviation and Returns
w_abs_m <- data.frame(matrix(ncol = length(names) + 1, nrow = n_months - 1))
colnames(w_abs_m) <- c("date", names)
w_abs_m <- w_abs_m %>% mutate(date = FF_monthly$date[-1])

weights_m <- data.frame(matrix(ncol = length(names) + 1, nrow = n_months - 1))
colnames(weights_m) <- c("date", names)
weights_m <- weights_m %>% mutate(date = FF_monthly$date[-1])

returns_m <- data.frame(matrix(ncol = length(names) + 4, nrow = n_months - 1))
colnames(returns_m) <- c("date", "mkt", "rf", "mkt-rf", names)
returns_m <- returns_m %>% mutate(date = FF_monthly$date[-1], 
                                  mkt = FF_monthly$mkt[-1], 
                                  rf = FF_monthly$rf[-1],
                                  "mkt-rf" = FF_monthly$`mkt-rf`[-1])

for (i in 1:length(names)) {
  weights_m[,names[i]] <- as.numeric(c_m[i]) / var_m[-n_months,var_names[i]]
  returns_m[,names[i]] <- weights_m[,names[i]] * returns_m$`mkt-rf` + returns_m$rf
  w_abs_m[1,names[i]] <- 1
  w_abs_m[-1,names[i]] <- weights_m[,names[i]] %>% diff() %>% abs()
}

# Run Regressions to Determine Alpha, Beta, etc.
reg_mkt_m <- vector(mode = "list", length = length(names))
reg_mkt_se_m <- vector(mode = "list", length = length(names))
reg_FF3_m <- vector(mode = "list", length = length(names))
reg_FF3_se_m <- vector(mode = "list", length = length(names))

b_m <- trading_months * returns_m$`mkt-rf`
b1_m <- trading_months * (FF_monthly$SMB[-1])
b2_m <- trading_months * (FF_monthly$HML[-1])

a_m <- as.data.frame(matrix(nrow = n_months - 1, ncol = length(names)))
colnames(a_m) <- names

for (i in 1:length(names)) {
  a_m[,names[i]] <- trading_months * (returns_m[, names[i]] - returns_m$rf)
  reg_mkt_m[[i]] <- lm(a_m[,names[i]] ~ b_m)
  reg_FF3_m[[i]] <- lm(a_m[,names[i]] ~ b_m + b1_m + b2_m)
  reg_mkt_se_m[[i]] <- reg_mkt_m[[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
  reg_FF3_se_m[[i]] <- reg_FF3_m[[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
}

# Create Ouput Table
output_names <- c("alpha_mkt", "beta_mkt", "RMSE", "SR", "AR", "alpha_FF3",
                  "alpha_FF3_se")
reg_output_m <- data.frame(matrix(ncol = length(names), nrow = length(output_names)))
colnames(reg_output_m) <- names
rownames(reg_output_m) <- output_names

for (i in 1:length(names)) {
  reg_output_m["alpha_mkt", i] <- reg_mkt_m[[i]]$coefficients[1]
  reg_output_m["beta_mkt", i] <- reg_mkt_m[[i]]$coefficients[2]
  reg_output_m["RMSE", i] <- sqrt(mean(residuals(reg_mkt_m[[i]])^2))
  reg_output_m["SR", i] <- sqrt(trading_months) * 
    (mean(returns_m[,names[i]] - returns_m$rf)) / sd(returns_m[,names[i]])
  reg_output_m["AR", i] <- sqrt(trading_months) * 
    reg_output_m["alpha_mkt", i] / reg_output_m["RMSE", i]
  reg_output_m["alpha_FF3", i] <- reg_FF3_m[[i]]$coefficients[1]
  reg_output_m["alpha_FF3_se", i] <- coeftest(reg_FF3_m[[i]], 
                                              vcovHC(reg_FF3_m[[i]], 
                                                     type = "HC"))[1,2]
}

alpha_stars_m <- vector(length = length(names))
for (i in 1:length(names)) {
  alpha_stars_m[i] <- ff3_alpha_stars(reg_FF3_m[[i]])
}

stargazer(reg_mkt_m[[1]], reg_mkt_m[[2]], reg_mkt_m[[3]], reg_mkt_m[[4]],
          se = list(reg_mkt_se_m[[1]], reg_mkt_se_m[[2]],
                    reg_mkt_se_m[[3]], reg_mkt_se_m[[4]]),
          type = "text", 
          keep.stat = c("n", "rsq", "ser"), df = FALSE, 
          column.labels = names_clean,
          table.layout = "-#c-t-s-a-",
          covariate.labels = c("Mkt-RF", "Alpha (&#945;)"),
          add.lines = list(c("Sharpe Ratio", 
                             sprintf("%.2f", as.numeric(reg_output_m["SR",]))),
                           c("Appr Ratio", 
                             sprintf("%.2f", as.numeric(reg_output_m["AR",]))),
                           c("Alpha FF3", 
                             paste(sprintf("%.2f", as.numeric(reg_output_m["alpha_FF3",])),
                                   alpha_stars_m, sep = "")),
                           c("", 
                             paste("(", 
                                   sprintf("%.2f", as.numeric(reg_output_m["alpha_FF3_se",])),
                                   ")", sep = ""))),
          digits = 2, out = "table.htm")

stargazer(reg_FF3_m[[1]], reg_FF3_m[[2]], reg_FF3_m[[3]], reg_FF3_m[[4]],
          se = list(reg_FF3_se_m[[1]], reg_FF3_se_m[[2]],
                    reg_FF3_se_m[[3]], reg_FF3_se_m[[4]]),
          type = "text", 
          keep.stat = c("n", "rsq"), df = FALSE, 
          column.labels = names_clean,
          table.layout = "-#c-t-s-",
          covariate.labels = c("Mkt-RF", "SMB", "HML", "Alpha (&#945;)"),
          digits = 2, out = "table.htm")

# Calculate and Plot Cumulative Return
cum_ret_m <- data.frame(matrix(ncol = length(names) + 2, nrow = n_months))
colnames(cum_ret_m) <- c("date", "mkt", names)
cum_ret_m <- cum_ret_m %>% mutate(date = FF_monthly$date)

cum_ret_m[1,-1] <- 1

for (i in 2:n_months) {
  cum_ret_m$mkt[i] <- cum_ret_m$mkt[i-1] * (1 + (returns_m$mkt[i-1] / 100))
  for (j in 1:length(names)) {
    cum_ret_m[i,names[j]] <- cum_ret_m[i-1, names[j]] * 
      (1 + (returns_m[i-1, names[j]] / 100))
    } 
}

cum_ret_m[n_months,]

ggplot(cum_ret_m, aes(x = date)) +
  geom_line(aes(y=mkt, color = "Buy and Hold"), size = line_size) +
  geom_line(aes(y=var_managed, color = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, color = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, color = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, color = "GARCH"), size = line_size) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks('log10', function(x) 10^x),
                     minor_breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
                                      1,2,3,4,5,6,7,8,9,
                                      10,20,30,40,50,60,70,80,90,
                                      100,200,300,400,500,600,700,800,900,
                                      1000,2000,3000,4000,5000,6000,7000,8000,9000,
                                      10000,20000,30000,40000,50000,
                                      60000,70000,80000,90000,100000),
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
        axis.text.y = element_text(size = 11)) +
  ggtitle("Cumulative Performance") + 
  xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = "#000000",
                                "Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Buy and Hold", "Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# Calculate Rolling Average Returns
returns_rolling_m <- data.frame(matrix(nrow = n_months - 1, ncol = length(names) + 2))
colnames(returns_rolling_m) <- c("date", "mkt", names)
returns_rolling_m <- returns_rolling_m %>% mutate(date = returns_m$date)
rolling_period <- 12

for (i in 1:(n_months -1 )) {
  returns_rolling_m[i,"mkt"] <- mean(returns_m[(i-min(i,rolling_period-1)):i,"mkt"])
  for (j in 1:length(names)) {
    returns_rolling_m[i,names[j]] <- mean(returns_m[(i-min(i,rolling_period-1)):i,names[j]])
  }
}

ggplot(returns_rolling_m, aes(x=date)) +
  geom_line(aes(y = ARIMA_var_managed, color = "ARIMA"), size = line_size) +
  geom_line(aes(y = ESA_var_managed, color = "ESA"), size = line_size) +
  geom_line(aes(y = GARCH_var_managed, color = "GARCH"), size = line_size) +
  geom_line(aes(y = mkt, color = "Buy and Hold"), size = line_size) +
  geom_line(aes(y = var_managed, color = "Realized Variance"), size = line_size) +
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
  ggtitle("Rolling One-Year Return") +
  xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = "#000000",
                                "Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Buy and Hold", "Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# Calculate Drawdowns
drawdown_m <- data.frame(matrix(nrow = n_months - 1, ncol = length(names) + 2))
colnames(drawdown_m) <- c("date", "mkt", names)
drawdown_m <- drawdown_m %>% mutate(date = returns_m$date)
drawdown_m[1,-1] <- 0

for (i in 2:(n_months - 1)) {
  drawdown_m$mkt[i] <- min(0, (cum_ret_m$mkt[i]-max(cum_ret_m$mkt[1:i]))/
                             max(cum_ret_m$mkt[1:i]))
  for (j in 1:length(names)) {
    drawdown_m[i,names[j]] <- min(0, (cum_ret_m[i,names[j]]-max(cum_ret_m[1:i,names[j]]))/
                                    max(cum_ret_m[1:i,names[j]]))
  }
}

ggplot(drawdown_m, aes(x = date)) +
  geom_line(aes(y=mkt, color = "Buy and Hold"), size = line_size) +
  geom_line(aes(y=var_managed, color = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, color = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, color = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, color = "GARCH"), size = line_size) +
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
  ggtitle("Drawdown") + 
  xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = "#000000",
                                "Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Buy and Hold", "Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# Analyze Distribution of Returns and Weigths
stat_ret_names <- c("Mean", "SD", "Skewness", "Kurtosis", "Min", "P25", "P50",
                    "P75", "Max")
stat_ret_m <- data.frame(matrix(ncol = length(names), 
                                nrow = length(stat_ret_names),
                                dimnames = list(stat_ret_names, names_clean)))

stat_weight_names <- c("Mean", "SD", "Skewness", "Kurtosis","P50","P75", "P90", 
                       "P99", "Max")
stat_weight_m <- data.frame(matrix(ncol = length(names),
                                   nrow = length(stat_weight_names),
                                   dimnames = list(stat_weight_names, names_clean)))

stat_ret_m[1,] <- apply(returns_m[,names], 2, mean)
stat_ret_m[2,] <- apply(returns_m[,names], 2, sd)
stat_ret_m[3,] <- apply(returns_m[,names], 2, skewness)
stat_ret_m[4,] <- apply(returns_m[,names], 2, kurtosis)
stat_ret_m[5:9,] <- as.data.frame(apply(returns_m[,names], 2, quantile))

stat_weight_m[1,] <- apply(weights_m[,names], 2, mean)
stat_weight_m[2,] <- apply(weights_m[,names], 2, sd)
stat_weight_m[3,] <- apply(weights_m[,names], 2, skewness)
stat_weight_m[4,] <- apply(weights_m[,names], 2, kurtosis)
stat_weight_m[5:9,] <- as.data.frame(apply(weights_m[,names], 2, quantile, 
                                           probs = c(0.5, 0.75, 0.9, 0.99, 1)))

stargazer(stat_ret_m, type = "text",summary = FALSE, 
          title = "Panel A: Return Distribution",
          out = "table.htm", digits = 2)

stargazer(stat_weight_m, type = "text",summary = FALSE, 
          title = "Panel B: Weight Distribution",
          out = "table.htm", digits = 2)

# Plot Density Function of Returns
ggplot(returns_m) +
  stat_density(aes(x=mkt, color = "Buy and Hold"), 
               adjust = 1, size = line_size, geom = "line") +
  stat_density(aes(x=var_managed, color = "Realized Variance"), 
               adjust = 1, size = line_size, geom = "line") +
  stat_density(aes(x=ARIMA_var_managed, color = "ARIMA"), 
               adjust = 1, size = line_size, geom = "line") +
  stat_density(aes(x=ESA_var_managed, color = "ESA"), 
               adjust = 1, size = line_size, geom = "line") +
  stat_density(aes(x=GARCH_var_managed, color = "GARCH"), 
               adjust = 1, size = line_size, geom = "line") +
  xlim(-40,40) +
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
  ylim(0,0.15) +
  ggtitle("Density Function of Returns") +
  ylab("") + xlab("") +
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = "#000000",
                                "Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Buy and Hold", "Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

ggplot(weights_m) +
  stat_density(aes(x=var_managed, color = "Realized Variance"), 
               adjust = 1, size = line_size, geom = "line") +
  stat_density(aes(x=ARIMA_var_managed, color = "ARIMA"), 
               adjust = 1, size = line_size, geom = "line") +
  stat_density(aes(x=ESA_var_managed, color = "ESA"), 
               adjust = 1, size = line_size, geom = "line") +
  stat_density(aes(x=GARCH_var_managed, color = "GARCH"), 
               adjust = 1, size = line_size, geom = "line") +
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
  ylim(0,1) +
  ggtitle("Density Function of Weights") +
  ylab("") + xlab("") +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Buy and Hold", "Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

# Compute Impact of Trading Costs and Breakeven Costs
breakeven_function_m <- function(cost, strategy) {
  returns_be_m <- c(1:(n_months-1))
  returns_be_m <- returns_m[,names[strategy]] - w_abs_m[,names[strategy]] * cost
  a <- trading_months * (returns_be_m - returns_m$rf)
  return(lm(a ~ b_m)$coefficients[1])
}

cost_m <- data.frame(matrix(nrow = length(names), ncol = 6))
colnames(cost_m) <- c("Alpha", "|Delta w|", "1bps", "10bps", "14bps", "Break-Even")
rownames(cost_m) <- names_clean

for (i in 1:length(names)) {
  cost_m[i,1] <- breakeven_function_m(0, i)
  cost_m[i,2] <- mean(w_abs_m[,names[i]])
  cost_m[i,3] <- breakeven_function_m(0.01, i)
  cost_m[i,4] <- breakeven_function_m(0.1, i)
  cost_m[i,5] <- breakeven_function_m(0.14, i)
  cost_m[i,6] <- paste(
    sprintf("%.0f", 
            uniroot(breakeven_function_m, strategy = i, lower = 0, upper = 100)$root *100),
    "bps", sep = "")
}

stargazer(cost_m, type = "text", summary = FALSE, out = "table.htm", digits = 2)

# Compute Impact of Leverage Constraints
leverage_function_m <- function(leverage, strategy) {
  returns_le_m <- c(1:(n_months-1))
  for (i in 1:(n_months-1)) {
    returns_le_m[i] <- min(weights_m[i, names[strategy]], 1 + leverage) *
      returns_m$`mkt-rf`[i] + returns_m$rf[i]
  }
  a <- trading_months * (returns_le_m - returns_m$rf)
  return(lm(a ~ b_m)$coefficients[1])
}

leverage_m <- data.frame(matrix(nrow = length(names), ncol = 4))
colnames(leverage_m) <- c("No Constraint", "100% Leverage", "50% Leverage", "0% Leverage")
rownames(leverage_m) <- names_clean

for (i in 1:length(names)) {
  leverage_m[i,1] <- breakeven_function_m(0, i)
  leverage_m[i,2] <- leverage_function_m(1, i)
  leverage_m[i,3] <- leverage_function_m(0.5, i)
  leverage_m[i,4] <- leverage_function_m(0, i)
}

stargazer(leverage_m, type = "text", summary = FALSE, out = "table.htm", digits = 2)

# Replicate Regressions Controlling for Recessions
reg_rec_m <- vector(mode = "list", length = length(names))
reg_rec_se_m <- vector(mode = "list", length = length(names))

rec_m <- recession$indicator[-1]
b_rec_m <- trading_months * returns_m$`mkt-rf` * recession$indicator[-1]

for (i in 1:length(names)) {
  reg_rec_m[[i]] <- lm(a_m[[i]] ~ b_m + b_rec_m + rec_m)
  reg_rec_se_m[[i]] <- reg_rec_m[[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
}

stargazer(reg_rec_m[[1]], reg_rec_m[[2]],
          reg_rec_m[[3]], reg_rec_m[[4]],
          se = list(reg_rec_se_m[[1]], reg_rec_se_m[[2]],
                    reg_rec_se_m[[3]], reg_rec_se_m[[4]]),
          column.labels = names_clean,
          table.layout = "-#c-t-s-a-",
          covariate.labels = c("Mkt-RF", "Mkt-RF * Recession",
                               "Recession", "Alpha (&#945;)"),
            type = "text", out = "table.htm", 
          keep.stat = c("n", "rsq"), digits = 2)

# Determine Alpha of Alternative Strategies When Controlling For Moreira & Muir
reg_strat_m <- vector(mode = "list", length = length(names) - 1)
reg_strat_se_m <- vector(mode = "list", length = length(names) - 1)

for (i in 1:(length(names) - 1)) {
  reg_strat_m[[i]] <- lm(a_m[,names[i+1]] ~ a_m[,names[i]] + b_m)
  reg_strat_se_m[[i]] <- reg_strat_m[[i]] %>% 
    vcovHC(type = "HC") %>% sqrt() %>% diag()
}


stargazer(reg_strat_m[[1]], reg_strat_m[[2]], reg_strat_m[[3]],
          se = list(reg_strat_se_m[[1]], reg_strat_se_m[[2]],
                    reg_strat_se_m[[3]]),
          type = "text", 
          column.labels = names_clean[-1],
          table.layout = "-#c-t-s-a-",
          dep.var.labels = names_clean[-1], 
          covariate.labels = c("Mkt-RF", "Var", "Alpha (&#945;)"), 
          out = "table.htm",
          keep.stat = c("n", "rsq"), digits = 2)

# Repeat Regression Controlling For Cost
alternative_regression_function <- function(cost, strategy) {
  returns_var_m <- c(1:(n_months-1))
  returns_var_m <- returns_m[,names[1]] - w_abs_m[,names[1]] * cost
  returns_alt_m <- c(1:(n_months-1))
  returns_alt_m <- returns_m[,names[strategy]] - w_abs_m[,names[strategy]] * cost
  a <- trading_months * (returns_alt_m - returns_m$rf)
  m <- trading_months * (returns_var_m - returns_m$rf)
  return(lm(a ~ m + b_m))
}

cost_vector <- c(0, 0.01, 0.1, 0.14, 0.5)
strat_cost_m <- data.frame(matrix(nrow = length(cost_vector) * 2, 
                                  ncol = length(names_clean[-1])))
colnames(strat_cost_m) <- names_clean[-1]

q <- 1
for (i in 1:(length(cost_vector) * 2)) {
  if (i %% 2 == 1) {
    rownames(strat_cost_m)[i] <- paste(cost_vector[q] * 100, "bps")
  }
  else {
    rownames(strat_cost_m)[i] <- paste("se", q, sep = "")
    q <- q + 1
  }
}

q <- 1
for (i in 1:(length(cost_vector) * 2)) {
  for (j in 1:(length(names[-1]))) {
    reg <- alternative_regression_function(cost_vector[q],j+1)
    if (i %% 2 == 1){
      strat_cost_m[i,j] <- 
        coefficients(reg)[1] %>% round(2) %>% format(nsmall = 2)
    }
    else {
      strat_cost_m[i,j] <- 
        paste(" (",
              coeftest(reg, vcovHC(reg, type = "HC"))[1,2] %>% 
                round(2) %>% format(nsmall = 2),
              ")", 
              sep = "")
    }
  }  
  if(i %% 2 == 0) q <- q + 1
}

stargazer(strat_cost_m, type = "text", summary = FALSE, out = "table.htm")

# Calculate Final Return, Considering Costs
cum_return_function_m <- function(cost, strategy) {
  cum_returns_m <- c(1:(n_months))
  cum_returns_m[1] <- 1
  for (i in 2:n_months) {
    ret <- returns_m[i-1, names[strategy]] - w_abs_m[i-1, names[strategy]] * cost
    cum_returns_m[i] <- cum_returns_m[i-1] * (1 + (ret / 100))
  }
  return(cum_returns_m[n_months])
}

round(cum_return_function_m(0.1, 1) / cum_return_function_m(0, 1) - 1, 2)
round(cum_return_function_m(0.1, 2) / cum_return_function_m(0, 2) - 1, 2)
round(cum_return_function_m(0.1, 3) / cum_return_function_m(0, 3) - 1, 2)
round(cum_return_function_m(0.1, 4) / cum_return_function_m(0, 4) - 1, 2)

# Analysis of Bad Periods With Delta w
bd_times <- returns_m
bd_times <- bd_times %>% mutate(weight_Var = weights_m$var_managed,
                                weight_ARIMA = weights_m$ARIMA_var_managed,
                                weight_ESA = weights_m$ESA_var_managed, 
                                weight_GARCH = weights_m$GARCH_var_managed)

extreme_periods <- list()
extreme_periods[[1]] <- bd_times %>% filter(mkt < 0 & weight_Var > 0)
extreme_periods[[2]] <- bd_times %>% filter(mkt < -5 & weight_Var > 1)
extreme_periods[[3]] <- bd_times %>% filter(var_managed < -10 & weight_Var > 1)
extreme_periods[[4]] <- bd_times %>% filter(var_managed < -15 & weight_Var > 1)
bad_times_names2 <- c("w > 1", "Mkt < 0% w > 1", "Mkt < -5% w > 1",
                      "Var < -10% w > 1", "Var < -15% w > 1")

delta_m <- data.frame(matrix(nrow = length(names),
                             ncol = length(bad_times_names2)))
rownames(delta_m) <- names_clean
colnames(delta_m) <- bad_times_names2

weight_names <- c(1:length(names))

for (i in 1:length(names)) {
  weight_names[i] <- paste("weight", names_clean[i], sep = "_")
  delta_m[i, 1] <- mean(filter(bd_times, weight_Var > 1)[,weight_names[i]])
  for (j in 1:length(extreme_periods)) {
    delta_m[i, j + 1] <- mean(extreme_periods[[j]][,weight_names[i]])
  }
}

rel_delta_m <- data.frame(matrix(nrow = length(names),
                                 ncol = length(bad_times_names2)))
rownames(rel_delta_m) <- names_clean
colnames(rel_delta_m) <- bad_times_names2

for (i in 1:length(names)) {
  rel_delta_m[i, 1] <-
    mean(filter(bd_times, weight_Var > 1)[,weight_names[i]] - 
           filter(bd_times, weight_Var > 1)[,weight_names[1]]) /
    mean(filter(bd_times, weight_Var > 1)[,weight_names[1]])
  for (j in 1:length(extreme_periods)) {
    rel_delta_m[i, j + 1] <-
      mean(extreme_periods[[j]][,weight_names[i]] - 
             extreme_periods[[j]][,weight_names[1]]) /
      mean(extreme_periods[[j]][,weight_names[1]])
  }
}

rel_delta_m <- rel_delta_m %>% round(digits = 2) %>% format(nsmall = 2)
stargazer(delta_m, summary = FALSE, type = "text", 
          title = "Panel A: Average Absolute Weights",
          out = "table.htm", digits = 2)
stargazer(rel_delta_m, summary = FALSE, type = "text", 
          title = "Panel B: Relative Deviation from Realized Variance Strategy",
          out = "table.htm", digits = 2)

# Analysis of Bad Times With Correlation 

cor_names_m <- c("Observations", "Mkt-Var ", "Mkt-ARIMA ", "Mkt-ESA ", 
                 "Mkt-GARCH ", "Var-ARIMA ", "Var-ESA ", "Var-GARCH ")
bad_times_names <- c("Whole Sample", "Mkt < 0% w > 1", "Mkt < -5% w > 1",
                     "Var < -10% w > 1", "Var < -15% w > 1")

cor_m <- data.frame(matrix(nrow = length(cor_names_m),
                           ncol = length(bad_times_names)))

rownames(cor_m) <- cor_names_m
colnames(cor_m) <- bad_times_names

# Compute Respective Statistics for Whole Sample
cor_m[1,1] <- nrow(returns_m)

for (i in 1:(length(cor_names_m) - 1)) {
  if (i <= length(names)) {
    cor_m[i + 1,1] <- cor(returns_m$mkt, returns_m[,names[i]])
  } 
  else {
    cor_m[i + 1,1] <- cor(returns_m$var_managed, returns_m[,names[i-length(names)+1]])
  }
}

# Compute Respective Statistics for Extreme Periods
for (i in 1:length(extreme_periods)) {
  cor_m[1, i + 1] <- nrow(extreme_periods[[i]])
  for (j in 1:(length(cor_names_m) - 1)) {
    if (j <= length(names)) {
      cor_m[j + 1, i + 1] <- cor(extreme_periods[[i]]$mkt, 
                                 extreme_periods[[i]][,names[j]])
    } 
    else {
      cor_m[j + 1,i + 1] <- cor(extreme_periods[[i]]$var_managed, 
                            extreme_periods[[i]][,names[j-length(names)+1]])
    }
  }
}

stargazer(cor_m, summary = FALSE, type = "text", out = "table.htm", digits = 2)

################################################################################
#****************************** Three Subperiods ******************************#

# Create Subperiods and Repeat Steps from above to Derive Regression Table
sub_dates <- data.frame(matrix(nrow = 3, ncol = 2, 
                               dimnames = list(c("Period 1", "Period 2", "Period 3"),
                                               c("First Day", "Last Day"))))
output_names1 <- c("alpha_mkt", "beta_mkt", "RMSE", "SD", "SR", "AR", "alpha_FF3",
                  "alpha_FF3_se")

sub_dates["Period 1",] <- c(first_day, 19570730)
sub_dates["Period 2",] <- c(19570801, 19880730)
sub_dates["Period 3",] <- c(19880801, last_day)

reg_mkt_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))
reg_mkt_se_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))
reg_FF3_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))

a_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))
b_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))
b1_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))
b2_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))

reg_output_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))
alpha_stars_m_sub <- vector(mode = "list", length = as.numeric(count(sub_dates)))

for (period in 1:as.numeric(count(sub_dates))){
  first_day_sub <- ymd(sub_dates[period, 1])
  last_day_sub <- ymd(sub_dates[period, 2])
  
  FF_monthly_sub <- FF_monthly %>%
    subset(subset = date < last_day_sub & date >= first_day_sub)
  
  var_m_sub <- var_m %>%
    subset(subset = date < last_day_sub & date >= first_day_sub)
  
  returns_m_sub <- returns_m %>% 
    subset(subset = date < last_day_sub & date > first_day_sub)
  
  weights_m_sub <- weights_m %>% 
    subset(subset = date < last_day_sub & date > first_day_sub)
  
  n_months_sub <- var_m_sub %>% count() %>% as.numeric()
  
  c_m_sub <- data.frame(matrix(ncol = length(names)))
  colnames(c_m_sub) <- names
  
  a_qe_m_sub <- c(1:length(names))
  b_qe_m_sub <- c(1:length(names))
  c_qe_m_sub <- c(1:length(names))
  
  # Usage of Population Variance, since the prefactor cancels out on aggregate
  for (i in 1:length(names)) {
    a_qe_m_sub[i] <- var(FF_monthly_sub$`mkt-rf`[-1] / 
                           var_m_sub[-n_months_sub,var_names[i]])
    b_qe_m_sub[i] <- 2 * cov(FF_monthly_sub$`mkt-rf`[-1] / 
                               var_m_sub[-n_months_sub,var_names[i]],
                             FF_monthly_sub$rf[-1])
    c_qe_m_sub[i] <- var(FF_monthly_sub$rf[-1]) - var(FF_monthly_sub$mkt[-1])
    
    c_m_sub[i] <- 1 / (2 * a_qe_m_sub[i]) * 
      (-b_qe_m_sub[i] + sqrt((b_qe_m_sub[i])^2 - 4 * a_qe_m_sub[i] * c_qe_m_sub[i]))
  }
  
  for (i in 1:length(names)) {
    weights_m_sub[,names[i]] <- 
      as.numeric(c_m_sub[i]) / var_m_sub[-n_months_sub,var_names[i]]
    returns_m_sub[,names[i]] <- 
      weights_m_sub[,names[i]] * returns_m_sub$`mkt-rf` + returns_m_sub$rf
  }
  
  # Run Regressions to Determine Alpha, Beta, etc.
  reg_mkt_m_sub[[period]] <- vector(mode = "list", length = length(names))
  reg_mkt_se_m_sub[[period]] <- vector(mode = "list", length = length(names))
  reg_FF3_m_sub[[period]] <- vector(mode = "list", length = length(names))
  
  b_m_sub[[period]] <- trading_months * FF_monthly_sub$`mkt-rf`[-1]
  b1_m_sub[[period]] <- trading_months * FF_monthly_sub$SMB[-1]
  b2_m_sub[[period]] <- trading_months * FF_monthly_sub$HML[-1]
  
  a_m_sub[[period]] <- as.data.frame(matrix(nrow = n_months_sub - 1, ncol = length(names)))
  colnames(a_m_sub[[period]]) <- names
  
  reg_output_m_sub[[period]] <- data.frame(matrix(ncol = length(names),
                                                  nrow = length(output_names1),
                                                  dimnames = list(output_names1,
                                                                  names)))
  
  for (i in 1:length(names)) {
    a_m_sub[[period]][,names[i]] <- 
      trading_months * (returns_m_sub[, names[i]] - returns_m_sub$rf)
    
    reg_output_m_sub[[period]]["SD", i] <- 
      sqrt(trading_months) * sd(returns_m_sub[,names[i]])
    reg_output_m_sub[[period]]["SR", i] <- sqrt(trading_months) * 
      (mean(returns_m_sub[,names[i]] - returns_m_sub$rf)) / 
      sd(returns_m_sub[,names[i]])
  }
}  

for (i in 1:length(names)){
  reg_mkt_m_sub[[1]][[i]] <- lm(a_m_sub[[1]][,names[i]] ~ b_m_sub[[1]])
  reg_mkt_m_sub[[2]][[i]] <- lm(a_m_sub[[2]][,names[i]] ~ b_m_sub[[2]])
  reg_mkt_m_sub[[3]][[i]] <- lm(a_m_sub[[3]][,names[i]] ~ b_m_sub[[3]])
  
  names(reg_mkt_m_sub[[1]][[i]]$coefficients) <- 
    c(names(reg_mkt_m_sub[[1]][[i]]$coefficients[1]), "Mkt-RF")
  names(reg_mkt_m_sub[[2]][[i]]$coefficients) <- 
    c(names(reg_mkt_m_sub[[2]][[i]]$coefficients[1]), "Mkt-RF")
  names(reg_mkt_m_sub[[3]][[i]]$coefficients) <- 
    c(names(reg_mkt_m_sub[[3]][[i]]$coefficients[1]), "Mkt-RF")
  
  reg_FF3_m_sub[[1]][[i]] <- 
    lm(a_m_sub[[1]][,names[i]] ~ b_m_sub[[1]] + b1_m_sub[[1]] + b2_m_sub[[1]])
  reg_FF3_m_sub[[2]][[i]] <- 
    lm(a_m_sub[[2]][,names[i]] ~ b_m_sub[[2]] + b1_m_sub[[2]] + b2_m_sub[[2]])
  reg_FF3_m_sub[[3]][[i]] <- 
    lm(a_m_sub[[3]][,names[i]] ~ b_m_sub[[3]] + b1_m_sub[[3]] + b2_m_sub[[3]])
  
  reg_mkt_se_m_sub[[1]][[i]] <- 
    reg_mkt_m_sub[[1]][[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
  reg_mkt_se_m_sub[[2]][[i]] <- 
    reg_mkt_m_sub[[2]][[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
  reg_mkt_se_m_sub[[3]][[i]] <- 
    reg_mkt_m_sub[[3]][[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
}

for (period in 1:as.numeric(count(sub_dates))) {
  for (i in 1:length(names)) {
    reg_output_m_sub[[period]]["alpha_mkt", i] <- 
      reg_mkt_m_sub[[period]][[i]]$coefficients[1]
    reg_output_m_sub[[period]]["beta_mkt", i] <- 
      reg_mkt_m_sub[[period]][[i]]$coefficients[2]
    reg_output_m_sub[[period]]["RMSE", i] <- 
      sqrt(mean(residuals(reg_mkt_m_sub[[period]][[i]])^2))
    reg_output_m_sub[[period]]["AR", i] <- sqrt(trading_months) * 
      reg_output_m_sub[[period]]["alpha_mkt", i] / 
      reg_output_m_sub[[period]]["RMSE", i]
    reg_output_m_sub[[period]]["alpha_FF3", i] <- 
      reg_FF3_m_sub[[period]][[i]]$coefficients[1]
    reg_output_m_sub[[period]]["alpha_FF3_se", i] <- 
      coeftest(reg_FF3_m_sub[[period]][[i]], 
               vcovHC(reg_FF3_m_sub[[period]][[i]], 
                      type = "HC"))[1,2]
    alpha_stars_m_sub[[period]][i] <- ff3_alpha_stars(reg_FF3_m_sub[[period]][[i]])
  }
}

stargazer(reg_mkt_m_sub[[1]][[1]], reg_mkt_m_sub[[1]][[2]],
          reg_mkt_m_sub[[1]][[3]], reg_mkt_m_sub[[1]][[4]],
          reg_mkt_m_sub[[2]][[1]], reg_mkt_m_sub[[2]][[2]],
          reg_mkt_m_sub[[2]][[3]], reg_mkt_m_sub[[2]][[4]],
          reg_mkt_m_sub[[3]][[1]], reg_mkt_m_sub[[3]][[2]],
          reg_mkt_m_sub[[3]][[3]], reg_mkt_m_sub[[3]][[4]],
          se = list(reg_mkt_se_m_sub[[1]][[1]], reg_mkt_se_m_sub[[1]][[2]],
                    reg_mkt_se_m_sub[[1]][[3]], reg_mkt_se_m_sub[[1]][[4]],
                    reg_mkt_se_m_sub[[2]][[1]], reg_mkt_se_m_sub[[2]][[2]],
                    reg_mkt_se_m_sub[[2]][[3]], reg_mkt_se_m_sub[[2]][[4]],
                    reg_mkt_se_m_sub[[3]][[1]], reg_mkt_se_m_sub[[3]][[2]],
                    reg_mkt_se_m_sub[[3]][[3]], reg_mkt_se_m_sub[[3]][[4]]),
          type = "text", 
          keep.stat = c("n", "rsq", "ser"), df = FALSE, 
          dep.var.labels = c("1926 - 1957", "1957 - 1988", "1988 - 2019"),
          column.labels = rep(names_clean, as.numeric(count(sub_dates))),
          table.layout = "-d#c-t-s-a-",
          covariate.labels = c("Mkt-RF", "Alpha (&#945;)"),
          add.lines = list(c("Standard Deviation", 
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[1]]["SD",])),
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[2]]["SD",])),
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[3]]["SD",]))),
                           c("Sharpe Ratio", 
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[1]]["SR",])),
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[2]]["SR",])),
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[3]]["SR",]))),
                           c("Appr Ratio", 
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[1]]["AR",])),
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[2]]["AR",])),
                             sprintf("%.2f", as.numeric(reg_output_m_sub[[3]]["AR",]))),
                           c("Alpha FF3", 
                             paste(sprintf("%.2f", 
                                           as.numeric(reg_output_m_sub[[1]]["alpha_FF3",])),
                                   alpha_stars_m_sub[[1]], sep = ""),
                             paste(sprintf("%.2f", 
                                           as.numeric(reg_output_m_sub[[2]]["alpha_FF3",])),
                                   alpha_stars_m_sub[[2]], sep = ""),
                             paste(sprintf("%.2f", 
                                           as.numeric(reg_output_m_sub[[3]]["alpha_FF3",])),
                                   alpha_stars_m_sub[[3]], sep = "")),
                           c("", 
                             paste("(", 
                                   sprintf("%.2f", 
                                           as.numeric(reg_output_m_sub[[1]]["alpha_FF3_se",])),
                                   ")", sep = ""),
                             paste("(", 
                                   sprintf("%.2f", 
                                           as.numeric(reg_output_m_sub[[2]]["alpha_FF3_se",])),
                                   ")", sep = ""),
                             paste("(", 
                                   sprintf("%.2f", 
                                           as.numeric(reg_output_m_sub[[3]]["alpha_FF3_se",])),
                                   ")", sep = ""))),
          digits = 2, out = "table.htm")

################################################################################
#***************************** Custom Rebalancing *****************************#

# Set Up List to Store Outputs from Different Frequencies
min_frequ <- 5
max_frequ <- 50
trading_cost <- 0.1

reg_output_c <- vector(mode = "list", length = max_frequ - min_frequ + 1)
reg_output_cost_c <- vector(mode = "list", length = max_frequ - min_frequ + 1)

final_return_c <- data.frame(matrix(nrow = max_frequ - min_frequ + 1, 
                                    ncol = length(names),
                                    dimnames = list(min_frequ:max_frequ, names)))

final_return_cost_c <- data.frame(matrix(nrow = max_frequ - min_frequ + 1, 
                                         ncol = length(names),
                                         dimnames = list(min_frequ:max_frequ, names)))

# Loop to Test Different Frequencies
for (frequ in min_frequ:max_frequ) {
  # Calculate Custom Variances
  n_custom <- as.integer(n_days / frequ)
  var_c <- data.frame(matrix(ncol = 2, nrow = n_custom))
  colnames(var_c) <- c("date", "var")
  
  var_c <- var_c %>% mutate(date = FF_daily$date[1:n_custom])
  
  for (i in 1:n_days) {
    if(i %% frequ == 0) {
      j <- i / frequ
      var_c$date[j] <- FF_daily$date[i]
      var_c$var[j] <- var(FF_daily$`mkt-rf`[(i - frequ + 1):i]) * (frequ - 1) / frequ
    }
  }
  
  # Set Up Return Vector and Calculate Market Return
  returns_c <- data.frame(matrix(ncol = length(names) + 4, nrow = (n_custom-1)))
  colnames(returns_c) <- c("date", "mkt", "rf", "mkt-rf", names)
  returns_c$date <- var_c$date[-1]
  
  for (i in 1:(n_days - frequ)) {
    if(i %% frequ == 0) {
      j <- i / frequ
      ret_temp <- 0
      rf_temp <- 0
      for (a in 1:frequ) {
        ret_temp <- ((1 + FF_daily$mkt[i + a] / 100) * (1 + ret_temp / 100) - 1) * 100
        rf_temp <- ((1 + FF_daily$rf[i + a] / 100) * (1 + rf_temp / 100) - 1) * 100
      }
      returns_c$mkt[j] <- ret_temp
      returns_c$rf[j] <- rf_temp
    }
  }
  
  # Calculate Squared Log Return
  returns_c <- returns_c %>% mutate("mkt-rf" = mkt - rf)
  u_sq_c <- c(1:n_custom)
  ret_temp <- 0
  rf_temp <- 0
  for (a in 1:frequ) {
    ret_temp <- ((1 + FF_daily$mkt[a] / 100) * (1 + ret_temp / 100) - 1) * 100
    rf_temp <- ((1 + FF_daily$rf[a] / 100) * (1 + rf_temp / 100) - 1) * 100
  }
  u_sq_c[1] <- log(1 + (ret_temp - rf_temp) / 100)^2
  u_sq_c[2:n_custom] <- log( 1 + returns_c$`mkt-rf` / 100)^2
  
  # Calculate ARIMA Variance
  variance_ts_c <- xts(var_c$var, order.by = var_c$date)
  ARIMA_model_c <- auto.arima(variance_ts_c, stepwise = FALSE, approximation = FALSE,
                              max.p = 10, max.q = 10, ic = "aic", test = "kpss")
  
  var_c$ARIMA_var <- c(fitted(ARIMA_model_c)[-1], forecast(ARIMA_model_c, h = 1)$mean)
  
  # Optimize ESA Parameter
  ESA_function_c <- function(lambda)
  {
    ESA_var <- c(1:n_custom)
    ESA_var[1] <- u_sq_c[1]
    for(i in 2:n_custom) {
      ESA_var[i] <- lambda * ESA_var[i-1] + (1 - lambda) * u_sq_c[i]
    }
    ESA_likelihood <- c(1:(n_custom-1))
    for(i in 1:(n_custom-1)) {
      ESA_likelihood[i] <- -log(ESA_var[i])-u_sq_c[i+1]/ESA_var[i]
    }
    return(sum(ESA_likelihood))
  }
  ESA_max_c <- optimize(ESA_function_c, interval = c(0, 1), maximum = TRUE, 
                         tol = 0.000000000000001)
  lambda_c <- ESA_max_c$maximum
  
  # Calculate ESA Variance
  var_c$ESA_var <- c(1:n_custom)
  var_c$ESA_var[1] <- var_c$var[1]
  for (i in 2:n_custom) {
    var_c$ESA_var[i] <- lambda_c * var_c$ESA_var[i - 1] + (1 - lambda_c) * 
      var_c$var[i]
  }

  # Optimize GARCH Parameters
  GARCH_function_c <- function(alpha, beta)
  {
    omega <- max(0,mean(u_sq_c)*(1-alpha-beta))
    GARCH_var <- c(1:n_custom)
    GARCH_var[1] <- u_sq_c[1]
    for(i in 2:n_custom) {
      GARCH_var[i] <- omega + beta*GARCH_var[i-1] + alpha*u_sq_c[i]
    }
    GARCH_likelihood <- c(1:(n_custom - 1))
    for(i in 1:(n_custom-1)) {
      GARCH_likelihood[i] <- -log(GARCH_var[i])-u_sq_c[i+1]/GARCH_var[i]
    }
    return (sum(GARCH_likelihood))
  }
  
  GARCH_max_c <- optimx(c(0.1, 0.9), function(x) GARCH_function_c(x[1], x[2]), 
                      method = "Nelder-Mead", control = list(maximize = TRUE))
  alpha_c <- GARCH_max_c$p1
  beta_c <- GARCH_max_c$p2
  omega_c <- max(0, mean(u_sq_c) * (1 - alpha_c-beta_c))
  
  # Calculate GARCH Variance
  var_c$GARCH_var <- c(1:n_custom)
  var_c$GARCH_var[1] <- u_sq_c[1]
  for (i in 2:n_custom) {
    var_c$GARCH_var[i] <- omega_c + alpha_c * u_sq_c[i] + beta_c * var_c$GARCH_var[i-1]
  }
  
  # Calculate c with Midnight Formula
  c_c <- data.frame(matrix(ncol = length(names)))
  colnames(c_c) <- names
  
  a_qe_c <- c(1:length(names))
  b_qe_c <- c(1:length(names))
  c_qe_c <- c(1:length(names))
  
  # Usage of Population Variance, since the prefactor cancels out on aggregate
  for (i in 1:length(names)) {
    a_qe_c[i] <- var(returns_c$`mkt-rf` / var_c[-n_custom,var_names[i]])
    b_qe_c[i] <- 2 * cov(returns_c$`mkt-rf` / var_c[-n_custom,var_names[i]], returns_c$rf)
    c_qe_c[i] <- var(returns_c$rf) - var(returns_c$mkt)
    
    c_c[i] <- 1 / (2 * a_qe_c[i]) * 
      (-b_qe_c[i] + sqrt((b_qe_c[i])^2 - 4 * a_qe_c[i] * c_qe_c[i]))
  }
  
  # Calculate Weights and Returns
  w_abs_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
  colnames(w_abs_c) <- names
  
  returns_cost_c <- data.frame(matrix(ncol = length(names) + 1, nrow = (n_custom-1)))
  colnames(returns_cost_c) <- c("date", names)
  returns_cost_c <- returns_cost_c %>% mutate(date = FF_daily$date[1:(n_custom-1)])
  
  weights_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
  colnames(weights_c) <- names
  
  for (i in 1:length(names)) {
    weights_c[,names[i]] <- as.numeric(c_c[i]) / var_c[-n_custom,var_names[i]]
    returns_c[,names[i]] <- weights_c[,names[i]] * returns_c$`mkt-rf` + returns_c$rf
    w_abs_c[1,names[i]] <- 1
    w_abs_c[-1,names[i]] <- abs(diff(weights_c[,names[i]]))
    returns_cost_c[,names[i]] <- returns_c[,names[i]] - w_abs_c[,names[i]] * trading_cost
  }

  # Calculate Cumulative Return
  cum_ret_c <- data.frame(matrix(ncol = length(names) + 2, nrow = n_custom))
  colnames(cum_ret_c) <- c("date", "mkt", names)
  cum_ret_c <- cum_ret_c %>% mutate(date = var_c$date)
  
  cum_ret_cost_c <- data.frame(matrix(ncol = length(names) + 2, nrow = n_custom))
  colnames(cum_ret_cost_c) <- c("date", "mkt", names)
  cum_ret_cost_c <- cum_ret_cost_c %>% mutate(date = var_c$date)
  
  cum_ret_c[1,-1] <- 1
  cum_ret_cost_c[1,-1] <- 1
  
  for (i in 2:n_custom) {
    cum_ret_c$mkt[i] <- cum_ret_c$mkt[i-1] * (1 + (returns_c$mkt[i-1] / 100))
    cum_ret_cost_c$mkt[i] <- cum_ret_c$mkt[i]
    for (j in 1:length(names)) {
      cum_ret_c[i,names[j]] <- cum_ret_c[i-1, names[j]] * 
        (1 + (returns_c[i-1, names[j]]/100))
      cum_ret_cost_c[i,names[j]] <- cum_ret_cost_c[i-1, names[j]] * 
        (1 + (returns_cost_c[i-1, names[j]] / 100))
    }
  }
  
  final_return_c[frequ-min_frequ+1,] <- cum_ret_c[n_custom,3:6]
  final_return_cost_c[frequ-min_frequ+1,] <- cum_ret_cost_c[n_custom,3:6]
  
  # Regressions to Determine Alpha, Beta, etc.
  reg_mkt_c <- vector(mode = "list", length = length(names))
  reg_mkt_cost_c <- vector(mode = "list", length = length(names))
  b_c <- trading_year / frequ * (returns_c$`mkt-rf`)
  
  for (i in 1:length(names)) {
    a <- trading_year / frequ * (returns_c[, names[i]] - returns_c$rf)
    reg_mkt_c[[i]] <- lm(a ~ b_c)
    a <- trading_year / frequ * (returns_cost_c[, names[i]] - returns_c$rf)
    reg_mkt_cost_c[[i]] <- lm(a ~ b_c)
  }
  
  output_names_c <- c("alpha_mkt", "RMSE", "Appr_Ratio")
  mod_num <- frequ - min_frequ + 1
  
  reg_output_c[[mod_num]] <- data.frame(matrix(ncol = length(names), 
                                               nrow = length(output_names_c),
                                               dimnames = list(output_names_c, names)))
  
  reg_output_cost_c[[mod_num]] <- data.frame(matrix(ncol = length(names), 
                                               nrow = length(output_names_c),
                                               dimnames = list(output_names_c, names)))
  
  for (i in 1:length(names)) {
    reg_output_c[[mod_num]]["alpha_mkt", i] <- reg_mkt_c[[i]]$coefficients[1]
    reg_output_c[[mod_num]]["RMSE", i] <- sqrt(mean(residuals(reg_mkt_c[[i]])^2))
    reg_output_c[[mod_num]]["Appr_Ratio", i] <- sqrt(trading_year / frequ) * 
      reg_output_c[[mod_num]]["alpha_mkt", i] / reg_output_c[[mod_num]]["RMSE", i]
  }
  
  for (i in 1:length(names)) {
    reg_output_cost_c[[mod_num]]["alpha_mkt", i] <- reg_mkt_cost_c[[i]]$coefficients[1]
    reg_output_cost_c[[mod_num]]["RMSE", i] <- sqrt(mean(residuals(reg_mkt_cost_c[[i]])^2))
    reg_output_cost_c[[mod_num]]["Appr_Ratio", i] <- sqrt(trading_year / frequ) * 
      reg_output_cost_c[[mod_num]]["alpha_mkt", i] / reg_output_cost_c[[mod_num]]["RMSE", i]
  }
  print(paste("Progress: ", frequ - min_frequ + 1, "/",
              max_frequ - min_frequ + 1, sep = ""))
}

# Consolidate Alphas and Appraisal Ratios
alpha_c <- data.frame(matrix(ncol = length(names),
                             nrow = max_frequ - min_frequ + 1,
                             dimnames = list(min_frequ:max_frequ, names)))

alpha_cost_c <- data.frame(matrix(ncol = length(names), 
                                  nrow = max_frequ - min_frequ + 1,
                                  dimnames = list(min_frequ:max_frequ, names)))

appr_c <- data.frame(matrix(ncol = length(names), 
                            nrow = max_frequ - min_frequ + 1,
                            dimnames = list(min_frequ:max_frequ, names)))

appr_cost_c <- data.frame(matrix(ncol = length(names), 
                                 nrow = max_frequ - min_frequ + 1,
                                 dimnames = list(min_frequ:max_frequ, names)))

# Plot Alphas, Appraisal Ratios, and Cumulative Performance (and group graphs)
for (i in 1:(max_frequ - min_frequ + 1)) {
  for (j in 1:length(names)) {
    alpha_c[i,j] <- reg_output_c[[i]][1,j]
    appr_c[i,j] <- reg_output_c[[i]][3,j]
    alpha_cost_c[i,j] <- reg_output_cost_c[[i]][1,j]
    appr_cost_c[i,j] <- reg_output_cost_c[[i]][3,j]
  }
}

fig04 <- ggplot(alpha_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, col = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, col = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH"), size = line_size) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "none", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Alpha Depending on Holding Period") +
  xlab("") + ylab("") +
  ylim(0, 6) +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

fig15 <- ggplot(appr_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, col = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, col = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH"), size = line_size) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "none", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Appraisal Ratio Depending on Holding Period") +
  xlab("") + ylab("") +
  ylim(0, 0.5) +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

fig05 <- ggplot(alpha_cost_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, col = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, col = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH"), size = line_size) +
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
  ggtitle("Alpha Depending on Holding Period (10bps Transactions Costs)") +
  xlab("") + ylab("") +
  ylim(0, 6) +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

fig16 <- ggplot(appr_cost_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, col = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, col = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH"), size = line_size) +
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
  ggtitle("Appraisal Ratio Depending on Holding Period (10bps Transaction Costs)") +
  xlab("") + ylab("") +
  ylim(0, 0.5) +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

fig19 <- ggplot(final_return_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, col = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, col = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH"), size = line_size) +
  theme_bw(base_family = "Times New Roman") +
  theme(legend.position = "none", 
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(25000,50000,75000), limits = c(0,90000)) +
  ggtitle("Cumulative Performance Depending on Holding Period") +
  xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

fig20 <- ggplot(final_return_cost_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance"), size = line_size) +
  geom_line(aes(y=ARIMA_var_managed, col = "ARIMA"), size = line_size) +
  geom_line(aes(y=ESA_var_managed, col = "ESA"), size = line_size) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH"), size = line_size) +
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
  scale_y_continuous(breaks = c(25000,50000,75000), limits = c(0,90000)) +
  ggtitle("Cumulative Performance Depending on Holding Period (10 bps Trading Cost)") +
  xlab("") + ylab("") +
  scale_color_manual(name = "", 
                     values = c("Realized Variance" = blue_col[4],
                                "ARIMA" = blue_col[6],
                                "ESA" = blue_col[8],
                                "GARCH" = blue_col[9]),
                     breaks = c("Realized Variance",
                                "ARIMA", "ESA", "GARCH")) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

grid.arrange(fig04, fig05, nrow = 2, heights = c(0.46, 0.54))
grid.arrange(fig15, fig16, nrow = 2, heights = c(0.46, 0.54))
grid.arrange(fig19, fig20, nrow = 2, heights = c(0.46, 0.54))
