library(ggplot2)
library(ggfortify)
library(tidyverse)
library(lubridate)
library(patchwork)
library(TSA)
library(lmtest)
library(forecast)
library(fUnitRoots)
library(tseries)
library(dynlm)
library(vars)



fedfunds <- read.csv("FEDFUNDS.csv")
mortgage <- read.csv("MORTGAGE30US.csv")
houses <- read.csv("HOUST1F.csv")

head(fedfunds)
head(mortgage)
head(houses)

# convert to date type
fedfunds$date <- as.Date(fedfunds$observation_date)
houses$date   <- as.Date(houses$observation_date)
mortgage$date <- as.Date(mortgage$observation_date)


# convert to monthly with mean
mortgage <- mortgage %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(MORTGAGE30US = mean(MORTGAGE30US, na.rm = TRUE),.groups = "drop")

head(mortgage)

# datefiltering
startDate <- as.Date("1986-01-01")
endDate <- as.Date("2026-01-01")

fedfunds <- fedfunds %>% filter(date >= startDate, date <= endDate) %>% arrange(date)
mortgage <- mortgage %>% filter(date >= startDate, date <= endDate) %>% arrange(date)
houses <- houses %>% filter(date >= startDate, date <= endDate) %>% arrange(date)

fedfunds_ts <- ts(fedfunds$FEDFUNDS, start = c(1986,1), frequency = 12)
mortgage_ts <- ts(mortgage$MORTGAGE30US, start = c(1986,1), frequency = 12)
houses_ts <- ts(houses$HOUST1F, start = c(1986,1), frequency = 12)

fedfunds_ts <- window(fedfunds_ts, end = c(2025, 12))
houses_ts <- window(houses_ts, end = c(2025, 12))


#### Initial Plots, differencing, logs, and stationarity checks ####

autoplot(fedfunds_ts) # not multiplicative, so not log, but mean in changing
# take diff
fedfunds_diff = diff(fedfunds_ts, differences = 1)
autoplot(fedfunds_diff)
# Stationary
adfTest(fedfunds_diff, type="nc")
kpss.test(fedfunds_diff)


# mortgage ts
autoplot(mortgage_ts)
# take log diff
mortgage_ld = diff(log(mortgage_ts), differences = 1)
autoplot(mortgage_ld)
# stationary
adfTest(mortgage_ld, type="nc")
kpss.test(mortgage_ld)

#intervention for the 2008 housing crash
houses_time = time(houses_ts)
interventionYear = 2007
interventionMonth = 3
resesTime = interventionYear + (interventionMonth - 1) / 12
interventionIndex = which.min(abs(houses_time - resesTime))   # robust to float error
recesDummy = ifelse(seq_along(houses_time) >= interventionIndex, 1, 0)
rampDummy = ifelse(seq_along(houses_time) >= interventionIndex, seq_along(houses_time) - interventionIndex + 1, 0)
recessionPredictor = cbind(step = recesDummy, ramp = rampDummy)

# housing tf
autoplot(houses_ts)
# log diff
houses_ld = diff(log(houses_ts), differences = 1)
autoplot(houses_ld)
# Stationary
adfTest(houses_ld, type="nc")
kpss.test(houses_ld)


# all stationary after differencing/log differencing


# plot together with patchwork
fedfunds_plot = autoplot(fedfunds_diff, color="firebrick") +
  labs(title="Differenced: Federal Funds Rate", y="", x="") +
  theme_minimal()
mortgage_plot = autoplot(mortgage_ld, color="steelblue") +
  labs(title="Log-Differenced: Mortgage Rate", y="", x="") +
  theme_minimal()
houses_plot = autoplot(houses_ld, color="darkgreen") +
  labs(title="Log-Differenced: Housing Starts", y="", x="") +
  theme_minimal()
fedfunds_plot / mortgage_plot / houses_plot



##### acf #####

fedfunds_acf = autoplot(Acf(fedfunds_diff)) +
  labs(title="Federal Funds Rate ACF", x="Lag", y="")
# largest spike at lag 1, some potential seasonality


mortgage_acf = autoplot(Acf(mortgage_ld)) +
  labs(title="Mortgage Rate ACF", x="Lag", y="")
# large spike at lag1, another at 21, probably spurious - MA(1)


houses_acf = autoplot(Acf(houses_ld)) +
  labs(title="Housing Starts ACF", x="Lag", y="")
# MA, maybe some seasonality


# combine with patchwork
fedfunds_acf / mortgage_acf / houses_acf


#### pacf ####
fedfunds_pacf = autoplot(Pacf(fedfunds_diff)) +
  labs(title="Federal Funds Rate PACF", x="Lag", y="")
# ARMA process

mortgage_pacf = autoplot(Pacf(mortgage_ld)) +
  labs(title="Mortgage Rate PACF", x="Lag", y="")
# probably MA

houses_pacf = autoplot(Pacf(houses_ld)) +
  labs(title="Housing Starts PACF", x="Lag", y="")
# probably MA

fedfunds_pacf / mortgage_pacf / houses_pacf


# show both together
fedfunds_acf / fedfunds_pacf
mortgage_acf / mortgage_pacf
houses_acf / houses_pacf




##### eacf #####
eacf(fedfunds_diff)
# inconclusive, maybe ARIMA(3,1,3)

eacf(mortgage_ld)
# (0,1,1) / (1,1,1) ?

eacf(houses_ld)
# (1,1,2) / (0,1,3) ?


###### arima ######

# federal funds
fed_fit1 <- Arima(fedfunds_ts, order = c(2,1,0))
fed_fit1
coeftest(fed_fit1) # both terms significant
checkresiduals(fed_fit1) # plots don't look good, clearly nonstationary
Box.test(residuals(fed_fit1), lag=24, type="L") # confirms residuals are not white noise
jarque.bera.test(residuals(fed_fit1))   #j-b test, not normal

# fed_fit1 NOT good

fed_fit2 <- Arima(fedfunds_ts, order=c(2,1,1))
fed_fit2
coeftest(fed_fit2) # all significant, MA slightly less
checkresiduals(fed_fit2) # not good
Box.test(residuals(fed_fit2), lag=24, type="L") # residuals white noise
jarque.bera.test(residuals(fed_fit2))   #j-b test fails badly


# fed_fit2 better but not normal

fed_fit3 <- Arima(fedfunds_ts, order=c(0,1,1))
fed_fit3
coeftest(fed_fit3) # significant
checkresiduals(fed_fit3) # not good
Box.test(residuals(fed_fit3), lag=24, type="L") # residuals not white noise
jarque.bera.test(residuals(fed_fit3))   #j-b test 


# fed_fit3 NOT good

fed_fit_auto <- auto.arima(fedfunds_ts)
fed_fit_auto
coeftest(fed_fit_auto) # significant
checkresiduals(fed_fit_auto) # still not great
Box.test(residuals(fed_fit_auto), lag=24, type="Ljung-Box") # residuals not white noise
jarque.bera.test(residuals(fed_fit_auto))   #j-b test 

fed_fit_autob <- auto.arima(fedfunds_ts, ic="bic")
fed_fit_autob
coeftest(fed_fit_autob) # significant
checkresiduals(fed_fit_autob) # still not great
Box.test(residuals(fed_fit_autob), lag=24, type="Ljung-Box")
jarque.bera.test(residuals(fed_fit_autob))   #j-b test 

# same as AIC

# AIC/BIC comparison
AIC(fed_fit1, fed_fit2, fed_fit3, fed_fit_auto, fed_fit_autob)
BIC(fed_fit1, fed_fit2, fed_fit3, fed_fit_auto, fed_fit_autob)

# fedfit2 white noise, but not normal but still best

# mortgage rates
mort_fit1 <- Arima(log(mortgage_ts), order=c(2,1,0))
mort_fit1
coeftest(mort_fit1) # significant
checkresiduals(mort_fit1) # pretty good
Box.test(residuals(mort_fit1), lag=24, type="Ljung-Box") # white noise
jarque.bera.test(residuals(mort_fit1))   #j-b test 


# mort_fit1 is good

mort_fit2 <- Arima(log(mortgage_ts), order=c(3,1,0))
mort_fit2
coeftest(mort_fit2) # AR3 less significant
checkresiduals(mort_fit2) # decent
Box.test(residuals(mort_fit2), lag=24, type="Ljung-Box") # white noise
jarque.bera.test(residuals(mort_fit2))   #j-b test 


# mort_fit2 is good


mort_fit3 <- Arima(log(mortgage_ts), order=c(0,1,1)) 
mort_fit3
coeftest(mort_fit3) #significant
checkresiduals(mort_fit3) # looks good
Box.test(residuals(mort_fit3), lag=24, type="Ljung-Box") # white noise
jarque.bera.test(residuals(mort_fit3))   #j-b test 


# mort_fit3 is good

mort_auto <- auto.arima(log(mortgage_ts))
mort_auto
coeftest(mort_auto) # also went with ARIMA(0,1,1)
checkresiduals(mort_auto)
Box.test(residuals(mort_auto), lag=24, type="Ljung-Box")
jarque.bera.test(residuals(mort_auto))   #j-b test 


mort_autob <- auto.arima(log(mortgage_ts), ic="bic")
mort_autob
coeftest(mort_autob) # also went with ARIMA(0,1,1)
checkresiduals(mort_autob)
Box.test(residuals(mort_autob), lag=24, type="Ljung-Box")
jarque.bera.test(residuals(mort_autob))   #j-b test 


# AIC/BIC comparison
AIC(mort_fit1, mort_fit2, mort_fit3, mort_auto, mort_autob)
BIC(mort_fit1, mort_fit2, mort_fit3, mort_auto, mort_autob)

# Auto.arima chose same model as mort_fit_3. Based on parsiomy and results, looks good


# housing
houses_fit1 <- Arima(log(houses_ts), order=c(0,1,1), xreg=recessionPredictor)
houses_fit1
coeftest(houses_fit1) # significant
checkresiduals(houses_fit1) # decent but misses seasonal structure
Box.test(residuals(houses_fit1), lag=24, type="Ljung-Box") # borderline
jarque.bera.test(residuals(houses_fit1))   #j-b test 

# good, could be better

houses_fit2 <- Arima(log(houses_ts), order=c(1,1,0), xreg=recessionPredictor)
houses_fit2
coeftest(houses_fit2) # significant
checkresiduals(houses_fit2) # misses seasonal structure
Box.test(residuals(houses_fit2), lag=24, type="Ljung-Box") # fails
jarque.bera.test(residuals(houses_fit2))   #j-b test 

# not white noise, NOT good

houses_fit3 <- Arima(log(houses_ts), order=c(1,1,1), xreg=recessionPredictor)
houses_fit3
coeftest(houses_fit3) # check significance
checkresiduals(houses_fit3) # still missing seasonal structure
Box.test(residuals(houses_fit3), lag=24, type="Ljung-Box") # borderline
jarque.bera.test(residuals(houses_fit3))   #j-b test 

# could be better

# add seasonal terms
houses_sarima1 <- Arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,0), period=12),
                        xreg=recessionPredictor)
houses_sarima1
coeftest(houses_sarima1) # seasonal term not significant
checkresiduals(houses_sarima1) # improved with seasonal term
Box.test(residuals(houses_sarima1), lag=24, type="Ljung-Box")
jarque.bera.test(residuals(houses_sarima1))   #j-b test 


# still not great

houses_sarima2 <- Arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,1), period=12),
                        xreg=recessionPredictor)
houses_sarima2
coeftest(houses_sarima2) # all terms significant
checkresiduals(houses_sarima2) # good
Box.test(residuals(houses_sarima2), lag=24, type="Ljung-Box")
jarque.bera.test(residuals(houses_sarima2))   #j-b test 

# best yet

recessionPredictor_scaled <- recessionPredictor
recessionPredictor_scaled[, "ramp"] <- recessionPredictor[, "ramp"] / 227

houses_sarima3 <- Arima(log(houses_ts), order = c(1,1,1),
                        seasonal = list(order = c(1,0,1), period = 12),
                        xreg = recessionPredictor_scaled,
                        method = "CSS")
houses_sarima3
coeftest(houses_sarima3) # AR and MA both insignificant
checkresiduals(houses_sarima3)
Box.test(residuals(houses_sarima3), lag=24, type="Ljung-Box") # passes
jarque.bera.test(residuals(houses_sarima3))   #j-b test 

houses_sarima4 <- Arima(log(houses_ts), order=c(2,1,1),
                        seasonal=list(order=c(1,0,1), period=12),
                        xreg=recessionPredictor)
houses_sarima4
coeftest(houses_sarima4) # ma1 marginally significant
checkresiduals(houses_sarima4) # good
Box.test(residuals(houses_sarima4), lag=24, type="Ljung-Box") # passes
jarque.bera.test(residuals(houses_sarima4))   #j-b test 

# best residuals yet, less significant coeftest

houses_sarima5 <- Arima(log(houses_ts), order=c(1,1,2),
                        seasonal=list(order=c(1,0,1), period=12),
                        xreg=recessionPredictor)
houses_sarima5
coeftest(houses_sarima5) # check significance
checkresiduals(houses_sarima5) # good
Box.test(residuals(houses_sarima5), lag=24, type="Ljung-Box") # passes
jarque.bera.test(residuals(houses_sarima5))   #j-b test 


# good

houses_auto <- auto.arima(log(houses_ts), xreg=recessionPredictor)
houses_auto
coeftest(houses_auto) # more complex model, a lot of insignificant terms
checkresiduals(houses_auto) # good
Box.test(residuals(houses_auto), lag=24, type="Ljung-Box") # passes
jarque.bera.test(residuals(houses_auto))   #j-b test 


houses_autob <- auto.arima(log(houses_ts), ic="bic", xreg=recessionPredictor)
houses_autob
coeftest(houses_autob) # more complex model, a lot of insignificant terms
checkresiduals(houses_autob) # good
Box.test(residuals(houses_autob), lag=24, type="Ljung-Box") 
jarque.bera.test(residuals(houses_autob))   #j-b test 

# best box test

# AIC/BIC comparison
AIC(houses_fit1, houses_fit2, houses_fit3,
    houses_sarima1, houses_sarima2,
    houses_sarima4, houses_sarima5, houses_auto, houses_autob)

BIC(houses_fit1, houses_fit2, houses_fit3,
    houses_sarima1, houses_sarima2,
    houses_sarima4, houses_sarima5, houses_auto, houses_autob)

# houses_sarima_2 the best by all metrics


##### Backtesting #####

source("backtest.R")

# federal funds
backtest(fed_fit1, fedfunds_ts, h=1, orig=.8*length(fedfunds_ts))
backtest(fed_fit2, fedfunds_ts, h=1, orig=.8*length(fedfunds_ts))
backtest(fed_fit3, fedfunds_ts, h=1, orig=.8*length(fedfunds_ts))
#backtest(fed_fit_auto, fedfunds_ts, h=1, orig=floor(.8*length(fedfunds_ts))) # doesn't work, over parameterized
backtest(fed_fit_autob, fedfunds_ts, h=1, orig=floor(.8*length(fedfunds_ts)))

# fed_fit2 and autob best results, choose fed_fit2 as best model

# mortgage rates (log transformed)
backtest(mort_fit1, log(mortgage_ts), h=1, orig=.8*length(mortgage_ts))
backtest(mort_fit2, log(mortgage_ts), h=1, orig=.8*length(mortgage_ts))
backtest(mort_fit3, log(mortgage_ts), h=1, orig=.8*length(mortgage_ts))
backtest(mort_auto, log(mortgage_ts), h=1, orig=.8*length(mortgage_ts)) 
backtest(mort_autob, log(mortgage_ts), h=1, orig=.8*length(mortgage_ts))

# mort_fit3 and autoarima models best - matches above

# housing starts (log transformed)
backtest(houses_fit1, log(houses_ts), h=1, orig=.8*length(houses_ts))
backtest(houses_fit2, log(houses_ts), h=1, orig=.8*length(houses_ts))
backtest(houses_fit3, log(houses_ts), h=1, orig=.8*length(houses_ts))
backtest(houses_sarima1, log(houses_ts), h=1, orig=.8*length(houses_ts))
backtest(houses_sarima2, log(houses_ts), h=1, orig=.8*length(houses_ts))
#backtest(houses_sarima3, log(houses_ts), h=1, orig=.8*length(houses_ts)) # BT doesn't work
backtest(houses_sarima4, log(houses_ts), h=1, orig=.8*length(houses_ts))
#backtest(houses_sarima5, log(houses_ts), h=1, orig=.8*length(houses_ts)) # BT doesn't work
#backtest(houses_auto, log(houses_ts), h=1, orig=.8*length(houses_ts)) # BT doesn't work
backtest(houses_autob, log(houses_ts), h=1, orig=.8*length(houses_ts))

# sarima model 2 best on MAE, competitive with RMSE, best overall

##### Forecasting #####


H <- 24

## future regressors for the housing models
last_ramp <- tail(recessionPredictor[, "ramp"], 1)
future_xreg <- cbind(step = rep(1, H),
                     ramp = (last_ramp + 1):(last_ramp + H))

# scaled version for houses_sarima3 (fit with ramp / 227)
future_xreg_scaled <- cbind(step = rep(1, H),
                            ramp = ((last_ramp + 1):(last_ramp + H)) / 227)

## logs back to original uits
nolog_forecast <- function(forecast) 
  {forecast$mean <- exp(forecast$mean)
  forecast$lower <- exp(forecast$lower)
  forecast$upper <- exp(forecast$upper)
  forecast$x <- exp(forecast$x)
  if (!is.null(forecast$fitted)) forecast$fitted <- exp(forecast$fitted)
  forecast}

## Federal funds (no transform, no xreg)
autoForecastFed <- forecast(fed_fit_auto, h = H); autoplot(autoForecastFed)
autoForecastFedB <- forecast(fed_fit_autob, h = H); autoplot(autoForecastFedB)

fedForecast1 <- forecast(fed_fit1, h = H); autoplot(fedForecast1)
fedForecast2 <- forecast(fed_fit2, h = H); autoplot(fedForecast2)   # chosen model
fedForecast3 <- forecast(fed_fit3, h = H); autoplot(fedForecast3)

## Mortgage
autoForecastMortgage <- nolog_forecast(forecast(mort_auto,  h = H))
autoplot(autoForecastMortgage)
autoForecastMortgageB <- nolog_forecast(forecast(mort_autob, h = H))
autoplot(autoForecastMortgageB)

mortgageForecast1 <- nolog_forecast(forecast(mort_fit1, h = H)); autoplot(mortgageForecast1)
mortgageForecast2 <- nolog_forecast(forecast(mort_fit2, h = H)); autoplot(mortgageForecast2)
mortgageForecast3 <- nolog_forecast(forecast(mort_fit3, h = H)); autoplot(mortgageForecast3)  # chosen

## Housing
autoForecastHouses  <- nolog_forecast(forecast(houses_auto,  xreg = future_xreg))
autoplot(autoForecastHouses)

autoForecastHousesB <- nolog_forecast(forecast(houses_autob, xreg = future_xreg))
autoplot(autoForecastHousesB)

houseForecast1 <- nolog_forecast(forecast(houses_fit1, xreg = future_xreg)); autoplot(houseForecast1)
houseForecast2 <- nolog_forecast(forecast(houses_fit2, xreg = future_xreg)); autoplot(houseForecast2)
houseForecast3 <- nolog_forecast(forecast(houses_fit3, xreg = future_xreg)); autoplot(houseForecast3)

sHouseForecast1 <- nolog_forecast(forecast(houses_sarima1, xreg = future_xreg)); autoplot(sHouseForecast1)
sHouseForecast2 <- nolog_forecast(forecast(houses_sarima2, xreg = future_xreg)); autoplot(sHouseForecast2)  # chosen
sHouseForecast3 <- nolog_forecast(forecast(houses_sarima3, xreg = future_xreg_scaled)); autoplot(sHouseForecast3)
sHouseForecast4 <- nolog_forecast(forecast(houses_sarima4, xreg = future_xreg)); autoplot(sHouseForecast4)
sHouseForecast5 <- nolog_forecast(forecast(houses_sarima5, xreg = future_xreg)); autoplot(sHouseForecast5)

# prewhitening/ccf

# fed to mortgage

# raw ccf
ccf((fedfunds_ts),(log(mortgage_ts)))
# all lags significant, highly autocorr

# prewhiten - using fed_fit_autob because it has no MA terms
prewhiten(fedfunds_ts, log(mortgage_ts), x.model = fed_fit_autob)

# peak lag
k <- 1
n <- length(fedfunds_ts)

# output
y_align <- subset(fedfunds_ts, start = k+1)
# input
x_align <- subset(log(mortgage_ts), end = n-k)

# check fit
ccf_fit1 <- auto.arima(y_align, xreg = x_align)
coeftest(ccf_fit1)
sqrt(mean(residuals(ccf_fit1)^2))
Acf(residuals(ccf_fit1))
Box.test(residuals(ccf_fit1), lag = 24, type = "L")

# drop insignificant seasonal terms
ccf_fit2 <- Arima(y_align, order = c(2,1,1), xreg = x_align)
# all significant
coeftest(ccf_fit2)
# lower AIC
AIC(ccf_fit1, ccf_fit2)

# compare backtesting with and without xreg
ccf_base <- Arima(y_align, order = c(2,1,1))
backtest(ccf_base, y_align, orig=(0.8 * length(y_align)), h = 1)
backtest(ccf_fit2, y_align, orig=(0.8 * length(y_align)), h = 1, xre = x_align)

# backtesting comparable, if not a little worse with xreg. There is a lead relationship
# but it is not strong enough to help forecasting

# fed to housing starts:

# raw CCF
ccf(fedfunds_ts, log(houses_ts))
prewhiten(fedfunds_ts, log(houses_ts), x.model = fed_fit_autob)

# negative lag = fed funds lead housing
k <- 2
n <- length(houses_ts)
y_align2 <- subset(log(houses_ts), start = k + 1)
x_align2 <- subset(fedfunds_ts, end = n - k)

# transfer model + diagnostics
ccf_fit_FH1 <- auto.arima(y_align2, xreg = x_align2)
ccf_fit_FH1
coeftest(ccf_fit_FH1)
Box.test(residuals(ccf_fit_FH1), lag = 24, type = "L")

# not signifiant xreg

# mortgage to housing starts
fedfunds_ts <- window(fedfunds_ts, end = c(2025, 12))
houses_ts <- window(houses_ts, end = c(2025, 12))

ccf(log(mortgage_ts), log(houses_ts))
prewhiten(log(mortgage_ts), log(houses_ts), x.model = mort_fit1)

k <- 1
n <- length(houses_ts)
y_align3 <- subset(log(mortgage_ts), start = k + 1)
x_align3 <- subset(log(houses_ts), end = n - k)

ccf_fit_HM <- auto.arima(y_align3, xreg = x_align3)
coeftest(ccf_fit_HM)

base_hm <- Arima(y_align3, order = c(0,1,1))
backtest(base_hm,    y_align3, orig = (0.8 * length(y_align3)), h = 1)
backtest(ccf_fit_HM, y_align3, orig = (0.8 * length(y_align3)), h = 1, xre = x_align3)

# doesn't help, not a good regression
# Fed to mortgage: significant lead, but reversed (mortgage leads fed) and no out-of-sample value
# Fed to housing: no significant direct link
# Mortgage to housing: no clean significant lead
# Housing to mortgage: significant in-sample but worsens forecasting


#  VAR:

v <- cbind(fed = fedfunds_diff, mortgage = mortgage_ld, housing = houses_ld)

autoplot(v, facets = T)

# auto- and cross-correlation exploration
ccf(fedfunds_diff, mortgage_ld)
ccf(fedfunds_diff, houses_ld)
ccf(mortgage_ld, houses_ld)

# Order selection
s <- VARselect(v, lag.max = 12, type = "const")
s$selection

# lowest p=1
fit1 <- VAR(v, p = 1, type = "const")
serial.test(fit1, lags.pt = 24, type = "PT.asymptotic")
# doesnt pass

# p=2
fit2 <- VAR(v, p = 2, type = "const")
serial.test(fit2, lags.pt = 24, type = "PT.asymptotic")
# still doesn't pass

# highest p=3
fit3 <- VAR(v, p = 3, type = "const")
serial.test(fit3, lags.pt = 24, type = "PT.asymptotic")
# barely passes

# p=4
fit4 <- VAR(v, p = 4, type = "const")
serial.test(fit4, lags.pt = 24, type = "PT.asymptotic")

# actually passes
summary(fit4)
coeftest(fit4)
autoplot(forecast(fit4, h = 24))

############################################################
# Train/test split + overlay observed test values (his validation block)
############################################################
n <- nrow(v)
n_test <- 24

train <- window(v, end = time(v)[n - n_test])
test<- window(v, start = time(v)[n - n_test + 1])

fit_train <- VAR(train, p = 4, type = "const")
f <- forecast(fit_train, h = n_test)

# Granger - does each series help predict the others?
causality(fit4, cause = "fed")$Granger
causality(fit4, cause = "mortgage")$Granger
causality(fit4, cause = "housing")$Granger

plot(irf(fit4, impulse = "mortgage", response = c("fed", "housing"), n.ahead = 24, boot = TRUE))
