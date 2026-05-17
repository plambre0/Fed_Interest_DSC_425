library(ggplot2)
library(ggfortify)
library(tidyverse)
library(lubridate)
library(patchwork)
library(TSA)
library(lmtest)
library(forecast)
library(fUnitRoots)

fedfunds <- read.csv("FEDFUNDS.csv")
mortgage <- read.csv("MORTGAGE30US.csv")
houses <- read.csv("HOUST1F.csv")

head(fedfunds)
head(mortgage)
head(houses)

# convert to date type
mortgage$weeks = as.Date(mortgage$observation_date)

# convert to monthly with mean
mortgage = mortgage %>%
  mutate(date = floor_date(weeks, "month")) %>%
  group_by(date) %>%
  summarise(MORTGAGE30US = mean(MORTGAGE30US, na.rm = TRUE))

head(mortgage)


#### Initial Plots, differencing, logs, and stationarity checks ####

# ff ts object
fedfunds_ts = ts(fedfunds$FEDFUNDS, start=c(1976,1), frequency=12)
autoplot(fedfunds_ts) # not multiplicative, so not log, but mean in changing
# take diff
fedfunds_diff = diff(fedfunds_ts, differences = 1)
autoplot(fedfunds_diff)
# Stationary
adfTest(fedfunds_diff, type="nc")
kpss.test(fedfunds_diff)


# mortgage ts
mortgage_ts = ts(mortgage$MORTGAGE30US, start=c(1976,1), frequency=12)
autoplot(mortgage_ts)
# take log diff
mortgage_ld = diff(log(mortgage_ts), differences = 1)
autoplot(mortgage_ld)
# stationary
adfTest(mortgage_ld, type="nc")
kpss.test(mortgage_ld)


houses_ts = ts(houses$HOUST1F, start=c(1976,1), frequency=12)
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
fedfunds_plot
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

# fed_fit1 NOT good

fed_fit2 <- Arima(fedfunds_ts, order=c(2,1,1))
fed_fit2
coeftest(fed_fit2) # all significant, MA slightly less
checkresiduals(fed_fit2) # not good
Box.test(residuals(fed_fit2), lag=24, type="L") # residuals not white noise

# fed_fit2 NOT good

fed_fit3 <- Arima(fedfunds_ts, order=c(0,1,1))
fed_fit3
coeftest(fed_fit3) # significant
checkresiduals(fed_fit3) # not good
Box.test(residuals(fed_fit3), lag=24, type="L") # residuals not white noise

# fed_fit3 NOT good

fed_fit_auto <- auto.arima(fedfunds_ts)
fed_fit_auto
coeftest(fed_fit_auto) # significant
checkresiduals(fed_fit_auto) # still not great
Box.test(residuals(fed_fit_auto), lag=24, type="Ljung-Box") # residuals not white noise

fed_fit_autob <- auto.arima(fedfunds_ts, ic="bic")
fed_fit_autob
coeftest(fed_fit_autob) # significant
checkresiduals(fed_fit_autob) # still not great
Box.test(residuals(fed_fit_autob), lag=24, type="Ljung-Box")
# same as AIC

# AIC/BIC comparison
AIC(fed_fit1, fed_fit2, fed_fit3, fed_fit_auto, fed_fit_autob)
BIC(fed_fit1, fed_fit2, fed_fit3, fed_fit_auto, fed_fit_autob)

# auto.arima was the best fit by all metrics, but residuals did not reach white noise


# mortgage rates
mort_fit1 <- Arima(log(mortgage_ts), order=c(2,1,0))
mort_fit1
coeftest(mort_fit1) # significant
checkresiduals(mort_fit1) # pretty good
Box.test(residuals(mort_fit1), lag=24, type="Ljung-Box") # white noise

# mort_fit1 is good

mort_fit2 <- Arima(log(mortgage_ts), order=c(3,1,0))
mort_fit2
coeftest(mort_fit2) # AR3 less significant
checkresiduals(mort_fit2) # decent
Box.test(residuals(mort_fit2), lag=24, type="Ljung-Box") # white noise

# mort_fit2 is good


mort_fit3 <- Arima(log(mortgage_ts), order=c(0,1,1)) 
mort_fit3
coeftest(mort_fit3) #significant
checkresiduals(mort_fit3) # looks good
Box.test(residuals(mort_fit3), lag=24, type="Ljung-Box") # white noise

# mort_fit3 is good

mort_auto <- auto.arima(log(mortgage_ts))
mort_auto
coeftest(mort_auto) # also went with ARIMA(0,1,1)
checkresiduals(mort_auto)
Box.test(residuals(mort_auto), lag=24, type="Ljung-Box")

mort_autob <- auto.arima(log(mortgage_ts), ic="bic")
mort_autob
coeftest(mort_autob) # also went with ARIMA(0,1,1)
checkresiduals(mort_autob)
Box.test(residuals(mort_autob), lag=24, type="Ljung-Box")

# AIC/BIC comparison
AIC(mort_fit1, mort_fit2, mort_fit3, mort_auto, mort_autob)
BIC(mort_fit1, mort_fit2, mort_fit3, mort_auto, mort_autob)

# Auto.arima chose same model as mort_fit_3. Based on parsiomy and results, looks good


# housing
houses_fit1 <- Arima(log(houses_ts), order=c(0,1,1))
houses_fit1
coeftest(houses_fit1) # significant
checkresiduals(houses_fit1) # decent but misses seasonal structure
Box.test(residuals(houses_fit1), lag=24, type="Ljung-Box") # borderline

# good, could be better

houses_fit2 <- Arima(log(houses_ts), order=c(1,1,0))
houses_fit2
coeftest(houses_fit2) # significant
checkresiduals(houses_fit2) # misses seasonal structure
Box.test(residuals(houses_fit2), lag=24, type="Ljung-Box") # fails

# not white noise, NOT good

houses_fit3 <- Arima(log(houses_ts), order=c(1,1,1))
houses_fit3
coeftest(houses_fit3) # check significance
checkresiduals(houses_fit3) # still missing seasonal structure
Box.test(residuals(houses_fit3), lag=24, type="Ljung-Box") # borderline

# could be better

# add seasonal terms
houses_sarima1 <- Arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,0), period=12))
houses_sarima1
coeftest(houses_sarima1) # seasonal term not significant
checkresiduals(houses_sarima1) # improved with seasonal term
Box.test(residuals(houses_sarima1), lag=24, type="Ljung-Box")

# still not great

houses_sarima2 <- Arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima2
coeftest(houses_sarima2) # all terms significant
checkresiduals(houses_sarima2) # good
Box.test(residuals(houses_sarima2), lag=24, type="Ljung-Box")

# best yet

houses_sarima3 <- Arima(log(houses_ts), order=c(1,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima3
coeftest(houses_sarima3) # AR and MA both insignificant
checkresiduals(houses_sarima3)
Box.test(residuals(houses_sarima3), lag=24, type="Ljung-Box") # passes

houses_sarima4 <- Arima(log(houses_ts), order=c(2,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima4
coeftest(houses_sarima4) # ma1 marginally significant
checkresiduals(houses_sarima4) # good
Box.test(residuals(houses_sarima4), lag=24, type="Ljung-Box") # passes

# best residuals yet, less significant coeftest

houses_sarima5 <- Arima(log(houses_ts), order=c(1,1,2),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima5
coeftest(houses_sarima5) # check significance
checkresiduals(houses_sarima5) # good
Box.test(residuals(houses_sarima5), lag=24, type="Ljung-Box") # passes

# good

houses_auto <- auto.arima(log(houses_ts))
houses_auto
coeftest(houses_auto) # more complex model, a lot of insignificant terms
checkresiduals(houses_auto) # good
Box.test(residuals(houses_auto), lag=24, type="Ljung-Box") # passes

houses_autob <- auto.arima(log(houses_ts), ic="bic")
houses_autob
coeftest(houses_autob) # more complex model, a lot of insignificant terms
checkresiduals(houses_autob) # good
Box.test(residuals(houses_autob), lag=24, type="Ljung-Box") 

# best box test

# AIC/BIC comparison
AIC(houses_fit1, houses_fit2, houses_fit3,
    houses_sarima1, houses_sarima2, houses_sarima3,
    houses_sarima4, houses_sarima5, houses_auto, houses_autob)
BIC(houses_fit1, houses_fit2, houses_fit3,
    houses_sarima1, houses_sarima2, houses_sarima3,
    houses_sarima4, houses_sarima5, houses_auto, houses_autob)

# houses_sarima_2 the best by all metrics


##### Forecasting #####
detach("package:aTSA", unload=TRUE)

# auto.arima forecasts
autoForecastFed <- forecast(fed_fit_auto, h=24)
autoplot(autoForecastFed)

autoForecastFedB <- forecast(fed_fit_autob, h=24)
autoplot(autoForecastFedB)

autoForecastMortgage <- forecast(mort_auto, h=24)
autoplot(autoForecastMortgage)

autoForecastMortgageB <- forecast(mort_autob, h=24)
autoplot(autoForecastMortgageB)

autoForecastHouses <- forecast(houses_auto, h=24)
autoplot(autoForecastHouses)

autoForecastHousesB <- forecast(houses_autob, h=24)
autoplot(autoForecastHousesB)

# federal funds forecasts
fedForecast1 <- forecast(fed_fit1, h=24)
autoplot(fedForecast1)

fedForecast2 <- forecast(fed_fit2, h=24)
autoplot(fedForecast2)

fedForecast3 <- forecast(fed_fit3, h=24)
autoplot(fedForecast3)

# mortgage forecasts
mortgageForecast1 <- forecast(mort_fit1, h=24)
autoplot(mortgageForecast1)

mortgageForecast2 <- forecast(mort_fit2, h=24)
autoplot(mortgageForecast2)

mortgageForecast3 <- forecast(mort_fit3, h=24)
autoplot(mortgageForecast3)

# housing starts forecasts
houseForecast1 <- forecast(houses_fit1, h=24)
autoplot(houseForecast1)

houseForecast2 <- forecast(houses_fit2, h=24)
autoplot(houseForecast2)

houseForecast3 <- forecast(houses_fit3, h=24)
autoplot(houseForecast3)

# housing starts sarima forecasts
sHouseForecast1 <- forecast(houses_sarima1, h=24)
autoplot(sHouseForecast1)

sHouseForecast2 <- forecast(houses_sarima2, h=24)
autoplot(sHouseForecast2)

sHouseForecast3 <- forecast(houses_sarima3, h=24)
autoplot(sHouseForecast3)

sHouseForecast4 <- forecast(houses_sarima4, h=24)
autoplot(sHouseForecast4)

sHouseForecast5 <- forecast(houses_sarima5, h=24)
autoplot(sHouseForecast5)
