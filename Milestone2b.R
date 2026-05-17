library(ggplot2)
library(ggfortify)
library(tidyverse)
library(lubridate)
library(patchwork)
library(TSA)
library(lmtest)
library(forecast)
library(fUnitRoots)

fedfunds <- read.csv("FEDFUNDS-2.csv")
mortgage <- read.csv("MORTGAGE30US-2.csv")
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

# ff pacf
pacf(as.numeric(fedfunds_diff))
# strong spikes at lag1 and lag2
# AR2 process
# some more spikes later, but less significant
# possible MA

# mortgage pacf
pacf(as.numeric(mortgage_diff))
# AR2 or AR3

# houses pacf
pacf(as.numeric(houses_log_diff))
# AR1 or AR2
# maybe some MA
# signicant spikes at 12 and 24. seasonality?


fedfunds_pacf = autoplot(pacf(as.numeric(fedfunds_diff), lag.max=24, plot=FALSE)) +
  labs(title="Federal Funds Rate PACF", x="Lag", y="") +
  theme_minimal()
mortgage_pacf = autoplot(pacf(as.numeric(mortgage_ld), lag.max=24, plot=FALSE)) +
  labs(title="Mortgage Rate PACF", x="Lag", y="") +
  theme_minimal()
houses_pacf = autoplot(pacf(as.numeric(houses_ld), lag.max=24, plot=FALSE)) +
  labs(title="Housing Starts PACF", x="Lag", y="") +
  theme_minimal()

fedfunds_pacf / mortgage_pacf / houses_pacf



##### eacf #####
library(TSA)
TSA::eacf(mortgage_ts)
TSA::eacf(mortgage_diff)
TSA::eacf(fedfunds_ts)
TSA::eacf(fedfunds_diff)
TSA::eacf(houses_ts)
TSA::eacf(houses_log_diff)


eacf(as.numeric(fedfunds_diff))
eacf(as.numeric(mortgage_ld))
eacf(as.numeric(houses_ld))



###### arima ######

# federal funds
fed_fit1 <- arima(fedfunds_ts, order=c(2,1,0))
fed_fit1a = Arima(fedfunds_ts, order = c(2,1,0))
fed_fit2 <- arima(fedfunds_ts, order=c(2,1,1))
fed_fit2a <- Arima(fedfunds_ts, order=c(2,1,1))
fed_fit3 <- arima(fedfunds_ts, order=c(0,1,1))
fed_fit3a <- Arima(fedfunds_ts, order=c(0,1,1))
fed_fit_auto <- auto.arima(fedfunds_ts)
fed_fit_auto

# AIC comparison
AIC(fed_fit1, fed_fit2, fed_fit3, fed_fit_auto)

# diagnostics
checkresiduals(fed_fit_auto)
Box.test(residuals(fed_fit_auto), lag=24, type="Ljung-Box")
coeftest(fed_fit_auto)

# mortgage rates
mort_fit1 <- arima(log(mortgage_ts), order=c(2,1,0))
mort_fit1a <- Arima(log(mortgage_ts), order=c(2,1,0))
mort_fit2 <- arima(log(mortgage_ts), order=c(3,1,0))
mort_fit2a <- Arima(log(mortgage_ts), order=c(3,1,0))
mort_fit3 <- arima(log(mortgage_ts), order=c(0,1,1))
mort_fit3a <- Arima(log(mortgage_ts), order=c(0,1,1))
mort_auto <- auto.arima(log(mortgage_ts))
mort_auto

# AIC comparison
AIC(mort_fit1, mort_fit2, mort_fit3, mort_auto)

# diagnostics
checkresiduals(mort_auto)
Box.test(residuals(mort_auto), lag=24, type="Ljung-Box")
coeftest(mort_auto)

# housing starts
houses_fit1 <- arima(log(houses_ts), order=c(0,1,1))
houses_fit1a <- Arima(log(houses_ts), order=c(0,1,1))
houses_fit2 <- arima(log(houses_ts), order=c(1,1,0))
houses_fit2a <- Arima(log(houses_ts), order=c(1,1,0))
houses_fit3 <- arima(log(houses_ts), order=c(1,1,1))
houses_fit3a <- Arima(log(houses_ts), order=c(1,1,1))
houses_sarima1 <- arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,0), period=12))
houses_sarima1a <- Arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,0), period=12))
houses_sarima2 <- arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima2a <- Arima(log(houses_ts), order=c(0,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima3 <- arima(log(houses_ts), order=c(1,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima3a <- Arima(log(houses_ts), order=c(1,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima4 <- arima(log(houses_ts), order=c(2,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima4a <- Arima(log(houses_ts), order=c(2,1,1),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima5 <- arima(log(houses_ts), order=c(1,1,2),
                        seasonal=list(order=c(1,0,1), period=12))
houses_sarima5a <- Arima(log(houses_ts), order=c(1,1,2),
                        seasonal=list(order=c(1,0,1), period=12))
houses_auto <- auto.arima(log(houses_ts))
houses_auto

# AIC comparison
AIC(houses_fit1, houses_fit2, houses_fit3,
    houses_sarima1, houses_sarima2, houses_sarima3,
    houses_sarima4, houses_sarima5, houses_auto)

# Ljung-Box for all
cat("\n--- Housing Starts Ljung-Box p-values ---\n")
cat("ARIMA(0,1,1):", Box.test(residuals(houses_fit1), lag=24, type="Ljung-Box")$p.value, "\n")
cat("ARIMA(1,1,0):", Box.test(residuals(houses_fit2), lag=24, type="Ljung-Box")$p.value, "\n")
cat("ARIMA(1,1,1):", Box.test(residuals(houses_fit3), lag=24, type="Ljung-Box")$p.value, "\n")
cat("ARIMA(0,1,1)(1,0,0)[12]:", Box.test(residuals(houses_sarima1), lag=24, type="Ljung-Box")$p.value, "\n")
cat("ARIMA(0,1,1)(1,0,1)[12]:", Box.test(residuals(houses_sarima2), lag=24, type="Ljung-Box")$p.value, "\n")
cat("ARIMA(1,1,1)(1,0,1)[12]:", Box.test(residuals(houses_sarima3), lag=24, type="Ljung-Box")$p.value, "\n")
cat("ARIMA(2,1,1)(1,0,1)[12]:", Box.test(residuals(houses_sarima4), lag=24, type="Ljung-Box")$p.value, "\n")
cat("ARIMA(1,1,2)(1,0,1)[12]:", Box.test(residuals(houses_sarima5), lag=24, type="Ljung-Box")$p.value, "\n")
cat("Auto ARIMA:", Box.test(residuals(houses_auto), lag=24, type="Ljung-Box")$p.value, "\n")

# confirm sarima2 is preferred over sarima4
coeftest(houses_sarima4)

# best model diagnostics
checkresiduals(houses_sarima2)
Box.test(residuals(houses_sarima2), lag=24, type="Ljung-Box")
coeftest(houses_sarima2)

#forecast models
autoForecastFed = forecast(fed_fit_auto, h=24)
autoplot(autoForecastFed)
autoForecastMortgage = forecast(mort_auto, h=24)
autoplot(autoForecastMortgage)
autoForecastHouses = forecast(houses_auto, h=24)
autoplot(autoForecastHouses)
fedForecast1 = forecast(fed_fit1a, 24)
autoplot(fedForecast1)
fedForecast2 = forecast(fed_fit2a, 24)
autoplot(fedForecast2)
fedForecast3 = forecast(fed_fit3a, 24)
autoplot(fedForecast3)
mortgageForecast1 = forecast(mort_fit1a, 24)
autoplot(mortgageForecast1)
mortgageForecast2 = forecast(mort_fit2a, 24)
autoplot(mortgageForecast2)
mortgageForecast3 = forecast(mort_fit3a, 24)
autoplot(mortgageForecast3)
houseForecast1 = forecast(houses_fit1a, 24)
autoplot(houseForecast1)
houseForecast2 = forecast(houses_fit2a, 24)
autoplot(houseForecast2)
houseForecast3 = forecast(houses_fit3a, 24)
autoplot(houseForecast3)
sHouseForecast1 = forecast(houses_sarima1a, 24)
autoplot(sHouseForecast1)
sHouseForecast2 = forecast(houses_sarima2a, 24)
autoplot(sHouseForecast2)
sHouseForecast3 = forecast(houses_sarima3a, 24)
autoplot(sHouseForecast3)
sHouseForecast4 = forecast(houses_sarima4a, 24)
autoplot(sHouseForecast4)
sHouseForecast5 = forecast(houses_sarima5a, 24)
autoplot(sHouseForecast5)
