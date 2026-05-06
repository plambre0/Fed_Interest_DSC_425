library(ggplot2)
library(ggfortify)
library(tidyverse)
library(lubridate)

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

# ts objects
fedfunds_ts = ts(fedfunds$FEDFUNDS, start=c(1976,1), frequency=12)
autoplot(fedfunds_ts)
# take diff
fedfunds_diff = diff(fedfunds_ts, differences = 1)
# high volatility
autoplot(fedfunds_diff)

mortgage_ts = ts(mortgage$MORTGAGE30US, start=c(1976,1), frequency=12)
autoplot(mortgage_ts)
# take diff
mortgage_diff = diff(mortgage_ts, differences = 1)
autoplot(mortgage_diff)

houses_ts = ts(houses$HOUST1F, start=c(1976,1), frequency=12)
autoplot(houses_ts)
# log diff
houses_log_diff = diff(log(houses_ts), differences = 1)
autoplot(houses_log_diff)


# acf

# eacf
library(TSA)
TSA::eacf(mortgage_ts)
TSA::eacf(mortgage_diff)
TSA::eacf(fedfunds_ts)
TSA::eacf(fedfunds_diff)
TSA::eacf(houses_ts)
TSA::eacf(houses_log_diff)
# pacf
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

# arima
