library(zoo)
library(forecast)
library(tseries)
library(TSA)

# Just set your dir here.
setwd("C:/Users/Paolo/Desktop/DSC425_Project")

################################################################################
# Load FEDFUNDS
################################################################################
FEDFUNDS <- read.csv("FEDFUNDS.csv")
FEDFUNDS_ts <- window(ts(FEDFUNDS$FEDFUNDS, start = c(1976,1), frequency = 12), c(1976, 1), c(2026, 1))
FEDFUNDS_ts_diff = diff(FEDFUNDS_ts, differences = 1)
################################################################################
# Load Regressors
################################################################################
# Personal Consumption Expenditures (PCE)
PCE <- read.csv("PCE.csv")
PCE_ts <- ts(PCE$PCE, start = c(1959, 1), frequency = 12)
# 100*(Real Gross Domestic Product-Real Potential Gross Domestic Product)/Real Potential Gross Domestic Product) (GDPC1_GDPPOT)
GDPC1_GDPPOT <- read.csv("GDPC1_GDPPOT.csv")
GDPC1_GDPPOT_ts <- ts(GDPC1_GDPPOT$GDPC1_GDPPOT, start = c(1949, 1), frequency = 4)
# Moody's Seasoned Baa Corporate Bond Yield Relative to Yield on 10-Year Treasury Constant Maturity (BAA10Y)
BAA10Y <- readxl::read_xlsx("BAA10Y.xlsx", sheet = "Daily")
BAA10Y_ts <- ts(BAA10Y$BAA10Y, start = c(1986, 1), frequency = 365)
# Loans and Leases in Bank Credit, All Commercial Banks (LOANS)
LOANS <- read.csv("LOANS.csv")
LOANS_ts <- ts(LOANS$LOANS, start = c(1947, 1), frequency = 12)
# Market Yield on U.S. Treasury Securities at 10-Year Constant Maturity, Quoted on an Investment Basis (DGS10)
DGS10 <- readxl::read_xlsx("DGS10.xlsx", sheet = "Daily")
DGS10_ts <- ts(DGS10$DGS10, start = c(1962, 1), frequency = 365)
# St. Louis Fed Financial Stress Index (STLFSI4)
STLFSI4 <- readxl::read_xlsx("STLFSI4.xlsx", sheet = "Weekly, Ending Friday")
STLFSI4_ts <- ts(STLFSI4$STLFSI4, start = c(1993, 12), frequency = 52)
# CBOE Volatility Index: VIX (VIXCLS)
VIXCLS <- readxl::read_xlsx("VIXCLS.xlsx", sheet = "Daily, Close")
VIXCLS_ts <- ts(VIXCLS$VIXCLS, start = c(1990, 1), frequency = 365)
################################################################################
# Exploratory Analysis on Regressors
################################################################################
autoplot(PCE_ts)
forecast::BoxCox.lambda(PCE_ts)
PCE_ts_ln <- log(PCE_ts)
autoplot(PCE_ts_ln)
adf.test(PCE_ts_ln)
kpss.test(PCE_ts_ln)
PCE_ts_ln_diff <- diff(PCE_ts_ln)
autoplot(PCE_ts_ln_diff)
t.test(PCE_ts_ln_diff)

autoplot(GDPC1_GDPPOT_ts)
forecast::BoxCox.lambda(GDPC1_GDPPOT_ts)
adf.test(GDPC1_GDPPOT_ts)
kpss.test(GDPC1_GDPPOT_ts)
GDPC1_GDPPOT_ts_diff <- diff(GDPC1_GDPPOT_ts)
autoplot(GDPC1_GDPPOT_ts_diff)
t.test(GDPC1_GDPPOT_ts_diff)

autoplot(BAA10Y_ts)
forecast::BoxCox.lambda(BAA10Y_ts)
BAA10Y_ts_ln <- log(BAA10Y_ts)
BAA10Y_ts_ln_diff <- diff(BAA10Y_ts_ln)
adf.test(na.remove(BAA10Y_ts_ln_diff))
kpss.test(na.remove(BAA10Y_ts_ln_diff))
autoplot(BAA10Y_ts_ln_diff)
t.test(BAA10Y_ts_ln_diff)

autoplot(LOANS_ts)
forecast::BoxCox.lambda(LOANS_ts)
LOANS_ts_ln <- log(LOANS_ts)
adf.test(LOANS_ts_ln)
kpss.test(LOANS_ts_ln)
LOANS_ts_ln_diff <- diff(LOANS_ts_ln)
autoplot(LOANS_ts_ln_diff)
t.test(LOANS_ts_ln_diff)

autoplot(DGS10_ts)
forecast::BoxCox.lambda(DGS10_ts)
adf.test(na.remove(DGS10_ts))
kpss.test(na.remove(DGS10_ts))
DGS10_ts_diff <- diff(DGS10_ts)
autoplot(DGS10_ts_diff)
t.test(DGS10_ts_diff)

autoplot(STLFSI4_ts)
forecast::BoxCox.lambda(STLFSI4_ts)
adf.test(na.remove(STLFSI4_ts))
kpss.test(na.remove(STLFSI4_ts))
STLFSI4_ts_diff <- diff(STLFSI4_ts)
autoplot(STLFSI4_ts_diff)
t.test(STLFSI4_ts_diff)

autoplot(VIXCLS_ts)
forecast::BoxCox.lambda(VIXCLS_ts)
LOANS_ts_ln <- log(VIXCLS_ts)
adf.test(na.remove(VIXCLS_ts))
kpss.test(na.remove(VIXCLS_ts))
VIXCLS_ts_ln_diff <- diff(VIXCLS_ts)
autoplot(VIXCLS_ts_ln_diff)
t.test(VIXCLS_ts_ln_diff)

################################################################################
# Transform Regressors Time Scale
################################################################################
library(tempdisagg)
library(tsbox)
library(imputeTS)
library(astsa)

DGS10_ts_mnth <- aggregate(as.zoo(DGS10_ts), as.yearmon, function(x) x[1])
DGS10_ts_mnth_imp <- na_seadec(DGS10_ts_mnth, algorithm = "interpolation")
GDPC1_GDPPOT_ts_wind <- window(GDPC1_GDPPOT_ts,
                               start = max(range(time(GDPC1_GDPPOT_ts))[1], range(time(DGS10_ts_mnth_imp))[1]),
                               end   = min(range(time(GDPC1_GDPPOT_ts))[2], range(time(DGS10_ts_mnth_imp))[2]))
DGS10_ts_mnth_imp_wind <- window(DGS10_ts_mnth_imp,
                                 start = max(range(time(GDPC1_GDPPOT_ts))[1], range(time(DGS10_ts_mnth_imp))[1]),
                                 end   = min(range(time(GDPC1_GDPPOT_ts))[2], range(time(DGS10_ts_mnth_imp))[2]))
GDPC1_GDPPOT_ts_disag <- td(GDPC1_GDPPOT_ts_wind ~ DGS10_ts_mnth_imp_wind,
                            method = "litterman-fixed",
                            fixed.rho = 0.99)
GDPC1_GDPPOT_ts_disag_fit <- predict(GDPC1_GDPPOT_ts_disag)

BAA10Y_ts_mnth <- aggregate(as.zoo(BAA10Y_ts), as.yearmon, function(x) x[1])
BAA10Y_ts_mnth_imp <- na_seadec(BAA10Y_ts_mnth, algorithm = "interpolation")

LOANS_ts_mnth <- aggregate(as.zoo(LOANS_ts), as.yearmon, function(x) x[1])
LOANS_ts_mnth_imp <- na_seadec(LOANS_ts_mnth, algorithm = "interpolation")

STLFSI4_ts_mnth <- aggregate(as.zoo(STLFSI4_ts), as.yearmon, function(x) x[1])

VIXCLS_ts_mnth <- aggregate(as.zoo(VIXCLS_ts), as.yearmon, function(x) x[1])
VIXCLS_ts_mnth_imp <- na_seadec(VIXCLS_ts_mnth, algorithm = "interpolation")

PCE_ts_ln
GDPC1_GDPPOT_ts_disag_fit
BAA10Y_ts_mnth_imp
LOANS_ts_ln_mnth_imp <- log(LOANS_ts_mnth_imp)
STLFSI4_ts_mnth
DGS10_ts_mnth_imp
VIXCLS_ts_ln_mnth_imp <- log(VIXCLS_ts_mnth_imp)

################################################################################
FEDFUNDS_ts <- as.ts(FEDFUNDS_ts)
PCE_ts_ln <- as.ts(PCE_ts_ln)
GDPC1_GDPPOT_ts_disag_fit <- as.ts(GDPC1_GDPPOT_ts_disag_fit)
BAA10Y_ts_mnth_imp <- as.ts(BAA10Y_ts_mnth_imp)
STLFSI4_ts_mnth <- as.ts(STLFSI4_ts_mnth)
DGS10_ts_mnth_imp <- as.ts(DGS10_ts_mnth_imp)
VIXCLS_ts_ln_mnth_imp <- as.ts(VIXCLS_ts_ln_mnth_imp)

series_list <- list(
  FEDFUNDS_ts,
  PCE_ts_ln,
  GDPC1_GDPPOT_ts_disag_fit,
  BAA10Y_ts_mnth_imp,
  STLFSI4_ts_mnth,
  DGS10_ts_mnth_imp,
  VIXCLS_ts_ln_mnth_imp)

common_start <- do.call(
  max,
  lapply(series_list, function(x) start(x)[1] + (start(x)[2]-1)/frequency(x)))

common_end <- do.call(
  min,
  lapply(series_list, function(x) end(x)[1] + (end(x)[2]-1)/frequency(x)))

series_common <- lapply(
  series_list,
  function(x) window(x, start = common_start, end = common_end))

FEDFUNDS_ts_common            <- series_common[[1]]
PCE_ts_ln_common             <- series_common[[2]]
GDPC1_GDPPOT_ts_common       <- series_common[[3]]
BAA10Y_ts_mnth_imp_common    <- series_common[[4]]
STLFSI4_ts_mnth_common       <- series_common[[5]]
DGS10_ts_mnth_imp_common     <- series_common[[6]]
VIXCLS_ts_ln_mnth_imp_common <- series_common[[7]]

common_data <- na.omit(ts.intersect(
  FEDFUNDS_ts,
  PCE_ts_ln,
  GDPC1_GDPPOT_ts_disag_fit,
  BAA10Y_ts_mnth_imp,
  STLFSI4_ts_mnth,
  DGS10_ts_mnth_imp,
  VIXCLS_ts_ln_mnth_imp
))
################################################################################
library(vars)
common_data_varselect <- as.data.frame(t(VARselect(common_data)$criteria))
common_data_varselect$AIC_rank <- rank(common_data_varselect[,1])
common_data_varselect$HQ_rank <- rank(common_data_varselect[,2])
common_data_varselect$SC_rank <- rank(common_data_varselect[,3])
common_data_varselect$FPE_rank <- rank(common_data_varselect[,4])
which.min(common_data_varselect$AIC_rank)
which.min(common_data_varselect$HQ_rank)
which.min(common_data_varselect$SC_rank)
which.min(common_data_varselect$FPE_rank)
plot(common_data_varselect[,1])
plot(common_data_varselect[,2])
plot(common_data_varselect[,3])
plot(common_data_varselect[,4])

common_data_VAR_4 <- VAR(common_data, p = 4)
common_data_jo <- ca.jo(common_data, spec = "longrun")
summary(common_data_jo)

library(tsDyn)
library(bvartools)

seas <- gen_vec(data = common_data, p = 4, r = 1, const = "unrestricted", seasonal = "unrestricted")
seas <- seas$data$X[, 22:33]

est_tsdyn <- VECM(common_data, lag = 4, r = 1, include = "none", estim = "ML")
summary(est_tsdyn)

causality(common_data_VAR_4, cause = c("FEDFUNDS_ts"))
causality(common_data_VAR_4, cause = c("PCE_ts_ln",
                                       "GDPC1_GDPPOT_ts_disag_fit",
                                       "BAA10Y_ts_mnth_imp",
                                       "STLFSI4_ts_mnth",
                                       "DGS10_ts_mnth_imp",
                                       "VIXCLS_ts_ln_mnth_imp"))

order <- c(
  "PCE_ts_ln",
  "GDPC1_GDPPOT_ts_disag_fit",
  "BAA10Y_ts_mnth_imp",
  "STLFSI4_ts_mnth",
  "VIXCLS_ts_ln_mnth_imp",
  "DGS10_ts_mnth_imp",
  "FEDFUNDS_ts"
)

svar_model <- VECM(common_data[, order], lag = 1)
svar_irf <- vars::irf(svar_model,
                      impulse = "PCE_ts_ln",
                      response = "FEDFUNDS_ts",
                      n.ahead = 24,
                      boot = TRUE)
