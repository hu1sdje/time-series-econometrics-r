# Econometric Analysis of Philippine Quarterly GDP (Q1 2000 – Q1 2026)

library(here)
library(readxl)
library(forecast)
library(rugarch)
library(FinTS)
library(urca)

# Loading data and forming a time series
# Data: Philippine quarterly GDP and household consumption.
filippines <- read_excel(here("philippines_gdp.xlsx"))

gdp_vec <- as.numeric(filippines$GDP)
gdp_vec <- gdp_vec[!is.na(gdp_vec)]
n <- length(gdp_vec)

y <- ts(gdp_vec, start = c(2000, 1), frequency = 4)

plot(y, main = "Philippine GDP, quarterly (Q1 2000 - Q1 2026)",
     ylab = "GDP (mln pesos, 2018 prices)", xlab = "Year", col = "darkblue")

# Trend and seasonality
t  <- 1:n
t2 <- t^2

mod_lin <- lm(y ~ t)
mod_qua <- lm(y ~ t + t2)
mod_exp <- lm(log(y) ~ t) 

cat("Linear trend: R2_adj =", summary(mod_lin)$adj.r.squared,
    " AIC =", AIC(mod_lin), "\n")
cat("Quadratic trend: R2_adj =", summary(mod_qua)$adj.r.squared,
    " AIC =", AIC(mod_qua), "\n")
cat("Exponential trend: R2_adj =", summary(mod_exp)$adj.r.squared,
    " AIC =", AIC(mod_exp), "\n")
# Exponential trend has the highest R2_adj and lowest AIC -> selected.

# Seasonality via dummy variables (reference quarter: Q4)
season <- factor(cycle(y), labels = c("Q1", "Q2", "Q3", "Q4"))
season <- relevel(season, ref = "Q4")

# Log-linear model: trend + seasonal dummies
mod_ts <- lm(log(y) ~ t + season)
print(summary(mod_ts))
print(anova(lm(log(y) ~ t), mod_ts))   # significance of adding seasonality

resid_ts <- ts(residuals(mod_ts), start = c(2000, 1), frequency = 4)

# Fitted model vs. actual series
y_fit <- exp(fitted(mod_ts))
plot(y, main = "Actual series vs. trend-seasonal model",
     ylab = "GDP", col = "black", lwd = 2)
lines(ts(y_fit, start = c(2000, 1), frequency = 4), col = "red", lwd = 2)
legend("topleft", c("Actual", "Model"), col = c("black", "red"), lwd = 2)

# Interval forecast 12 steps ahead (3 years)
h <- 12
t_fc <- (n + 1):(n + h)
# Continue the quarterly cycle correctly: after Q1 2026 comes Q2, Q3, Q4, Q1...
last_q <- cycle(y)[n]
seq_q  <- (last_q:(last_q + h - 1)) %% 4 + 1
season_fc <- factor(c("Q1", "Q2", "Q3", "Q4")[seq_q], levels = levels(season))

new_dat  <- data.frame(t = t_fc, season = season_fc)
pred_log <- predict(mod_ts, newdata = new_dat,
                    interval = "prediction", level = 0.95)
pred <- exp(pred_log)                   # back-transform log -> level
print(round(pred, 0))

plot(y, xlim = c(2000, 2030), ylim = c(min(y), max(pred[, "upr"]) * 1.05),
     main = "Interval GDP forecast (trend-seasonal model)",
     ylab = "GDP", lwd = 2)
fc_time <- as.numeric(time(ts(rep(NA, h), start = c(2026, 2), frequency = 4)))
lines(fc_time, pred[, "fit"], col = "red",  lwd = 2)
lines(fc_time, pred[, "lwr"], col = "blue", lty = 2)
lines(fc_time, pred[, "upr"], col = "blue", lty = 2)
legend("topleft", c("Actual", "Forecast", "95% CI"),
       col = c("black", "red", "blue"), lwd = c(2, 2, 1), lty = c(1, 1, 2))

# Testing residuals for trend (Spearman's test)
plot(resid_ts, main = "Model residuals (trend + seasonality)",
     ylab = "Residuals", xlab = "Year", col = "darkblue")

sp <- cor.test(as.numeric(resid_ts), 1:n, method = "spearman", exact = FALSE)
print(sp)   # p-value > 0.05 => no trend left in residuals

# 4. Smoothing (moving averages and Holt-Winters)
# Moving averages
m <- 4    # window = 1 year
p <- 8    # window = 2 years

ma_m <- stats::filter(y, rep(1 / m, m), sides = 2)
ma_p <- stats::filter(y, rep(1 / p, p), sides = 2)

plot(y, main = "Moving averages (windows of 4 and 8 quarters)",
     ylab = "GDP", lwd = 1.5)
lines(ma_m, col = "red",  lwd = 2)
lines(ma_p, col = "blue", lwd = 2)
legend("topleft", c("Actual", "MA(4)", "MA(8)"),
       col = c("black", "red", "blue"), lwd = c(1.5, 2, 2))

# Holt-Winters exponential smoothing (with seasonality)
hw <- HoltWinters(y)
plot(y, main = "Holt-Winters exponential smoothing",
     ylab = "GDP", xlab = "Year", col = "black", lwd = 1.5)
lines(fitted(hw)[, 1], col = "red", lwd = 2)
legend("topleft", c("Actual", "Holt-Winters"), col = c("black", "red"), lwd = 2)

# 5. Harmonic and spectral analysis
# Residuals after removing trend only (to detect seasonal frequencies)
mod_trend_only <- lm(log(y) ~ t)
e <- as.numeric(residuals(mod_trend_only))
N <- length(e)

# Fourier coefficients
fourier_coefs <- function(x) {
  N <- length(x); K <- floor(N / 2); k <- 1:K
  a <- sapply(k, function(kk) (2 / N) * sum(x * cos(2 * pi * kk * (1:N) / N)))
  b <- sapply(k, function(kk) (2 / N) * sum(x * sin(2 * pi * kk * (1:N) / N)))
  data.frame(k = k, freq = k / N, period = N / k,
             a = a, b = b, amplitude = sqrt(a^2 + b^2))
}
FC <- fourier_coefs(e)
print(head(FC[order(-FC$amplitude), ], 10))

# Periodogram
spec.pgram(e, log = "no", taper = 0, fast = FALSE, detrend = FALSE,
           demean = TRUE, main = "Periodogram of residuals")
abline(v = 1 / 4, col = "red", lty = 2)   # period of 4 quarters
abline(v = 1 / 2, col = "red", lty = 2)   # period of 2 quarters

# Model: trend + leading harmonics
K_sel <- 3
top_k <- FC[order(-FC$amplitude), ][1:K_sel, ]
harm_part <- rep(0, N)
for (j in 1:K_sel) {
  harm_part <- harm_part +
    top_k$a[j] * cos(2 * pi * top_k$k[j] * (1:N) / N) +
    top_k$b[j] * sin(2 * pi * top_k$k[j] * (1:N) / N)
}
y_model <- exp(fitted(mod_trend_only) + harm_part)
plot(y, main = "Actual vs. model (trend + harmonics)", ylab = "GDP", lwd = 2)
lines(ts(y_model, start = c(2000, 1), frequency = 4), col = "red", lwd = 2)

# Seasonal harmonics (periods of ~4 and ~2 quarters)
season_k <- FC[round(FC$period) %in% c(2, 4), ]
print(season_k)

harm <- rep(0, N)
for (j in 1:nrow(season_k)) {
  harm <- harm +
    season_k$a[j] * cos(2 * pi * season_k$k[j] * (1:N) / N) +
    season_k$b[j] * sin(2 * pi * season_k$k[j] * (1:N) / N)
}
plot(e, type = "l", col = "gray50", main = "Seasonal harmonics",
     ylab = "Residuals", xlab = "Quarter")
lines(harm, col = "red", lwd = 2)

# 6. Modeling residuals: ARMA / ARIMA
r <- as.numeric(resid_ts)
sp2 <- cor.test(r, 1:length(r), method = "spearman", exact = FALSE)
cat("Spearman's test on ARMA residuals: p-value =", sp2$p.value, "\n")

# Hold out the last 10 observations to evaluate the forecast
H <- 10
r_train <- r[1:(length(r) - H)]
r_test  <- r[(length(r) - H + 1):length(r)]

# Grid of ARMA(p, q) specifications
m_arma11 <- Arima(r_train, order = c(1, 0, 1))
m_arma22 <- Arima(r_train, order = c(2, 0, 2))
m_arma21 <- Arima(r_train, order = c(2, 0, 1))
m_arma23 <- Arima(r_train, order = c(2, 0, 3))
m_arma24 <- Arima(r_train, order = c(2, 0, 4))

cat("ARMA(1,1): AIC =", AIC(m_arma11), "\n")
cat("ARMA(2,2): AIC =", AIC(m_arma22), "\n")
cat("ARMA(2,1): AIC =", AIC(m_arma21), "\n")
cat("ARMA(2,3): AIC =", AIC(m_arma23), "\n")
cat("ARMA(2,4): AIC =", AIC(m_arma24), "\n")

# Automatic selection as a cross-check of the manual choice
m_auto <- auto.arima(r_train, d = 0, seasonal = FALSE, stationary = TRUE)
print(m_auto)

# Best by AIC: ARMA(1,1) and ARMA(2,3)
# Interval forecast H = 10 steps
fc1 <- forecast(m_arma23, h = H, level = 95)
fc2 <- forecast(m_arma11, h = H, level = 95)
print(fc1)
print(fc2)

# Compare forecasts with actual values
plot(fc1, main = "Forecast ARMA(2,3) vs. actual",
     xlim = c(length(r_train) - 30, length(r)),
     ylim = c(-0.15, 0.15))
lines((length(r_train) + 1):length(r), r_test, col = "red", lwd = 2, type = "o")
legend("topleft", c("Forecast", "95% CI", "Actual"),
       col = c("blue", "grey", "red"), lty = 1, lwd = c(2, 5, 2))

plot(fc2, main = "Forecast ARMA(1,1) vs. actual",
     xlim = c(length(r_train) - 30, length(r)))
lines((length(r_train) + 1):length(r), r_test, col = "red", lwd = 2, type = "o")
legend("topleft", c("Forecast", "95% CI", "Actual"),
       col = c("blue", "grey", "red"), lty = 1, lwd = c(2, 5, 2))

# Volatility modeling: ARMA-GARCH
print(ArchTest(residuals(m_arma11), lags = 10))

# ARMA(1,1) + GARCH(1,1) specification with normal errors
spec <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "norm"
)

fit_gh <- ugarchfit(spec, data = r_train, solver = "hybrid")
print(fit_gh)

# Ex post forecast 10 steps ahead (on residuals)
fc_gh   <- ugarchforecast(fit_gh, n.ahead = H)
mean_fc <- as.numeric(fitted(fc_gh))
sd_fc   <- as.numeric(sigma(fc_gh))

# 95% confidence interval
lwr_gh <- mean_fc - 1.96 * sd_fc
upr_gh <- mean_fc + 1.96 * sd_fc

# Move from the residual forecast back to the GDP level:
# log(GDP_hat) = trend + seasonality (from mod_ts) + ARMA-GARCH residual forecast
t_h      <- (n - H + 1):n
seas_h   <- season[t_h]
det_part <- predict(mod_ts, newdata = data.frame(t = t_h, season = seas_h))

y_fc <- exp(det_part + mean_fc)
y_lo <- exp(det_part + lwr_gh)
y_up <- exp(det_part + upr_gh)
y_actual_tail <- y[t_h]

# Ex post forecast chart
plot(t_h, y_actual_tail, type = "o", col = "black", lwd = 2,
     ylim = range(c(y_actual_tail, y_lo, y_up)),
     xlab = "Observation index", ylab = "GDP",
     main = "Ex post GDP forecast: trend + ARMA(1,1) + GARCH(1,1)")
lines(t_h, y_fc, col = "red",  lwd = 2)
lines(t_h, y_lo, col = "blue", lty = 2)
lines(t_h, y_up, col = "blue", lty = 2)
legend("topleft", c("Actual", "Forecast", "95% CI"),
       col = c("black", "red", "blue"), lwd = c(2, 2, 1), lty = c(1, 1, 2))

# Ex post forecast quality
cat("MAE  =", mean(abs(y_actual_tail - y_fc)), "\n")
cat("RMSE =", sqrt(mean((y_actual_tail - y_fc)^2)), "\n")
cat("MAPE =", 100 * mean(abs((y_actual_tail - y_fc) / y_actual_tail)), "%\n")
cat("Coverage of the 95% CI:",
    mean(y_actual_tail >= y_lo & y_actual_tail <= y_up), "\n")


# 8. Cointegration: Engle-Granger and ECM
# Two series that move together: GDP and household consumption.
# Consumption is ~70-75% of GDP by the expenditure approach => economic link.
x1 <- log(as.numeric(filippines$GDP))
x2 <- log(as.numeric(filippines$`household consumption`))
x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
n2 <- min(length(x1), length(x2))
x1 <- x1[1:n2]; x2 <- x2[1:n2]

X1 <- ts(x1, start = c(2000, 1), frequency = 4)
X2 <- ts(x2, start = c(2000, 1), frequency = 4)

plot(cbind(X1, X2), main = "log(GDP) and log(household consumption)")

# Order of integration via the augmented Dickey-Fuller test
adf_test_summary <- function(series, name) {
  cat("\n--- ADF for ", name, " (levels) ---\n", sep = "")
  print(summary(ur.df(series, type = "trend", selectlags = "AIC")))
  cat("\n--- ADF for diff ", name, " (first differences) ---\n", sep = "")
  print(summary(ur.df(diff(series), type = "drift", selectlags = "AIC")))
}
adf_test_summary(X1, "log_GDP")
adf_test_summary(X2, "log_HC")

# Engle-Granger cointegration test
# Step 1: cointegrating regression
coint_reg <- lm(X1 ~ X2)
print(summary(coint_reg))
u_hat <- residuals(coint_reg)

# Step 2: ADF on the residuals of the cointegrating regression
eg_adf <- ur.df(u_hat, type = "none", selectlags = "AIC")
print(summary(eg_adf))

# Error correction model (ECM)
dX1   <- diff(X1)
dX2   <- diff(X2)
u_lag <- u_hat[-length(u_hat)]   # u_{t-1}

# ECM without seasonality
ecm1 <- lm(dX1 ~ dX2 + u_lag)
print(summary(ecm1))
plot(residuals(ecm1), type = "l", main = "ECM residuals", ylab = "resid")
abline(h = 0, col = "red")
Acf(residuals(ecm1), main = "ACF of ECM residuals")

# ECM with seasonality
season_ecm <- factor(cycle(X1)[-1], labels = c("Q1", "Q2", "Q3", "Q4"))
season_ecm <- relevel(season_ecm, ref = "Q4")

ecm_season <- lm(dX1 ~ dX2 + u_lag + season_ecm)
print(summary(ecm_season))
plot(residuals(ecm_season), type = "l",
     main = "ECM residuals (with seasonality)", ylab = "resid")
abline(h = 0, col = "red")
Acf(residuals(ecm_season), main = "ACF of ECM residuals (with seasonality)")

