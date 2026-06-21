# time-series-econometrics-r
End-to-end time-series econometrics in R: trend, seasonality, ARMA/ARIMA with train/test validation, ARMA-GARCH volatility, and cointegration (ADF, Engle-Granger, ECM) on Philippine GDP.

## Overview
A complete econometric study of quarterly Philippine GDP 208 constant prices (Q1 2000 – Q1 2026). The project moves from descriptive trend and seasonality modelling through spectral analysis and ARMA/ARIMA forecasting to volatility modelling with ARMA-GARCH, and finishes with a cointegration analysis of the long-run link between GDP and household consumption.

The full analysis with all charts and output is in the rendered report: 

## Methods
- Trend modelling (linear, quadratic, exponential), selected by adjusted R² and AIC.
- Seasonality via quarterly dummy variables, tested with ANOVA.
- Residual diagnostics with Spearman's rank correlation.
- Smoothing: moving averages and Holt-Winters exponential smoothing.
- Harmonic and spectral analysis: Fourier coefficients and periodogram.
- ARMA/ARIMA model selection with a train/test split and AIC, cross-checked against `auto.arima`.
- Volatility modelling: ARCH LM test and an ARMA(1,1)-GARCH(1,1) model, evaluated out of sample.
- Cointegration: augmented Dickey-Fuller test, Engle-Granger two-step procedure, and an error correction model (ECM).

## Key results
- Growth is best described by an exponential trend; adding quarterly seasonality significantly improves the fit.
- The trend-and-seasonality residuals show ARCH effects, so an ARMA-GARCH model is used; the out-of-sample ex post forecast tracks actual GDP closely.
- GDP and household consumption share a long-run equilibrium: the ECM error-correction term is negative and significant, i.e. deviations are corrected over time.

## How to run
1. Clone the repository and open `time-series-econometrics-r.Rproj` in RStudio.
2. Install the required packages once:
   install.packages(c("here", "readxl", "forecast", "rugarch", "FinTS", "urca", "rmarkdown"))
3. Open `analysis.Rmd` and click Knit to reproduce `analysis.html` and the figures, or run `philippines_gdp_timeseries.R` line by line.

