############################################################################
## Construccion causal de la base para pronostico de cuantiles cambiarios.
############################################################################

future_return <- function(daily_return, horizon = 1L) {
  horizon <- as.integer(horizon)
  n <- length(daily_return)
  out <- rep(NA_real_, n)
  if (horizon < 1L || n <= horizon) return(out)
  for (time in seq_len(n - horizon)) {
    out[time] <- sum(daily_return[(time + 1L):(time + horizon)])
  }
  out
}

log_change <- function(x) {
  c(NA_real_, 100 * diff(log(x)))
}

asinh_change <- function(x) {
  scale <- stats::median(abs(x[is.finite(x)]), na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- 1
  c(NA_real_, 100 * diff(asinh(x / scale)))
}

rolling_stat <- function(x, window, fun, min_obs = window) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (time in seq_len(n)) {
    first <- max(1L, time - window + 1L)
    values <- x[first:time]
    values <- values[is.finite(values)]
    if (length(values) >= min_obs) out[time] <- fun(values)
  }
  out
}

rolling_sd <- function(x, window) {
  rolling_stat(x, window, stats::sd, min_obs = window)
}

rolling_mean_abs <- function(x, window) {
  rolling_stat(x, window, function(values) mean(abs(values)), min_obs = window)
}

rolling_skewness <- function(x, window) {
  rolling_stat(x, window, function(values) {
    scale <- stats::sd(values)
    if (!is.finite(scale) || scale == 0) return(0)
    mean(((values - mean(values)) / scale)^3)
  }, min_obs = window)
}

rolling_kurtosis <- function(x, window) {
  rolling_stat(x, window, function(values) {
    scale <- stats::sd(values)
    if (!is.finite(scale) || scale == 0) return(0)
    mean(((values - mean(values)) / scale)^4) - 3
  }, min_obs = window)
}

ewma_volatility <- function(x, lambda = 0.94) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  first <- which(is.finite(x))[1L]
  if (is.na(first)) return(out)
  initial <- stats::sd(x[seq_len(min(length(x), first + 19L))], na.rm = TRUE)
  if (!is.finite(initial)) initial <- abs(x[first])
  variance <- initial^2
  for (time in first:length(x)) {
    if (is.finite(x[time])) {
      variance <- lambda * variance + (1 - lambda) * x[time]^2
    }
    out[time] <- sqrt(variance)
  }
  out
}

build_qbart_dataset <- function(base, horizon = 1L) {
  required <- c(
    "fecha", "retorno_tc", "vix", "baa10y", "tasa_usa", "nfci",
    "yield_to_worst", "tp_1y", "tp_5y", "tp_10y", "embig_peru",
    "tasa_peru", "interv_bcrp", "rin", "oro", "cobre", "zinc",
    "plata", "wti", "spx", "mxef", "dxy", "jpy", "chf", "eur",
    "gbp", "clp", "cop", "mxn", "brl", "cny"
  )
  missing_columns <- setdiff(required, names(base))
  if (length(missing_columns) > 0L) {
    stop("Faltan columnas: ", paste(missing_columns, collapse = ", "))
  }

  base <- base[order(base$fecha), , drop = FALSE]
  base$retorno_tc <- 100 * as.numeric(base$retorno_tc)
  base$y_objetivo <- future_return(base$retorno_tc, horizon)

  level_variables <- c(
    "vix", "baa10y", "tasa_usa", "nfci", "yield_to_worst",
    "tp_1y", "tp_5y", "tp_10y", "embig_peru"
  )
  price_variables <- c(
    "rin", "oro", "cobre", "zinc", "plata", "spx", "mxef",
    "dxy", "jpy", "chf", "eur", "gbp", "clp", "cop", "mxn",
    "brl", "cny"
  )

  for (variable in price_variables) {
    base[[paste0("d_", variable)]] <- log_change(base[[variable]])
  }
  # WTI puede ser no positivo; asinh conserva el orden sin producir NaN.
  base$d_wti <- asinh_change(base$wti)
  base$d_vix <- log_change(base$vix)
  base$diff_tasas <- base$tasa_peru - base$tasa_usa
  base$interv_bcrp_missing <- as.integer(is.na(base$interv_bcrp))
  base$interv_bcrp[is.na(base$interv_bcrp)] <- 0

  for (lag in 1:5) {
    shift <- lag - 1L
    values <- if (shift == 0L) {
      base$retorno_tc
    } else {
      c(rep(NA_real_, shift), head(base$retorno_tc, -shift))
    }
    base[[paste0("ret_lag_", lag)]] <- values
  }
  base$abs_ret_lag_1 <- abs(base$ret_lag_1)
  base$sq_ret_lag_1 <- base$ret_lag_1^2
  base$vol_sd_5 <- rolling_sd(base$retorno_tc, 5L)
  base$vol_sd_20 <- rolling_sd(base$retorno_tc, 20L)
  base$vol_sd_60 <- rolling_sd(base$retorno_tc, 60L)
  base$vol_ewma_094 <- ewma_volatility(base$retorno_tc, 0.94)
  base$mean_abs_ret_20 <- rolling_mean_abs(base$retorno_tc, 20L)
  base$skew_ret_20 <- rolling_skewness(base$retorno_tc, 20L)
  base$kurt_ret_60 <- rolling_kurtosis(base$retorno_tc, 60L)
  base$vol_ratio_5_60 <- base$vol_sd_5 / pmax(base$vol_sd_60, 1e-8)

  predictors <- c(
    level_variables,
    "d_vix",
    "diff_tasas",
    "interv_bcrp",
    "interv_bcrp_missing",
    paste0("d_", price_variables),
    "d_wti",
    paste0("ret_lag_", 1:5),
    "abs_ret_lag_1",
    "sq_ret_lag_1",
    "vol_sd_5",
    "vol_sd_20",
    "vol_sd_60",
    "vol_ewma_094",
    "mean_abs_ret_20",
    "skew_ret_20",
    "kurt_ret_60",
    "vol_ratio_5_60"
  )
  model_columns <- c("fecha", "y_objetivo", predictors)
  model_data <- base[stats::complete.cases(base[, model_columns]), model_columns]

  list(data = model_data, predictors = predictors)
}
