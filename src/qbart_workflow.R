############################################################################
## Utilidades reproducibles para seleccion, calibracion y evaluacion QBART.
############################################################################

pinball_loss <- function(y, q, tau) {
  error <- as.numeric(y) - as.numeric(q)
  mean(error * (tau - as.numeric(error < 0)))
}

quantile_scores <- function(predictions, y, taus) {
  predictions <- as.matrix(predictions)
  stopifnot(ncol(predictions) == length(taus))
  vapply(
    seq_along(taus),
    function(j) pinball_loss(y, predictions[, j], taus[j]),
    numeric(1)
  )
}

quantile_coverage <- function(predictions, y, taus) {
  predictions <- as.matrix(predictions)
  stopifnot(ncol(predictions) == length(taus))
  vapply(
    seq_along(taus),
    function(j) mean(as.numeric(y) <= predictions[, j]),
    numeric(1)
  )
}

.qbart_grid_row <- function(grid, row) {
  out <- lapply(grid[row, , drop = FALSE], function(value) {
    value <- value[[1L]]
    if (is.factor(value)) as.character(value) else value
  })
  names(out) <- names(grid)
  out
}

tune_qbart_by_tau <- function(x_train, y_train, x_validation, y_validation,
                              taus, grid, hp_base, seeds = 20222834L,
                              verbose = TRUE) {
  if (nrow(grid) < 1L) stop("La grilla QBART esta vacia.")
  seeds <- as.integer(seeds)
  results <- vector("list", length(taus) * nrow(grid))
  result_index <- 0L

  for (tau in taus) {
    for (grid_index in seq_len(nrow(grid))) {
      config <- .qbart_grid_row(grid, grid_index)
      losses <- numeric(length(seeds))
      start_time <- Sys.time()

      for (seed_index in seq_along(seeds)) {
        hp <- utils::modifyList(hp_base, config)
        hp$seed <- seeds[seed_index]
        prediction <- .fit_qbart_una(
          x_train,
          y_train,
          x_validation,
          tau,
          hp
        )$pred
        losses[seed_index] <- pinball_loss(y_validation, prediction, tau)
      }

      result_index <- result_index + 1L
      result <- data.frame(
        tau = tau,
        as.data.frame(config, stringsAsFactors = FALSE),
        pinball = mean(losses),
        pinball_sd = if (length(losses) > 1L) stats::sd(losses) else 0,
        seeds = paste(seeds, collapse = ","),
        minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
        check.names = FALSE
      )
      results[[result_index]] <- result

      if (isTRUE(verbose)) {
        config_text <- paste(
          paste(names(config), unlist(config), sep = "="),
          collapse = ", "
        )
        message(
          sprintf(
            "[tau=%.2f, %d/%d] %s -> pinball=%.6f",
            tau,
            grid_index,
            nrow(grid),
            config_text,
            result$pinball
          )
        )
      }
    }
  }

  table <- do.call(rbind, results)
  table <- table[order(table$tau, table$pinball), , drop = FALSE]
  best_rows <- unlist(
    lapply(taus, function(tau) {
      candidates <- which(table$tau == tau)
      candidates[which.min(table$pinball[candidates])]
    }),
    use.names = FALSE
  )
  best <- table[best_rows, , drop = FALSE]
  rownames(table) <- NULL
  rownames(best) <- NULL

  list(
    table = table,
    best = best,
    config_columns = names(grid)
  )
}

apply_qbart_tuning <- function(hp_full, tuning) {
  by_tau <- list()
  for (row in seq_len(nrow(tuning$best))) {
    tau <- tuning$best$tau[row]
    config <- lapply(tuning$config_columns, function(column) {
      value <- tuning$best[[column]][row]
      if (is.factor(value)) as.character(value) else value
    })
    names(config) <- tuning$config_columns
    by_tau[[qbart_tau_name(tau)]] <- config
  }
  hp_full$by_tau <- by_tau
  hp_full
}

calibration_offsets <- function(y_history, q_history, taus,
                                window = 500L, min_obs = 100L) {
  q_history <- as.matrix(q_history)
  if (ncol(q_history) != length(taus)) {
    stop("q_history debe tener una columna por tau.")
  }
  valid_rows <- is.finite(y_history) & apply(is.finite(q_history), 1L, all)
  y_history <- as.numeric(y_history[valid_rows])
  q_history <- q_history[valid_rows, , drop = FALSE]

  if (length(y_history) < min_obs) return(rep(0, length(taus)))
  if (is.finite(window) && window > 0L && length(y_history) > window) {
    keep <- seq.int(length(y_history) - window + 1L, length(y_history))
    y_history <- y_history[keep]
    q_history <- q_history[keep, , drop = FALSE]
  }

  vapply(seq_along(taus), function(j) {
    stats::quantile(
      y_history - q_history[, j],
      probs = taus[j],
      names = FALSE,
      type = 8
    )
  }, numeric(1))
}

apply_quantile_calibration <- function(raw_predictions, offsets, taus,
                                       rearrange = TRUE) {
  calibrated <- sweep(as.matrix(raw_predictions), 2L, offsets, "+")
  if (isTRUE(rearrange)) reordenar_cuantiles(calibrated, taus) else calibrated
}

tune_qrf <- function(x_train, y_train, x_validation, y_validation,
                     taus, grid, seeds = 20222834L,
                     weights = NULL, verbose = TRUE) {
  if (!requireNamespace("quantregForest", quietly = TRUE)) {
    stop("Falta el paquete quantregForest.")
  }
  if (is.null(weights)) {
    weights <- ifelse(abs(taus - 0.5) > 0.25, 0.4, 0.2)
  }
  weights <- weights / sum(weights)
  baseline <- vapply(
    taus,
    function(tau) {
      q0 <- stats::quantile(y_train, tau, names = FALSE)
      pinball_loss(y_validation, rep(q0, length(y_validation)), tau)
    },
    numeric(1)
  )
  baseline <- pmax(baseline, .Machine$double.eps)
  results <- vector("list", nrow(grid))

  for (grid_index in seq_len(nrow(grid))) {
    config <- .qbart_grid_row(grid, grid_index)
    scores <- numeric(length(seeds))
    raw_losses <- matrix(NA_real_, nrow = length(seeds), ncol = length(taus))
    start_time <- Sys.time()

    for (seed_index in seq_along(seeds)) {
      set.seed(seeds[seed_index])
      fit <- quantregForest::quantregForest(
        x = as.matrix(x_train),
        y = as.numeric(y_train),
        ntree = as.integer(config$ntree),
        mtry = as.integer(config$mtry),
        nodesize = as.integer(config$nodesize),
        keep.inbag = FALSE
      )
      prediction <- predict(
        fit,
        newdata = as.matrix(x_validation),
        what = taus
      )
      raw_losses[seed_index, ] <- quantile_scores(
        prediction,
        y_validation,
        taus
      )
      scores[seed_index] <- sum(weights * raw_losses[seed_index, ] / baseline)
    }

    results[[grid_index]] <- data.frame(
      as.data.frame(config, stringsAsFactors = FALSE),
      relative_tail_score = mean(scores),
      score_sd = if (length(scores) > 1L) stats::sd(scores) else 0,
      setNames(
        as.list(colMeans(raw_losses)),
        paste0("pinball_", vapply(taus, qbart_tau_name, character(1)))
      ),
      minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
      check.names = FALSE
    )
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "[QRF %d/%d] ntree=%d, mtry=%d, nodesize=%d -> score=%.5f",
          grid_index,
          nrow(grid),
          config$ntree,
          config$mtry,
          config$nodesize,
          results[[grid_index]]$relative_tail_score
        )
      )
    }
  }

  table <- do.call(rbind, results)
  table <- table[order(table$relative_tail_score), , drop = FALSE]
  rownames(table) <- NULL
  list(table = table, best = table[1L, , drop = FALSE])
}

kupiec_test <- function(y, q, tau) {
  hit <- as.numeric(as.numeric(y) <= as.numeric(q))
  hit <- hit[is.finite(hit)]
  n <- length(hit)
  exceptions <- sum(hit)
  empirical <- exceptions / n
  log_term <- function(count, probability) {
    if (count == 0L) return(0)
    count * log(pmax(probability, .Machine$double.eps))
  }
  null_loglik <- log_term(exceptions, tau) +
    log_term(n - exceptions, 1 - tau)
  alternative_loglik <- log_term(exceptions, empirical) +
    log_term(n - exceptions, 1 - empirical)
  statistic <- -2 * (null_loglik - alternative_loglik)
  c(
    statistic = statistic,
    p_value = stats::pchisq(statistic, df = 1L, lower.tail = FALSE),
    coverage = empirical
  )
}

dq_test <- function(y, q, tau, lags = 4L) {
  hit <- as.numeric(as.numeric(y) <= as.numeric(q)) - tau
  valid <- is.finite(hit)
  hit <- hit[valid]
  n <- length(hit)
  if (n <= lags + 5L) {
    return(c(statistic = NA_real_, p_value = NA_real_))
  }

  response <- hit[(lags + 1L):n]
  regressors <- cbind(
    intercept = 1,
    vapply(
      seq_len(lags),
      function(lag) hit[(lags + 1L - lag):(n - lag)],
      numeric(n - lags)
    )
  )
  moment <- crossprod(regressors, response)
  moment_variance <- tau * (1 - tau) * crossprod(regressors)
  statistic <- tryCatch(
    as.numeric(crossprod(moment, qr.solve(moment_variance, moment))),
    error = function(error) NA_real_
  )
  c(
    statistic = statistic,
    p_value = stats::pchisq(
      statistic,
      df = ncol(regressors),
      lower.tail = FALSE
    )
  )
}

dm_pinball_test <- function(y, q_model, q_benchmark, tau, hac_lag = NULL) {
  y <- as.numeric(y)
  model_error <- y - as.numeric(q_model)
  benchmark_error <- y - as.numeric(q_benchmark)
  model_loss <- model_error * (tau - as.numeric(model_error < 0))
  benchmark_loss <- benchmark_error * (tau - as.numeric(benchmark_error < 0))
  difference <- model_loss - benchmark_loss
  difference <- difference[is.finite(difference)]
  n <- length(difference)
  if (n < 10L) return(c(statistic = NA_real_, p_value = NA_real_))
  if (is.null(hac_lag)) hac_lag <- max(1L, floor(n^(1 / 3)))
  hac_lag <- min(as.integer(hac_lag), n - 1L)

  centered <- difference - mean(difference)
  long_run_variance <- sum(centered^2) / n
  for (lag in seq_len(hac_lag)) {
    covariance <- sum(
      centered[(lag + 1L):n] * centered[seq_len(n - lag)]
    ) / n
    weight <- 1 - lag / (hac_lag + 1)
    long_run_variance <- long_run_variance + 2 * weight * covariance
  }
  standard_error <- sqrt(pmax(long_run_variance, 0) / n)
  statistic <- if (standard_error > 0) {
    mean(difference) / standard_error
  } else {
    NA_real_
  }
  c(
    statistic = statistic,
    p_value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
    mean_loss_difference = mean(difference)
  )
}
