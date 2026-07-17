############################################################################
## Adaptador QBART
##
## La corrida de investigacion usa por defecto el kernel C++ original de
## BayesQArt. El motor en R de notebooks/bayesqart.R se conserva solamente
## para desarrollo y debe solicitarse con engine = "local_experimental".
############################################################################

hp_qbart_default <- list(
  engine = "official",
  allow_experimental = FALSE,
  m = 200L,
  burn = 3000L,
  ndpost = 5000L,
  nc = 50L,
  min_obs = 5L,
  kappa = 2.0,
  maxdepth = 3L,
  aa = 2.0,
  bb = 3.0,
  pbd = 0.4,
  pb = 0.5,
  alpha = 0.95,
  beta = 2.0,
  cutpoint_strategy = "quantile_rank",
  seed = 20222834L
)

qbart_tau_name <- function(tau) {
  paste0("tau_", format(tau, scientific = FALSE, trim = TRUE))
}

qbart_hp_for_tau <- function(hp, tau) {
  stopifnot(length(tau) == 1L, is.finite(tau), tau > 0, tau < 1)
  by_tau <- hp$by_tau
  hp$by_tau <- NULL
  if (is.null(by_tau)) return(hp)

  keys <- c(qbart_tau_name(tau), format(tau, scientific = FALSE, trim = TRUE))
  hit <- keys[keys %in% names(by_tau)]
  if (length(hit) == 0L) return(hp)
  utils::modifyList(hp, by_tau[[hit[1L]]])
}

.qbart_validate_inputs <- function(x_train, y_train, x_pred, tau) {
  x_train <- as.matrix(x_train)
  x_pred <- as.matrix(x_pred)
  storage.mode(x_train) <- "double"
  storage.mode(x_pred) <- "double"
  y_train <- as.numeric(y_train)

  if (nrow(x_train) != length(y_train)) {
    stop("x_train y y_train tienen tamanos incompatibles.")
  }
  if (ncol(x_train) != ncol(x_pred)) {
    stop("x_train y x_pred deben tener las mismas columnas.")
  }
  if (anyNA(x_train) || anyNA(x_pred) || anyNA(y_train)) {
    stop("QBART no admite valores faltantes en la matriz de modelacion.")
  }
  if (!is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("tau debe pertenecer al intervalo (0, 1).")
  }
  if (is.null(colnames(x_train))) {
    colnames(x_train) <- paste0("X", seq_len(ncol(x_train)))
  }
  colnames(x_pred) <- colnames(x_train)
  list(x_train = x_train, y_train = y_train, x_pred = x_pred)
}

.qbart_rank_transform <- function(x_train, x_pred) {
  x_train_out <- x_train
  x_pred_out <- x_pred
  n <- nrow(x_train)

  for (j in seq_len(ncol(x_train))) {
    xj <- x_train[, j]
    if (length(unique(xj)) <= 1L) {
      x_train_out[, j] <- 0.5
      x_pred_out[, j] <- 0.5
      next
    }
    empirical_cdf <- stats::ecdf(xj)
    x_train_out[, j] <- pmin(pmax(empirical_cdf(xj), 1 / (n + 1)), n / (n + 1))
    x_pred_out[, j] <- pmin(
      pmax(empirical_cdf(x_pred[, j]), 1 / (n + 1)),
      n / (n + 1)
    )
  }
  list(x_train = x_train_out, x_pred = x_pred_out)
}

.qbart_transform_predictors <- function(x_train, x_pred, strategy) {
  strategy <- match.arg(strategy, c("quantile_rank", "none"))
  if (strategy == "none") {
    return(list(x_train = x_train, x_pred = x_pred))
  }
  .qbart_rank_transform(x_train, x_pred)
}

.qbart_source_local_engine <- function() {
  if (exists("bayesqart", mode = "function", inherits = TRUE)) return(invisible(TRUE))

  root <- Sys.getenv("QBART_PROJECT_ROOT", unset = getwd())
  candidates <- c(
    file.path(root, "notebooks", "bayesqart.R"),
    file.path(root, "bayesqart.R"),
    file.path(getwd(), "notebooks", "bayesqart.R"),
    file.path(getwd(), "bayesqart.R")
  )
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) {
    stop("No se encontro notebooks/bayesqart.R para el motor experimental.")
  }
  source(path, local = .GlobalEnv)
  invisible(TRUE)
}

qbart_engine_available <- function(engine = c("official", "local_experimental")) {
  engine <- match.arg(engine)
  if (engine == "official") {
    return(
      requireNamespace("BayesQArt", quietly = TRUE) &&
        requireNamespace("Rcpp", quietly = TRUE)
    )
  }
  TRUE
}

.qbart_resolve_engine <- function(hp) {
  engine <- hp$engine %||% "official"
  engine <- match.arg(engine, c("official", "auto", "local_experimental"))

  if (engine == "auto") {
    if (qbart_engine_available("official")) return("official")
    if (isTRUE(hp$allow_experimental)) return("local_experimental")
    stop(
      "BayesQArt oficial no esta instalado. La corrida no usara el sampler ",
      "experimental salvo que allow_experimental = TRUE."
    )
  }
  if (engine == "official" && !qbart_engine_available("official")) {
    stop(
      "Falta el paquete oficial BayesQArt. Instale el paquete desde ",
      "https://github.com/bpkindo/bayesqart antes de la corrida."
    )
  }
  if (engine == "local_experimental" && !isTRUE(hp$allow_experimental)) {
    stop(
      "El motor local es experimental. Para una prueba explicita use ",
      "allow_experimental = TRUE; no lo use para resultados finales."
    )
  }
  engine
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.fit_qbart_una <- function(x_train, y_train, x_pred, tau,
                           hp = hp_qbart_default) {
  hp <- qbart_hp_for_tau(hp, tau)
  checked <- .qbart_validate_inputs(x_train, y_train, x_pred, tau)
  transformed <- .qbart_transform_predictors(
    checked$x_train,
    checked$x_pred,
    hp$cutpoint_strategy %||% "quantile_rank"
  )
  x_train <- transformed$x_train
  x_pred <- transformed$x_pred
  y_train <- checked$y_train
  engine <- .qbart_resolve_engine(hp)

  seed <- as.integer((hp$seed %||% 20222834L) + round(10000 * tau))
  set.seed(seed)

  if (engine == "official") {
    suppressPackageStartupMessages(library(Rcpp))
    loadNamespace("BayesQArt")
    fit_fun <- getFromNamespace("BayesQArt", "BayesQArt")
    out <- fit_fun(
      y = y_train,
      X = x_train,
      Xtest = x_pred,
      quantile = tau,
      burn = as.integer(hp$burn),
      nd = as.integer(hp$ndpost),
      m = as.integer(hp$m),
      min_obs_node = as.integer(hp$min_obs),
      aa_parm = hp$aa,
      bb_parm = hp$bb,
      nc = as.integer(hp$nc),
      pbd = hp$pbd %||% 0.4,
      pb = hp$pb %||% 0.5,
      alpha = hp$alpha %||% 0.95,
      betap = hp$beta %||% 2.0,
      kappa = hp$kappa,
      maxdepth = as.integer(hp$maxdepth)
    )
    used <- as.integer(out$vars_used)
    used <- used[used >= 0L & used < ncol(x_train)]
    counts <- tabulate(used + 1L, nbins = ncol(x_train))
    imp <- if (sum(counts) > 0) counts / sum(counts) else counts
    imp <- setNames(imp, colnames(x_train))
  } else {
    .qbart_source_local_engine()
    out <- bayesqart(
      y = y_train,
      X = x_train,
      Xtest = x_pred,
      tau = tau,
      n_trees = as.integer(hp$m),
      burn = as.integer(hp$burn),
      n_draws = as.integer(hp$ndpost),
      nc = as.integer(hp$nc),
      min_obs = as.integer(hp$min_obs),
      alpha = hp$alpha %||% 0.95,
      beta = hp$beta %||% 2.0,
      kappa = hp$kappa,
      maxdepth = as.integer(hp$maxdepth),
      aa = hp$aa,
      bb = hp$bb,
      verbose = 0
    )
    imp <- out$var_importance
  }

  list(
    pred = as.numeric(out$pred_test),
    imp = imp[colnames(x_train)],
    engine = engine,
    tau = tau,
    seed = seed
  )
}

reordenar_cuantiles <- function(q, taus = NULL) {
  q <- as.matrix(q)
  if (ncol(q) <= 1L) return(q)
  out <- t(apply(q, 1L, sort))
  if (!is.null(taus)) colnames(out) <- vapply(taus, qbart_tau_name, character(1))
  out
}

predecir_qbart_grid <- function(x_train, y_train, x_pred, taus_q,
                                hp = hp_qbart_default,
                                rearrange = TRUE,
                                return_details = FALSE) {
  x_pred <- as.matrix(x_pred)
  fits <- lapply(
    taus_q,
    function(tt) .fit_qbart_una(x_train, y_train, x_pred, tt, hp)
  )
  q <- vapply(fits, function(fit) fit$pred, numeric(nrow(x_pred)))
  if (is.null(dim(q))) q <- matrix(q, nrow = nrow(x_pred))
  colnames(q) <- vapply(taus_q, qbart_tau_name, character(1))
  if (isTRUE(rearrange)) q <- reordenar_cuantiles(q, taus_q)

  if (isTRUE(return_details)) return(list(pred = q, fits = fits))
  q
}

expandir_a_grilla_fina <- function(q_coarse, taus_q, taus_fine) {
  t(apply(q_coarse, 1L, function(qc) {
    stats::approx(x = taus_q, y = qc, xout = taus_fine, rule = 2)$y
  }))
}

importancia_qbart_por_cuantil <- function(x, y, taus_q,
                                          hp = hp_qbart_default) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  x_dummy <- x[seq_len(min(5L, nrow(x))), , drop = FALSE]

  # BayesQArt registra la variable de cada nacimiento de rama aceptado,
  # incluido el burn-in. Es una medida exploratoria, no una PIP posterior.
  cols <- lapply(taus_q, function(tt) {
    fit <- .fit_qbart_una(x, y, x_dummy, tt, hp)
    fit$imp[colnames(x)]
  })
  importance <- do.call(cbind, cols)
  colnames(importance) <- vapply(taus_q, qbart_tau_name, character(1))
  rownames(importance) <- colnames(x)
  table <- data.frame(
    variable = rownames(importance),
    importance,
    row.names = NULL,
    check.names = FALSE
  )

  center <- qbart_tau_name(taus_q[which.min(abs(taus_q - 0.5))])
  for (column in grep("^tau_", names(table), value = TRUE)) {
    if (column == center) next
    denominator <- table[[center]]
    table[[paste0("ratio_", column, "_centro")]] <- ifelse(
      denominator > 0,
      table[[column]] / denominator,
      NA_real_
    )
  }
  table
}
