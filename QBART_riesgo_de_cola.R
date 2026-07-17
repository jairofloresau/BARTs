############################################################################
## QBART para riesgo de cola del tipo de cambio
##
## Uso:
##   Rscript QBART_riesgo_de_cola.R
##
## Variables de entorno opcionales:
##   QBART_MODE=prueba|exploracion|completa
##   QBART_WINDOW_OBS=2520       # 0 usa ventana expansiva
##   QBART_REESTIMATION_STEP=60
##   QBART_MAX_TEST_ORIGINS=0    # 0 usa todo el test
############################################################################

options(stringsAsFactors = FALSE)

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "quantreg", "quantregForest",
  "digest", "Rcpp", "BayesQArt"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Faltan paquetes requeridos: ",
    paste(missing_packages, collapse = ", "),
    ". BayesQArt se instala desde https://github.com/bpkindo/bayesqart."
  )
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(quantreg)
  library(quantregForest)
})

script_path <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_root <- if (length(script_path) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_path[1L]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
if (!file.exists(file.path(project_root, "notebooks", "qbart_adapter.R"))) {
  stop("Ejecute el script desde la raiz del repositorio BARTs.")
}
Sys.setenv(QBART_PROJECT_ROOT = project_root)
source(file.path(project_root, "notebooks", "qbart_adapter.R"))
source(file.path(project_root, "src", "qbart_features.R"))
source(file.path(project_root, "src", "qbart_workflow.R"))

seed_base <- 20222834L
set.seed(seed_base)
taus <- c(0.05, 0.50, 0.95)
horizon <- 1L
mode <- match.arg(
  Sys.getenv("QBART_MODE", unset = "prueba"),
  c("prueba", "exploracion", "completa")
)

data_dir <- file.path(project_root, "data", "datos", "BASE DIARIO")
output_dir <- Sys.getenv(
  "QBART_OUTPUT_DIR",
  unset = file.path(project_root, "data", "processed")
)
figure_dir <- Sys.getenv(
  "QBART_FIGURE_DIR",
  unset = file.path(project_root, "reports", "figures")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value)) default else value
}

if (mode == "prueba") {
  hp_full <- utils::modifyList(
    hp_qbart_default,
    list(m = 20L, burn = 50L, ndpost = 100L)
  )
  hp_tune <- utils::modifyList(hp_full, list(burn = 20L, ndpost = 30L))
  qbart_grid <- expand.grid(
    kappa = c(1, 2),
    m = c(10L, 20L),
    min_obs = 5L,
    maxdepth = 3L
  )
  tuning_seeds <- seed_base
  default_step <- 200L
  default_max_origins <- 200L
} else if (mode == "exploracion") {
  hp_full <- utils::modifyList(
    hp_qbart_default,
    list(m = 100L, burn = 1000L, ndpost = 2000L)
  )
  hp_tune <- utils::modifyList(hp_full, list(burn = 200L, ndpost = 300L))
  qbart_grid <- expand.grid(
    kappa = c(0.5, 1, 2),
    m = c(50L, 100L),
    min_obs = c(5L, 15L),
    maxdepth = 3L
  )
  tuning_seeds <- seed_base + c(0L, 1000L)
  default_step <- 60L
  default_max_origins <- 0L
} else {
  hp_full <- hp_qbart_default
  hp_tune <- utils::modifyList(hp_full, list(burn = 500L, ndpost = 1000L))
  qbart_grid <- expand.grid(
    kappa = c(0.5, 1, 2),
    m = c(100L, 200L),
    min_obs = c(5L, 15L),
    maxdepth = c(3L, 4L)
  )
  tuning_seeds <- seed_base + c(0L, 1000L, 2000L)
  default_step <- 60L
  default_max_origins <- 0L
}

window_obs <- read_env_integer("QBART_WINDOW_OBS", 2520L)
reestimation_step <- read_env_integer("QBART_REESTIMATION_STEP", default_step)
max_test_origins <- read_env_integer(
  "QBART_MAX_TEST_ORIGINS",
  default_max_origins
)
calibration_window <- read_env_integer("QBART_CALIBRATION_WINDOW", 500L)
calibration_min_obs <- read_env_integer("QBART_CALIBRATION_MIN_OBS", 100L)
suffix <- paste0("_", mode, "_v2")

message("=== QBART ", toupper(mode), " ===")
message("Motor: ", hp_full$engine, " | ventana: ", window_obs, " observaciones")

# -------------------------------------------------------------------------
# Datos y predictores causales
# -------------------------------------------------------------------------
data_path <- file.path(data_dir, "base_final_diaria.rds")
if (!file.exists(data_path)) stop("No se encontro: ", data_path)
base <- readRDS(data_path)
prepared <- build_qbart_dataset(base, horizon)
model_data <- prepared$data
predictors <- prepared$predictors

n_total <- nrow(model_data)
train_end <- floor(0.60 * n_total)
validation_end <- floor(0.80 * n_total)
index_train <- seq_len(train_end)
index_validation <- seq.int(train_end + 1L, validation_end)
index_test <- seq.int(validation_end + 1L, n_total)
if (max_test_origins > 0L) {
  index_test <- head(index_test, max_test_origins)
}

window_index <- function(last_index, maximum_length) {
  if (maximum_length <= 0L) return(seq_len(last_index))
  seq.int(max(1L, last_index - maximum_length + 1L), last_index)
}
index_tuning_train <- window_index(train_end, window_obs)

x_train <- as.matrix(model_data[index_tuning_train, predictors])
y_train <- model_data$y_objetivo[index_tuning_train]
x_validation <- as.matrix(model_data[index_validation, predictors])
y_validation <- model_data$y_objetivo[index_validation]

message(
  "Datos: train=", length(index_tuning_train),
  ", validacion=", length(index_validation),
  ", test=", length(index_test),
  ", predictores=", length(predictors)
)

# -------------------------------------------------------------------------
# Tuning QBART separado por tau
# -------------------------------------------------------------------------
message("=== Tuning QBART separado por percentil ===")
qbart_tuning <- tune_qbart_by_tau(
  x_train,
  y_train,
  x_validation,
  y_validation,
  taus,
  qbart_grid,
  hp_tune,
  seeds = tuning_seeds
)
hp_selected <- apply_qbart_tuning(hp_full, qbart_tuning)
utils::write.csv(
  qbart_tuning$table,
  file.path(output_dir, paste0("qbart_tuneo_por_tau", suffix, ".csv")),
  row.names = FALSE
)
print(qbart_tuning$best)

# -------------------------------------------------------------------------
# Tuning QRF justo: score relativo por tau con mayor peso en colas
# -------------------------------------------------------------------------
p <- length(predictors)
if (mode == "prueba") {
  qrf_grid <- data.frame(
    ntree = 100L,
    mtry = max(1L, floor(sqrt(p))),
    nodesize = 10L
  )
} else {
  qrf_grid <- expand.grid(
    ntree = if (mode == "completa") c(500L, 1000L) else c(300L, 500L),
    mtry = unique(pmax(1L, c(floor(sqrt(p)), floor(p / 3)))),
    nodesize = c(5L, 15L)
  )
}
message("=== Tuning QRF con score relativo y ponderado en colas ===")
qrf_tuning <- tune_qrf(
  x_train,
  y_train,
  x_validation,
  y_validation,
  taus,
  qrf_grid,
  seeds = tuning_seeds
)
utils::write.csv(
  qrf_tuning$table,
  file.path(output_dir, paste0("qrf_tuneo", suffix, ".csv")),
  row.names = FALSE
)
qrf_hp <- list(
  ntree = as.integer(qrf_tuning$best$ntree),
  mtry = as.integer(qrf_tuning$best$mtry),
  nodesize = as.integer(qrf_tuning$best$nodesize)
)

# La validacion constituye el historial inicial del calibrador. Se reestima
# con las cadenas finales y nunca usa observaciones del test futuro.
message("=== Prediccion de validacion para calibracion causal ===")
validation_raw <- predecir_qbart_grid(
  x_train,
  y_train,
  x_validation,
  taus,
  hp_selected,
  rearrange = TRUE
)

# -------------------------------------------------------------------------
# Evaluacion fuera de muestra
# -------------------------------------------------------------------------
n_origins <- length(index_test)
n_tau <- length(taus)
n_blocks <- ceiling(n_origins / reestimation_step)
y_test <- model_data$y_objetivo[index_test]
checkpoint_path <- file.path(
  output_dir,
  paste0("checkpoint_qbart", suffix, ".rds")
)
data_hash <- digest::digest(
  list(
    date = model_data$fecha,
    y = model_data$y_objetivo,
    predictors = predictors
  ),
  algo = "xxhash64"
)
run_config <- list(
  version = 2L,
  mode = mode,
  bayesqart_version = as.character(utils::packageVersion("BayesQArt")),
  data_hash = data_hash,
  taus = taus,
  hp_qbart = hp_selected,
  hp_qrf = qrf_hp,
  window_obs = window_obs,
  reestimation_step = reestimation_step,
  index_test = index_test,
  calibration_window = calibration_window,
  calibration_min_obs = calibration_min_obs,
  seed_base = seed_base
)

new_prediction_store <- function() {
  methods <- c(
    "random_walk", "quantile_reg", "quantile_rf",
    "qbart_raw", "qbart_calibrated"
  )
  setNames(
    lapply(methods, function(x) matrix(NA_real_, n_origins, n_tau)),
    methods
  )
}

predictions <- new_prediction_store()
calibration_log <- data.frame()
initial_block <- 1L
if (file.exists(checkpoint_path)) {
  checkpoint <- readRDS(checkpoint_path)
  if (isTRUE(all.equal(checkpoint$config, run_config))) {
    predictions <- checkpoint$predictions
    calibration_log <- checkpoint$calibration_log
    initial_block <- checkpoint$last_block + 1L
    if (!is.null(checkpoint$rng_state)) .Random.seed <- checkpoint$rng_state
    message("Checkpoint compatible: se reanuda en bloque ", initial_block)
  } else {
    message("Checkpoint incompatible: se inicia una corrida nueva.")
  }
}

formula_rq <- stats::as.formula(
  paste("y_objetivo ~", paste(predictors, collapse = " + "))
)
evaluation_start <- Sys.time()

if (initial_block <= n_blocks) {
  for (block in seq.int(initial_block, n_blocks)) {
    positions <- seq.int(
      (block - 1L) * reestimation_step + 1L,
      min(block * reestimation_step, n_origins)
    )
    time_indices <- index_test[positions]
    first_origin <- min(time_indices)
    training_indices <- window_index(first_origin - 1L, window_obs)
    training_data <- model_data[training_indices, , drop = FALSE]
    block_data <- model_data[time_indices, , drop = FALSE]
    x_history <- as.matrix(training_data[, predictors])
    y_history <- training_data$y_objetivo
    x_block <- as.matrix(block_data[, predictors])

    rw_quantiles <- stats::quantile(y_history, taus, names = FALSE)

    rq_fit <- quantreg::rq(formula_rq, tau = taus, data = training_data)
    rq_prediction <- reordenar_cuantiles(
      predict(rq_fit, newdata = block_data),
      taus
    )

    set.seed(seed_base + 50000L + block)
    qrf_fit <- quantregForest::quantregForest(
      x = x_history,
      y = y_history,
      ntree = qrf_hp$ntree,
      mtry = qrf_hp$mtry,
      nodesize = qrf_hp$nodesize,
      keep.inbag = FALSE
    )
    qrf_prediction <- predict(qrf_fit, newdata = x_block, what = taus)

    hp_block <- hp_selected
    hp_block$seed <- seed_base + 100000L * block
    qbart_raw <- predecir_qbart_grid(
      x_history,
      y_history,
      x_block,
      taus,
      hp_block,
      rearrange = TRUE
    )

    completed <- which(
      seq_len(n_origins) < min(positions) &
        apply(is.finite(predictions$qbart_raw), 1L, all)
    )
    calibration_y <- c(y_validation, y_test[completed])
    calibration_q <- rbind(
      validation_raw,
      predictions$qbart_raw[completed, , drop = FALSE]
    )
    offsets <- calibration_offsets(
      calibration_y,
      calibration_q,
      taus,
      window = calibration_window,
      min_obs = calibration_min_obs
    )
    qbart_calibrated <- apply_quantile_calibration(
      qbart_raw,
      offsets,
      taus
    )

    predictions$random_walk[positions, ] <- matrix(
      rw_quantiles,
      nrow = length(positions),
      ncol = n_tau,
      byrow = TRUE
    )
    predictions$quantile_reg[positions, ] <- rq_prediction
    predictions$quantile_rf[positions, ] <- qrf_prediction
    predictions$qbart_raw[positions, ] <- qbart_raw
    predictions$qbart_calibrated[positions, ] <- qbart_calibrated

    calibration_log <- rbind(
      calibration_log,
      data.frame(
        block = block,
        date = min(block_data$fecha),
        tau = taus,
        offset = offsets,
        calibration_n = length(calibration_y)
      )
    )
    saveRDS(
      list(
        predictions = predictions,
        calibration_log = calibration_log,
        last_block = block,
        config = run_config,
        rng_state = .Random.seed
      ),
      checkpoint_path
    )

    elapsed <- as.numeric(difftime(Sys.time(), evaluation_start, units = "mins"))
    remaining <- elapsed / (block - initial_block + 1L) * (n_blocks - block)
    message(
      sprintf(
        "[bloque %d/%d] %s | %.1f min, faltan ~%.1f min",
        block,
        n_blocks,
        as.character(max(training_data$fecha)),
        elapsed,
        remaining
      )
    )
  }
}

utils::write.csv(
  calibration_log,
  file.path(output_dir, paste0("qbart_calibracion", suffix, ".csv")),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Resultados
# -------------------------------------------------------------------------
score_table <- data.frame(tau = taus)
coverage_table <- data.frame(tau = taus)
for (method in names(predictions)) {
  score_table[[method]] <- quantile_scores(predictions[[method]], y_test, taus)
  coverage_table[[method]] <- quantile_coverage(
    predictions[[method]],
    y_test,
    taus
  )
}
relative_table <- score_table
relative_table[, -1L] <- sweep(
  relative_table[, -1L, drop = FALSE],
  1L,
  relative_table$random_walk,
  "/"
)

utils::write.csv(
  score_table,
  file.path(output_dir, paste0("qbart_quantile_score", suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  coverage_table,
  file.path(output_dir, paste0("qbart_cobertura", suffix, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  relative_table,
  file.path(output_dir, paste0("qbart_score_relativo", suffix, ".csv")),
  row.names = FALSE
)

diagnostic_rows <- list()
diagnostic_index <- 0L
for (method in names(predictions)) {
  for (tau_index in seq_along(taus)) {
    tau <- taus[tau_index]
    kupiec <- kupiec_test(y_test, predictions[[method]][, tau_index], tau)
    dq <- dq_test(y_test, predictions[[method]][, tau_index], tau)
    diagnostic_index <- diagnostic_index + 1L
    diagnostic_rows[[diagnostic_index]] <- data.frame(
      method = method,
      tau = tau,
      kupiec_statistic = kupiec["statistic"],
      kupiec_p_value = kupiec["p_value"],
      dq_statistic = dq["statistic"],
      dq_p_value = dq["p_value"]
    )
  }
}
diagnostic_table <- do.call(rbind, diagnostic_rows)
utils::write.csv(
  diagnostic_table,
  file.path(output_dir, paste0("qbart_backtests", suffix, ".csv")),
  row.names = FALSE
)

dm_rows <- list()
dm_index <- 0L
for (method in c("qbart_raw", "qbart_calibrated")) {
  for (tau_index in seq_along(taus)) {
    dm <- dm_pinball_test(
      y_test,
      predictions[[method]][, tau_index],
      predictions$quantile_rf[, tau_index],
      taus[tau_index]
    )
    dm_index <- dm_index + 1L
    dm_rows[[dm_index]] <- data.frame(
      method = method,
      benchmark = "quantile_rf",
      tau = taus[tau_index],
      statistic = dm["statistic"],
      p_value = dm["p_value"],
      mean_loss_difference = dm["mean_loss_difference"]
    )
  }
}
dm_table <- do.call(rbind, dm_rows)
utils::write.csv(
  dm_table,
  file.path(output_dir, paste0("qbart_dm_vs_qrf", suffix, ".csv")),
  row.names = FALSE
)

message("=== Quantile score ===")
print(round(score_table, 6), row.names = FALSE)
message("=== Cobertura ===")
print(round(coverage_table, 4), row.names = FALSE)

# -------------------------------------------------------------------------
# Determinantes por percentil
# -------------------------------------------------------------------------
determinant_indices <- window_index(validation_end, window_obs)
x_determinants <- as.matrix(model_data[determinant_indices, predictors])
y_determinants <- model_data$y_objetivo[determinant_indices]
importance_table <- importancia_qbart_por_cuantil(
  x_determinants,
  y_determinants,
  taus,
  hp_selected
)
utils::write.csv(
  importance_table,
  file.path(output_dir, paste0("qbart_importancia_por_tau", suffix, ".csv")),
  row.names = FALSE
)

tau_columns <- vapply(taus, qbart_tau_name, character(1))
heat_data <- importance_table |>
  dplyr::select(variable, dplyr::all_of(tau_columns)) |>
  tidyr::pivot_longer(
    -variable,
    names_to = "quantile",
    values_to = "importance"
  ) |>
  dplyr::group_by(variable) |>
  dplyr::mutate(
    relative_importance = ifelse(
      max(importance) > 0,
      importance / max(importance),
      0
    )
  ) |>
  dplyr::ungroup()

importance_plot <- ggplot(
  heat_data,
  aes(
    x = quantile,
    y = reorder(variable, relative_importance),
    fill = relative_importance
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#f7f3e8", high = "#a3261f") +
  labs(
    title = "QBART: uso de predictores por percentil",
    subtitle = paste(
      "Nacimientos de ramas aceptados (incluye burn-in),",
      "normalizados dentro de variable"
    ),
    x = "Percentil",
    y = NULL,
    fill = "Importancia\nrelativa"
  ) +
  theme_minimal(base_family = "serif")
ggsave(
  file.path(
    figure_dir,
    paste0("qbart_importancia_por_tau", suffix, ".png")
  ),
  importance_plot,
  width = 8,
  height = 10
)

message("=== Resumen ===")
message("Motor QBART: ", hp_selected$engine)
message("Ventana: ", window_obs, " | paso: ", reestimation_step)
for (tau in taus) {
  hp_tau <- qbart_hp_for_tau(hp_selected, tau)
  message(
    sprintf(
      "tau=%.2f: kappa=%g, m=%d, min_obs=%d, maxdepth=%d",
      tau,
      hp_tau$kappa,
      hp_tau$m,
      hp_tau$min_obs,
      hp_tau$maxdepth
    )
  )
}
message("Resultados guardados en: ", output_dir)
if (mode != "completa") {
  message(
    "La corrida no es final. Use QBART_MODE=completa despues de superar ",
    "las pruebas sinteticas."
  )
}
