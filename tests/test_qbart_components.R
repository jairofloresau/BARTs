args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1L]))
} else {
  normalizePath("tests/test_qbart_components.R")
}
project_root <- dirname(dirname(script_file))
Sys.setenv(QBART_PROJECT_ROOT = project_root)
source(file.path(project_root, "notebooks", "qbart_adapter.R"))
source(file.path(project_root, "src", "qbart_features.R"))
source(file.path(project_root, "src", "qbart_workflow.R"))

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert_true(hp_qbart_default$engine == "official", "El motor por defecto no es oficial.")
assert_true(hp_qbart_default$aa == 2, "El prior aa debe seguir la referencia oficial.")
assert_true(hp_qbart_default$bb == 3, "El prior bb debe seguir la referencia oficial.")

hp_tau <- hp_qbart_default
hp_tau$by_tau <- list(
  tau_0.05 = list(kappa = 0.5),
  tau_0.95 = list(kappa = 2)
)
assert_true(qbart_hp_for_tau(hp_tau, 0.05)$kappa == 0.5, "Fallo hp por tau.")
assert_true(qbart_hp_for_tau(hp_tau, 0.95)$kappa == 2, "Fallo hp por tau.")

x_train <- matrix(c(4, 1, 3, 2, 10, 20, 30, 40), ncol = 2)
x_test <- matrix(c(0, 2.5, 5, 15, 25, 50), ncol = 2)
ranked <- .qbart_rank_transform(x_train, x_test)
assert_true(all(is.finite(ranked$x_train)), "El rank transform produjo no finitos.")
assert_true(
  all(diff(ranked$x_train[order(x_train[, 1]), 1]) >= 0),
  "El rank transform no preserva el orden."
)

raw <- cbind(c(-1, -0.5, 0), c(0, 0.5, 1), c(1, 1.5, 2))
y_history <- c(-2, -1, 0)
offsets <- calibration_offsets(y_history, raw, c(0.05, 0.5, 0.95), min_obs = 3)
calibrated <- apply_quantile_calibration(raw, offsets, c(0.05, 0.5, 0.95))
assert_true(all(calibrated[, 1] <= calibrated[, 2]), "Cruce tras calibracion.")
assert_true(all(calibrated[, 2] <= calibrated[, 3]), "Cruce tras calibracion.")

set.seed(117)
misspecified_hits <- stats::rbinom(1000, size = 1L, prob = 0.15)
dq_misspecified <- dq_test(
  y = ifelse(misspecified_hits == 1L, -1, 1),
  q = rep(0, length(misspecified_hits)),
  tau = 0.05
)
assert_true(
  is.finite(dq_misspecified["p_value"]) &&
    dq_misspecified["p_value"] < 0.01,
  "El test DQ no detecta una cobertura incondicional incorrecta."
)

if (qbart_engine_available("official")) {
  set.seed(91)
  x <- matrix(runif(240), nrow = 120, ncol = 2)
  y <- sin(2 * pi * x[, 1]) + stats::rnorm(120, sd = 0.2)
  hp_smoke <- utils::modifyList(
    hp_qbart_default,
    list(m = 5L, burn = 10L, ndpost = 20L, nc = 20L, min_obs = 3L)
  )
  prediction <- predecir_qbart_grid(
    x[1:100, , drop = FALSE],
    y[1:100],
    x[101:120, , drop = FALSE],
    c(0.05, 0.5, 0.95),
    hp_smoke
  )
  assert_true(identical(dim(prediction), c(20L, 3L)), "Dimension QBART incorrecta.")
  assert_true(all(is.finite(prediction)), "QBART produjo predicciones no finitas.")
  assert_true(all(prediction[, 1] <= prediction[, 2]), "QBART produjo cruces.")
  assert_true(all(prediction[, 2] <= prediction[, 3]), "QBART produjo cruces.")
}

cat("OK: componentes QBART\n")
