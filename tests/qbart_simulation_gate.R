args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1L]))
} else {
  normalizePath("tests/qbart_simulation_gate.R")
}
project_root <- dirname(dirname(script_file))
Sys.setenv(QBART_PROJECT_ROOT = project_root)
source(file.path(project_root, "notebooks", "qbart_adapter.R"))
source(file.path(project_root, "src", "qbart_workflow.R"))

if (!qbart_engine_available("official")) {
  stop("El simulation gate requiere el paquete oficial BayesQArt.")
}
if (!requireNamespace("quantregForest", quietly = TRUE)) {
  stop("El simulation gate requiere quantregForest.")
}

mode <- match.arg(
  Sys.getenv("QBART_SIM_MODE", unset = "quick"),
  c("quick", "full")
)
strict <- identical(Sys.getenv("QBART_SIM_STRICT", unset = "0"), "1")
set.seed(7319)

n_fit <- if (mode == "full") 400L else 200L
n_validation <- if (mode == "full") 400L else 200L
n_test <- if (mode == "full") 1200L else 500L
p <- 5L
n <- n_fit + n_validation + n_test
x <- matrix(runif(n * p), nrow = n, ncol = p)
colnames(x) <- paste0("x", seq_len(p))

conditional_mean <- function(z) {
  1.5 * sin(pi * z[, 1] * z[, 2]) +
    1.2 * (z[, 3] - 0.5)^2 +
    0.5 * z[, 5]
}
conditional_sd <- function(z) {
  0.12 + 0.65 * z[, 4]
}

mu <- conditional_mean(x)
sigma <- conditional_sd(x)
y <- mu + sigma * rnorm(n)
fit_index <- seq_len(n_fit)
validation_index <- seq.int(n_fit + 1L, n_fit + n_validation)
test_index <- seq.int(n_fit + n_validation + 1L, n)
taus <- c(0.05, 0.50, 0.95)

hp <- utils::modifyList(
  hp_qbart_default,
  if (mode == "full") {
    list(m = 100L, burn = 1000L, ndpost = 1500L, nc = 50L)
  } else {
    list(m = 20L, burn = 100L, ndpost = 150L, nc = 30L)
  }
)
all_prediction_index <- c(validation_index, test_index)
qbart_all <- predecir_qbart_grid(
  x[fit_index, , drop = FALSE],
  y[fit_index],
  x[all_prediction_index, , drop = FALSE],
  taus,
  hp
)
qbart_validation <- qbart_all[seq_along(validation_index), , drop = FALSE]
qbart_test <- qbart_all[
  length(validation_index) + seq_along(test_index),
  ,
  drop = FALSE
]
offsets <- calibration_offsets(
  y[validation_index],
  qbart_validation,
  taus,
  window = Inf,
  min_obs = 100L
)
qbart_calibrated <- apply_quantile_calibration(qbart_test, offsets, taus)

set.seed(7319)
qrf <- quantregForest::quantregForest(
  x = x[fit_index, , drop = FALSE],
  y = y[fit_index],
  ntree = if (mode == "full") 1000L else 300L,
  mtry = 3L,
  nodesize = 10L,
  keep.inbag = FALSE
)
qrf_test <- predict(qrf, newdata = x[test_index, , drop = FALSE], what = taus)
oracle <- vapply(
  taus,
  function(tau) mu[test_index] + sigma[test_index] * stats::qnorm(tau),
  numeric(length(test_index))
)

methods <- list(
  qbart_raw = qbart_test,
  qbart_calibrated = qbart_calibrated,
  qrf = qrf_test,
  oracle = oracle
)
result <- do.call(rbind, lapply(names(methods), function(method) {
  data.frame(
    mode = mode,
    method = method,
    tau = taus,
    pinball = quantile_scores(methods[[method]], y[test_index], taus),
    coverage = quantile_coverage(methods[[method]], y[test_index], taus)
  )
}))
rownames(result) <- NULL
print(result, row.names = FALSE, digits = 5)

output_dir <- file.path(project_root, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  result,
  file.path(output_dir, paste0("qbart_simulation_gate_", mode, ".csv")),
  row.names = FALSE
)

calibrated_result <- result[result$method == "qbart_calibrated", ]
qrf_result <- result[result$method == "qrf", ]
coverage_tolerance <- if (mode == "full") 0.035 else 0.075
coverage_pass <- all(
  abs(calibrated_result$coverage - calibrated_result$tau) <= coverage_tolerance
)
relative_pinball <- calibrated_result$pinball / qrf_result$pinball
tail_index <- abs(calibrated_result$tau - 0.5) > 0.25
tail_performance_pass <- all(relative_pinball[tail_index] <= 1.10)
center_ratio <- relative_pinball[which.min(abs(calibrated_result$tau - 0.5))]

cat(
  sprintf(
    paste0(
      "Gate cobertura: %s | Gate colas QBART/QRF <= 1.10: %s | ",
      "max ratio colas: %.3f | ratio mediana: %.3f\n"
    ),
    if (coverage_pass) "PASS" else "FAIL",
    if (tail_performance_pass) "PASS" else "FAIL",
    max(relative_pinball[tail_index]),
    center_ratio
  )
)
if (strict && !(coverage_pass && tail_performance_pass)) {
  stop("QBART aun no supera el simulation gate.", call. = FALSE)
}
