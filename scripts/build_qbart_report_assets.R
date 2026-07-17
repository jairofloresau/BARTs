############################################################################
## Figuras y tablas reproducibles para el informe exploratorio QBART v2.
############################################################################

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1L]))
} else {
  normalizePath("scripts/build_qbart_report_assets.R")
}
root <- dirname(dirname(script_file))
source(file.path(root, "src", "qbart_features.R"))
source(file.path(root, "src", "qbart_workflow.R"))

processed <- file.path(root, "data", "processed")
figures <- file.path(root, "reports", "figures")
tables <- file.path(root, "reports", "tables")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)

score <- read.csv(
  file.path(processed, "qbart_quantile_score_exploracion_v2.csv"),
  check.names = FALSE
)
coverage <- read.csv(
  file.path(processed, "qbart_cobertura_exploracion_v2.csv"),
  check.names = FALSE
)
calibration <- read.csv(
  file.path(processed, "qbart_calibracion_exploracion_v2.csv"),
  check.names = FALSE
)
importance <- read.csv(
  file.path(processed, "qbart_importancia_por_tau_exploracion_v2.csv"),
  check.names = FALSE
)
simulation <- read.csv(
  file.path(root, "outputs", "qbart_simulation_gate_full.csv"),
  check.names = FALSE
)
checkpoint <- readRDS(
  file.path(processed, "checkpoint_qbart_exploracion_v2.rds")
)
base <- readRDS(
  file.path(root, "data", "datos", "BASE DIARIO", "base_final_diaria.rds")
)
model_data <- build_qbart_dataset(base, 1L)$data
test_index <- checkpoint$config$index_test
y_test <- model_data$y_objetivo[test_index]
test_dates <- as.Date(model_data$fecha[test_index])
predictions <- checkpoint$predictions
taus <- checkpoint$config$taus

method_labels <- c(
  random_walk = "Random Walk",
  quantile_reg = "Regresion cuantílica",
  quantile_rf = "QRF",
  qbart_raw = "QBART crudo",
  qbart_calibrated = "QBART calibrado"
)
method_colors <- c(
  "Random Walk" = "#6C757D",
  "Regresion cuantílica" = "#4F772D",
  "QRF" = "#0D3B66",
  "QBART crudo" = "#A23B3B",
  "QBART calibrado" = "#E09F3E"
)
theme_report <- theme_minimal(base_family = "serif", base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "#4A4A4A"),
    legend.position = "bottom"
  )

save_plot <- function(plot, filename, width = 7.2, height = 4.4) {
  ggsave(
    file.path(figures, filename),
    plot,
    width = width,
    height = height,
    device = "pdf"
  )
}

# Scores por percentil.
score_long <- score |>
  pivot_longer(-tau, names_to = "method", values_to = "pinball") |>
  mutate(
    method = unname(method_labels[method]),
    tau = factor(tau, levels = taus, labels = c("5", "50", "95"))
  )
score_plot <- ggplot(score_long, aes(tau, pinball, fill = method)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.72) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Pérdida pinball fuera de muestra",
    subtitle = "Menor valor implica mejor pronóstico cuantílico",
    x = "Percentil",
    y = "Pérdida pinball",
    fill = NULL
  ) +
  theme_report
save_plot(score_plot, "qbart_report_scores.pdf")

# Cobertura y objetivo nominal.
coverage_long <- coverage |>
  pivot_longer(-tau, names_to = "method", values_to = "coverage") |>
  mutate(
    method = unname(method_labels[method]),
    percentile = factor(tau, levels = taus, labels = c("5", "50", "95"))
  )
coverage_plot <- ggplot(
  coverage_long,
  aes(method, coverage, color = method)
) +
  geom_hline(
    aes(yintercept = tau),
    color = "#202020",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_point(size = 2.6) +
  facet_wrap(~percentile, scales = "free_y") +
  scale_color_manual(values = method_colors) +
  labs(
    title = "Cobertura empírica por percentil",
    subtitle = "La línea discontinua corresponde a la cobertura nominal",
    x = NULL,
    y = "Proporción observada",
    color = NULL
  ) +
  theme_report +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
save_plot(coverage_plot, "qbart_report_coverage.pdf")

# Resúmenes anuales y por bloque.
point_loss <- function(prediction, tau) {
  error <- y_test - prediction
  error * (tau - as.numeric(error < 0))
}
yearly_rows <- list()
block_rows <- list()
year <- format(test_dates, "%Y")
block <- ceiling(seq_along(y_test) / 60L)
row_index <- 0L
block_index <- 0L
for (j in seq_along(taus)) {
  qrf_loss <- point_loss(predictions$quantile_rf[, j], taus[j])
  qbart_loss <- point_loss(predictions$qbart_calibrated[, j], taus[j])
  for (current_year in unique(year)) {
    selected <- year == current_year
    row_index <- row_index + 1L
    yearly_rows[[row_index]] <- data.frame(
      year = as.integer(current_year),
      tau = taus[j],
      n = sum(selected),
      qrf = mean(qrf_loss[selected]),
      qbart = mean(qbart_loss[selected]),
      ratio = mean(qbart_loss[selected]) / mean(qrf_loss[selected]),
      coverage_qrf = mean(
        y_test[selected] <= predictions$quantile_rf[selected, j]
      ),
      coverage_qbart = mean(
        y_test[selected] <= predictions$qbart_calibrated[selected, j]
      )
    )
  }
  for (current_block in unique(block)) {
    selected <- block == current_block
    block_index <- block_index + 1L
    block_rows[[block_index]] <- data.frame(
      block = current_block,
      date = min(test_dates[selected]),
      tau = taus[j],
      qrf = mean(qrf_loss[selected]),
      qbart = mean(qbart_loss[selected]),
      loss_difference = mean(qbart_loss[selected] - qrf_loss[selected])
    )
  }
}
yearly <- do.call(rbind, yearly_rows)
block_summary <- do.call(rbind, block_rows)
write.csv(
  yearly,
  file.path(tables, "qbart_yearly_metrics.csv"),
  row.names = FALSE
)
write.csv(
  block_summary,
  file.path(tables, "qbart_block_metrics.csv"),
  row.names = FALSE
)

yearly_plot <- yearly |>
  mutate(percentile = factor(tau, taus, c("5", "50", "95"))) |>
  ggplot(aes(year, ratio, color = percentile, group = percentile)) +
  geom_hline(yintercept = 1, color = "#202020", linetype = "dashed") +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  scale_color_manual(values = c("5" = "#A23B3B", "50" = "#6C757D", "95" = "#0D3B66")) +
  scale_x_continuous(breaks = unique(yearly$year)) +
  labs(
    title = "Desempeño relativo anual de QBART",
    subtitle = "Ratio de pérdida QBART/QRF; valores menores que uno favorecen a QBART",
    x = "Año",
    y = "QBART calibrado / QRF",
    color = "Percentil"
  ) +
  theme_report
save_plot(yearly_plot, "qbart_report_yearly_ratio.pdf")

block_plot <- block_summary |>
  mutate(percentile = factor(tau, taus, c("5", "50", "95"))) |>
  ggplot(aes(date, loss_difference, fill = loss_difference > 0)) +
  geom_hline(yintercept = 0, color = "#202020", linewidth = 0.45) +
  geom_col(width = 65) +
  facet_wrap(~percentile, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = c("TRUE" = "#C65D3B", "FALSE" = "#2A7F62"),
    labels = c("TRUE" = "QRF gana", "FALSE" = "QBART gana")
  ) +
  labs(
    title = "Diferencia de pérdida por bloque",
    subtitle = "Positivo: QBART tiene mayor pérdida; cada bloque contiene hasta 60 observaciones",
    x = NULL,
    y = "Pinball QBART - QRF",
    fill = NULL
  ) +
  theme_report
save_plot(block_plot, "qbart_report_block_difference.pdf", height = 6.4)

# Correcciones de calibración.
calibration_plot <- calibration |>
  mutate(
    date = as.Date(date),
    percentile = factor(tau, taus, c("5", "50", "95"))
  ) |>
  ggplot(aes(date, offset, color = percentile)) +
  geom_hline(yintercept = 0, color = "#202020", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  facet_wrap(~percentile, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("5" = "#A23B3B", "50" = "#6C757D", "95" = "#0D3B66")) +
  labs(
    title = "Corrección causal de calibración",
    subtitle = "Cuantil móvil de los errores de pronóstico disponibles en cada origen",
    x = NULL,
    y = "Offset aditivo",
    color = "Percentil"
  ) +
  theme_report
save_plot(calibration_plot, "qbart_report_calibration.pdf", height = 6.2)

# Importancia exploratoria: nacimientos de ramas aceptados.
importance_long <- importance |>
  select(variable, all_of(c("tau_0.05", "tau_0.5", "tau_0.95"))) |>
  pivot_longer(-variable, names_to = "tau_name", values_to = "importance") |>
  mutate(
    percentile = recode(
      tau_name,
      "tau_0.05" = "5",
      "tau_0.5" = "50",
      "tau_0.95" = "95"
    )
  )
top_variables <- importance_long |>
  group_by(percentile) |>
  slice_max(importance, n = 10, with_ties = FALSE) |>
  pull(variable) |>
  unique()
importance_selected <- importance_long |>
  filter(variable %in% top_variables) |>
  mutate(
    variable = factor(
      variable,
      levels = rev(
        importance_long |>
          filter(variable %in% top_variables) |>
          group_by(variable) |>
          summarise(maximum = max(importance), .groups = "drop") |>
          arrange(maximum) |>
          pull(variable)
      )
    ),
    percentile = factor(percentile, c("5", "50", "95"))
  )
importance_heatmap <- ggplot(
  importance_selected,
  aes(percentile, variable, fill = importance)
) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#F4EFE6", high = "#8F1D2C") +
  labs(
    title = "Uso exploratorio de predictores por percentil",
    subtitle = "Proporción de nacimientos de ramas aceptados, incluido el burn-in",
    x = "Percentil",
    y = NULL,
    fill = "Frecuencia"
  ) +
  theme_report +
  theme(legend.position = "right")
save_plot(
  importance_heatmap,
  "qbart_report_importance_heatmap.pdf",
  width = 7.3,
  height = 7.8
)

importance_bars <- importance_selected |>
  group_by(percentile) |>
  slice_max(importance, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(variable = reorder(variable, importance)) |>
  ggplot(aes(importance, variable, fill = percentile)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~percentile, scales = "free_y") +
  scale_fill_manual(values = c("5" = "#A23B3B", "50" = "#6C757D", "95" = "#0D3B66")) +
  labs(
    title = "Principales señales por percentil",
    subtitle = "La medida es exploratoria y no equivale a una probabilidad posterior de inclusión",
    x = "Frecuencia relativa",
    y = NULL
  ) +
  theme_report
save_plot(importance_bars, "qbart_report_importance_top.pdf", width = 8.2, height = 5.4)

write.csv(
  importance_selected,
  file.path(tables, "qbart_importance_selected.csv"),
  row.names = FALSE
)

# Simulation gate.
simulation_plot_data <- simulation |>
  filter(method %in% c("qbart_raw", "qbart_calibrated", "qrf")) |>
  mutate(
    method = recode(
      method,
      qbart_raw = "QBART crudo",
      qbart_calibrated = "QBART calibrado",
      qrf = "QRF"
    ),
    percentile = factor(tau, taus, c("5", "50", "95"))
  )
simulation_plot <- ggplot(
  simulation_plot_data,
  aes(percentile, coverage, color = method, group = method)
) +
  geom_point(size = 2.4) +
  geom_line(linewidth = 0.8) +
  geom_point(
    aes(x = percentile, y = tau),
    inherit.aes = FALSE,
    data = data.frame(
      percentile = factor(c("5", "50", "95"), c("5", "50", "95")),
      tau = taus
    ),
    shape = 4,
    size = 3,
    stroke = 1
  ) +
  scale_color_manual(values = method_colors) +
  labs(
    title = "Simulation gate: cobertura",
    subtitle = "Las cruces indican la cobertura nominal conocida en la simulación",
    x = "Percentil",
    y = "Cobertura",
    color = NULL
  ) +
  theme_report
save_plot(simulation_plot, "qbart_report_simulation.pdf")

# Metadatos de la muestra para el documento.
n_total <- nrow(model_data)
train_end <- floor(0.60 * n_total)
validation_end <- floor(0.80 * n_total)
sample_summary <- data.frame(
  sample = c("Total", "Train inicial", "Validacion", "Test"),
  n = c(
    n_total,
    train_end,
    validation_end - train_end,
    n_total - validation_end
  ),
  start = as.character(c(
    min(model_data$fecha),
    min(model_data$fecha[seq_len(train_end)]),
    min(model_data$fecha[(train_end + 1L):validation_end]),
    min(model_data$fecha[(validation_end + 1L):n_total])
  )),
  end = as.character(c(
    max(model_data$fecha),
    max(model_data$fecha[seq_len(train_end)]),
    max(model_data$fecha[(train_end + 1L):validation_end]),
    max(model_data$fecha[(validation_end + 1L):n_total])
  ))
)
write.csv(
  sample_summary,
  file.path(tables, "qbart_sample_summary.csv"),
  row.names = FALSE
)

cat("Activos del informe generados en reports/figures y reports/tables.\n")
