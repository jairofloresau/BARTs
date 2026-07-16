# ============================================================================
# TESIS: Pronostico e Identificacion de Determinantes de los Retornos
#        Extremos del Tipo de Cambio en Peru
#
# SCRIPT: QBART - Bayesian Quantile Additive Regression Trees
#         (reemplaza al BART heteroscedastico de "BART riesgo de cola 2.R"
#          por pedido del asesor: QBART estima cada cuantil DIRECTAMENTE,
#          en vez de reconstruir la distribucion via location-scale)
# ============================================================================
#
# QUE CAMBIA RESPECTO A "BART riesgo de cola 2.R" (resumen; el detalle de
# cada cambio esta comentado en el lugar donde ocurre):
#
#   1. Modelo principal: QBART (bayesqart.R / qbart_adapter.R) en vez de
#      wbart() del paquete BART con reconstruccion media+log-varianza.
#      QBART ajusta UN modelo POR CADA cuantil tau -> es matematicamente mas
#      directo (no asume Normalidad condicional y_i | x_i con solo sigma(x)
#      variando), pero mas caro: nada de "un ajuste, todos los tau"; cada tau
#      es su propio MCMC de arboles.
#   2. Solo 3 cuantiles: tau = 0.05, 0.50, 0.95 (antes: 7 tau + grilla fina de
#      99 para el qwCRPS). Por eso el qwCRPS integrado en grilla fina
#      desaparece: no corresponde con solo 3 puntos. La evaluacion queda en
#      Quantile Score (pinball) por tau + cobertura empirica.
#   3. Sin BART package: no hace falta reconstruir_draws()/blindar_x()/
#      ajustar_bart_het() del script original -> predecir_qbart_grid() ya
#      devuelve la matriz de cuantiles pedida y maneja el caso de 1 sola fila
#      de x_pred internamente (via vapply + coercion a matriz).
#   4. Ventana expansiva CON CHECKPOINTS: dado el costo de QBART (ver mas
#      abajo), la evaluacion recursiva puede tardar horas. Se guarda el
#      progreso por bloque en "checkpoint_ventana.rds" y se reanuda solo si
#      la configuracion (paso_reest, hp_qbart) coincide con la guardada.
#   5. Determinantes: ya NO hace falta el desglose media/varianza + PDP con
#      reconstruir_draws() del script original (era el sustituto de
#      "importancia por cuantil" para un modelo location-scale). QBART da
#      importancia POR CUANTIL de forma NATIVA (frecuencia de uso de cada
#      variable en los arboles de CADA tau), via importancia_qbart_por_cuantil().
#
# COSTO COMPUTACIONAL (motor R puro; no esta instalado el paquete Rcpp
# 'BayesQArt' en esta maquina, asi que qbart_adapter.R cae al fallback de
# bayesqart.R): un benchmark a la escala real del problema (n~6200,
# p=30, m=100 arboles) dio ~0.96 seg POR ITERACION MCMC. Con
# m=100/burn=1000/ndpost=2000 y paso_reest=60 (23 bloques x 3 tau x 3000
# iteraciones = 69 ajustes completos), la ventana expansiva sola es
# ~55 horas; sumando tuneo (~4h) y determinantes (~2h), la corrida COMPLETA
# es del orden de ~60 horas (~2.5 dias) corriendo sin parar. De ahi que:
#   (a) el checkpoint por bloque sea obligatorio, no un lujo, y
#   (b) el script arranque SIEMPRE en MODO_PRUEBA (parametros minimos) hasta
#       que se cambie manualmente el flag de abajo.
# ============================================================================


# ============================================================================
# 1. SETUP   [reutilizado tal cual de "BART riesgo de cola 2.R", con UN
#             cambio: se quita el paquete "BART" de la lista porque el
#             modelo principal ya no lo usa -sigue en la maquina para quien
#             quiera comparar contra el script original, pero este script no
#             lo necesita-]
# ============================================================================

paquetes <- c(
  "dplyr", "tidyr", "purrr", "ggplot2",   # manipulacion y graficos
  "quantreg",                              # Koenker y Bassett (1978) - benchmark
  "quantregForest"                         # Meinshausen (2006) - benchmark
)
instalar_faltantes <- paquetes[!paquetes %in% installed.packages()[, "Package"]]
if (length(instalar_faltantes) > 0) install.packages(instalar_faltantes)
invisible(lapply(paquetes, library, character.only = TRUE))

set.seed(20222834)  # codigo de la autora, para reproducibilidad del MCMC

ruta_datos  <- "D:/BARTs/data/datos/BASE DIARIO"
ruta_out    <- "D:/BARTs/data/processed"
ruta_figs   <- "D:/BARTs/reports/figures"
dir.create(ruta_out,  recursive = TRUE, showWarnings = FALSE)
dir.create(ruta_figs, recursive = TRUE, showWarnings = FALSE)

h_pronostico <- 1

# --- CAMBIO: taus se redefine a solo 3 cuantiles ---
# El script original usaba taus = c(.01,.05,.10,.50,.90,.95,.99) + una grilla
# fina de 99 puntos para integrar el qwCRPS. QBART ajusta un modelo MCMC
# COMPLETO por cada tau (no comparte arboles entre cuantiles como el BART
# location-scale), asi que agrandar el vector de tau multiplica el costo de
# TODO el script proporcionalmente. Por pedido explicito: solo 0.05/0.50/0.95.
taus <- c(0.05, 0.50, 0.95)

# --- CAMBIO: paso_reest ahora se decide mas abajo, segun MODO_PRUEBA (§6) ---
# (en el script original era una constante fija; aqui depende del modo)


# ============================================================================
# 2. CARGAR BASE Y CONSTRUIR VARIABLE DEPENDIENTE + PREDICTORES
#    [reutilizado tal cual de "BART riesgo de cola 2.R", sin cambios]
# ============================================================================

base <- readRDS(file.path(ruta_datos, "base_final_diaria.rds")) %>%
  arrange(fecha)

base <- base %>% mutate(retorno_tc = 100 * retorno_tc)

retorno_futuro <- function(retorno_diario, h) {
  n <- length(retorno_diario)
  out <- rep(NA_real_, n)
  for (t in seq_len(n - h)) {
    out[t] <- sum(retorno_diario[(t + 1):(t + h)])
  }
  out
}

base <- base %>%
  mutate(y_objetivo = retorno_futuro(retorno_tc, h_pronostico))

variables_nivel <- c(
  "vix", "baa10y",
  "tasa_usa", "nfci", "yield_to_worst",
  "tp_1y", "tp_5y", "tp_10y",
  "embig_peru"
)

variables_precio <- c(
  "rin",
  "oro", "cobre", "zinc", "plata", "wti",
  "spx", "mxef",
  "dxy", "jpy", "chf", "eur", "gbp",
  "clp", "cop", "mxn", "brl", "cny"
)

base <- base %>%
  mutate(across(all_of(variables_precio), ~ 100 * c(NA, diff(log(.x))),
                .names = "d_{.col}")) %>%
  mutate(
    diff_tasas = tasa_peru - tasa_usa,
    interv_bcrp = ifelse(is.na(interv_bcrp), 0, interv_bcrp)
  )

predictores <- c(
  variables_nivel,
  "diff_tasas", "interv_bcrp",
  paste0("d_", variables_precio),
  "retorno_tc"
)

datos_modelo <- base %>%
  select(fecha, y_objetivo, all_of(predictores)) %>%
  na.omit()

cat("Observaciones disponibles:", nrow(datos_modelo), "\n")
cat("Predictores:", length(predictores), "\n")
cat("Rango de fechas:", as.character(min(datos_modelo$fecha)), "a",
    as.character(max(datos_modelo$fecha)), "\n\n")


# ============================================================================
# 3. SPLIT CRONOLOGICO EN TRES BLOQUES: TRAIN / VALIDACION / TEST
#    [reutilizado tal cual de "BART riesgo de cola 2.R", sin cambios]
# ============================================================================

n_total   <- nrow(datos_modelo)
corte_tr  <- floor(0.60 * n_total)
corte_val <- floor(0.80 * n_total)

idx_train <- 1:corte_tr
idx_val   <- (corte_tr + 1):corte_val
idx_test  <- (corte_val + 1):n_total

x_train <- as.matrix(datos_modelo[idx_train, predictores])
y_train <- datos_modelo$y_objetivo[idx_train]
x_val   <- as.matrix(datos_modelo[idx_val, predictores])
y_val   <- datos_modelo$y_objetivo[idx_val]

cat("Train:", length(idx_train), "obs. |",
    "Validacion (tuneo):", length(idx_val), "obs. |",
    "Test (evaluacion final):", length(idx_test), "obs.\n\n")

formula_rq <- as.formula(paste("y_objetivo ~", paste(predictores, collapse = " + ")))


# ============================================================================
# 4. METRICA: QUANTILE SCORE (PINBALL)
#    [reutilizado tal cual de "BART riesgo de cola 2.R"; se elimina el
#     qwCRPS integrado en grilla fina -no aplica con solo 3 tau puntuales,
#     ver nota en §1-]
# ============================================================================

pinball <- function(y_real, y_pred, tau) {
  u <- y_real - y_pred
  mean(u * (tau - as.numeric(u < 0)))
}

qs_metodo <- function(pred_matrix, y_real, taus) {
  sapply(seq_along(taus), function(i) pinball(y_real, pred_matrix[, i], taus[i]))
}


# ============================================================================
# 5. MOTOR QBART
# ============================================================================
# CAMBIO: reemplaza por completo a la §5 del script original (ajustar_bart_het
# + reconstruir_draws + blindar_x). Ya no hace falta ningun helper propio:
# predecir_qbart_grid() e importancia_qbart_por_cuantil() (definidas en
# qbart_adapter.R) devuelven directamente lo que el resto del script necesita.
#
# NOTA: tus archivos estan en D:/BARTs/notebooks/, no en D:/BARTs/ como en tu
# mensaje -se usan rutas absolutas para que no importe el working directory-.
#
# CAMBIO: se fuerza el motor en R puro (FORZAR_MOTOR_R <- TRUE) ANTES de
# source(qbart_adapter.R). Se instalo el paquete Rcpp 'BayesQArt' via Rtools,
# pero un benchmark a la escala real de este problema (n~6200, p=30, m=100)
# mostro que resulta ~1.7x MAS LENTO por iteracion que bayesqart.R (0.82 vs
# 0.43 s/iter; bajar 'nc' no lo mejora, el costo no esta ahi). Sin este
# interruptor, qbart_adapter.R elegiria el paquete compilado automaticamente
# por estar instalado -y seria la opcion mas lenta, no la mas rapida-.
FORZAR_MOTOR_R <- TRUE
source("D:/BARTs/notebooks/bayesqart.R")
source("D:/BARTs/notebooks/qbart_adapter.R")


# ============================================================================
# 6. MODO DE CORRIDA: prueba / exploracion / completa
# ============================================================================
# CAMBIO: se agrega un tercer modo, "exploracion", entre el smoke test
# (MODO_PRUEBA original) y la corrida completa (~60h). Con el motor en R puro
# confirmado como el mas rapido, "exploracion" usa la MISMA estructura de
# bloques que la corrida completa (paso_reest = 60, los 23 bloques reales) -
# para que el resultado sea directamente comparable/informativo- pero con
# menos arboles (m=50) y cadenas MCMC mas cortas (burn=300, ndpost=500).
# Estimado: ~26 min de tuneo + ~3.3h de ventana expansiva + ~8 min de
# determinantes =~ 3.9h, con margen bajo las 6h pedidas.
MODO <- "exploracion"   # "prueba" | "exploracion" | "completa"

if (MODO == "prueba") {
  hp_qbart          <- modifyList(hp_qbart_default,
                                  list(m = 15L, burn = 150L, ndpost = 200L, kappa = 2))
  paso_reest        <- 200L
  burn_tune         <- 15L
  ndpost_tune       <- 20L
  grilla_qbart      <- expand.grid(kappa = 2, m = 15L)   # grilla achicada para el smoke test
  cat("*** MODO = 'prueba': hiperparametros minimos, solo para validar el flujo ***\n")
  cat("*** Los resultados de esta corrida NO son la version final ***\n\n")
} else if (MODO == "exploracion") {
  hp_qbart          <- modifyList(hp_qbart_default,
                                  list(m = 50L, burn = 300L, ndpost = 500L, kappa = 2))
  paso_reest        <- 60L   # misma granularidad de bloques que la corrida completa
  burn_tune         <- 150L
  ndpost_tune       <- 250L
  grilla_qbart      <- expand.grid(kappa = c(1, 2, 3), m = c(50, 100))  # 6 configs
  cat("*** MODO = 'exploracion': m=50, burn=300, ndpost=500, paso_reest=60 ***\n")
  cat("*** Estimado ~3.9h. Resultados indicativos, no la version final. ***\n\n")
} else {
  hp_qbart          <- modifyList(hp_qbart_default,
                                  list(m = 100L, burn = 1000L, ndpost = 2000L, kappa = 2))
  paso_reest        <- 60L
  burn_tune         <- 250L
  ndpost_tune       <- 500L
  grilla_qbart      <- expand.grid(kappa = c(1, 2, 3), m = c(50, 100, 200))
}

sufijo_out <- paste0("_", MODO)


# ============================================================================
# 7. TUNEO DE HIPERPARAMETROS: GRILLA kappa x m
# ============================================================================
# CAMBIO: reemplaza a la §6 del script original (grilla k x ntree x sparse
# para wbart). Aqui: kappa (~ shrinkage, analogo a "k" en BART) x m (numero
# de arboles, analogo a "ntree"). No hay analogo directo a "sparse" (DART) en
# QBART tal como esta implementado -ver nota al final de qbart_adapter.R-.
# Seleccion: pinball PROMEDIO de los 3 tau sobre la validacion (no hay qwCRPS
# integrado con solo 3 puntos, asi que el promedio simple es el criterio).

cat("=== TUNEO QBART: evaluando", nrow(grilla_qbart), "configuraciones sobre la validacion ===\n\n")

resultados_tuneo <- vector("list", nrow(grilla_qbart))

for (i in seq_len(nrow(grilla_qbart))) {
  cfg <- grilla_qbart[i, ]
  t0  <- Sys.time()

  hp_cfg <- modifyList(hp_qbart_default,
                       list(kappa = cfg$kappa, m = cfg$m,
                            burn = burn_tune, ndpost = ndpost_tune))

  pred_val <- predecir_qbart_grid(x_train, y_train, x_val, taus, hp_cfg)

  qs_val <- qs_metodo(pred_val, y_val, taus)
  pinball_prom <- mean(qs_val)

  resultados_tuneo[[i]] <- data.frame(
    kappa = cfg$kappa, m = cfg$m,
    pinball_prom = pinball_prom,
    setNames(as.list(qs_val), paste0("pinball_tau_", taus)),
    minutos = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  )

  cat(sprintf("[%2d/%2d] kappa=%g  m=%3d  ->  pinball_prom=%.5f  (%.1f min)\n",
              i, nrow(grilla_qbart), cfg$kappa, cfg$m, pinball_prom,
              resultados_tuneo[[i]]$minutos))
}

tabla_tuneo_qbart <- do.call(rbind, resultados_tuneo)
tabla_tuneo_qbart <- tabla_tuneo_qbart[order(tabla_tuneo_qbart$pinball_prom), ]

cat("\n--- Resultados del tuneo QBART (ordenados por pinball promedio; menor = mejor) ---\n")
print(round(tabla_tuneo_qbart, 5), row.names = FALSE)

write.csv(tabla_tuneo_qbart,
          file.path(ruta_out, paste0("qbart_tuneo_grilla", sufijo_out, ".csv")),
          row.names = FALSE)

mejor <- tabla_tuneo_qbart[1, ]
kappa_opt <- mejor$kappa
m_opt     <- mejor$m

cat("\n============================================================\n")
cat("CONFIGURACION SELECCIONADA: kappa =", kappa_opt, "| m =", m_opt, "\n")
cat("============================================================\n\n")

# El resto del script usa hp_qbart con kappa/m TUNEADOS pero burn/ndpost
# COMPLETOS (los de §6), igual que el script original reestima con la cadena
# completa tras el tuneo.
hp_qbart <- modifyList(hp_qbart, list(kappa = kappa_opt, m = m_opt))


# ============================================================================
# 8. EVALUACION FUERA DE MUESTRA CON VENTANA EXPANSIVA (RECURSIVA POR BLOQUES,
#    CON CHECKPOINTS)
# ============================================================================
# CAMBIO: reemplaza a la §7 del script "BART riesgo de cola 2.R" (que ya
# habiamos optimizado a bloques). Aqui se agrega, ademas, CHECKPOINTING por
# bloque: dado que un solo bloque completo (3 tau x 3000 iteraciones) puede
# tardar ~40 minutos en modo completo, no reanudar desde cero ante un corte
# es indispensable. El checkpoint guarda tambien la configuracion usada
# (paso_reest, hp_qbart) y se descarta -no se reutiliza- si no coincide con
# la configuracion actual, para no mezclar una corrida de prueba con una
# completa (o dos configuraciones de tuneo distintas) a medio camino.

origenes <- idx_test
n_or  <- length(origenes)
n_tau <- length(taus)
n_bloques <- ceiling(n_or / paso_reest)

archivo_checkpoint <- file.path(ruta_out, "checkpoint_ventana.rds")

config_actual <- list(paso_reest = paso_reest, hp_qbart = hp_qbart, taus = taus,
                      modo = MODO)

bloque_inicial <- 1L
if (file.exists(archivo_checkpoint)) {
  ckpt <- readRDS(archivo_checkpoint)
  mismo_config <- isTRUE(all.equal(ckpt$config, config_actual))
  if (mismo_config) {
    pred_fina      <- ckpt$pred_fina
    bloque_inicial <- ckpt$ultimo_bloque + 1L
    cat("=== Checkpoint encontrado y compatible: reanudando desde el bloque",
        bloque_inicial, "de", n_bloques, "===\n\n")
  } else {
    cat("=== Checkpoint encontrado pero con OTRA configuracion (paso_reest/hp_qbart/",
        "modo distinto) -> se descarta y se arranca de cero ===\n\n")
  }
}

if (bloque_inicial == 1L) {
  pred_fina <- list(
    random_walk  = matrix(NA_real_, n_or, n_tau),
    quantile_reg = matrix(NA_real_, n_or, n_tau),
    quantile_rf  = matrix(NA_real_, n_or, n_tau),
    qbart        = matrix(NA_real_, n_or, n_tau)
  )
}

y_real <- datos_modelo$y_objetivo[origenes]

cat("=== EVALUACION RECURSIVA POR BLOQUES (QBART):", n_or, "origenes,",
    n_bloques, "bloques de hasta", paso_reest, "obs. ===\n\n")

t_inicio_eval <- Sys.time()

if (bloque_inicial <= n_bloques) {
  for (b in bloque_inicial:n_bloques) {

    idx_bloque <- ((b - 1) * paso_reest + 1):min(b * paso_reest, n_or)
    t_bloque   <- origenes[idx_bloque]
    t_ini      <- min(t_bloque)

    idx_hasta <- 1:(t_ini - 1)
    df_hasta  <- datos_modelo[idx_hasta, ]
    x_hasta   <- as.matrix(df_hasta[, predictores])
    y_hasta   <- df_hasta$y_objetivo

    df_b <- datos_modelo[t_bloque, , drop = FALSE]
    x_b  <- as.matrix(df_b[, predictores])

    # -- BENCHMARK 1: Random Walk -- (identico al script original)
    cuantiles_rw <- quantile(y_hasta, probs = taus, na.rm = TRUE)

    # -- BENCHMARK 2: Regresion cuantilica lineal --
    # Cruce de cuantiles posible (cada tau es un modelo separado en rq) ->
    # se reordena por fila (Chernozhukov et al., 2010), como en el original.
    modelo_rq  <- rq(formula_rq, tau = taus, data = df_hasta)
    pred_rq_b  <- predict(modelo_rq, newdata = df_b)
    pred_rq_b  <- t(apply(pred_rq_b, 1, sort))

    # -- BENCHMARK 3: Quantile Random Forest --
    modelo_qrf  <- quantregForest(x = x_hasta, y = y_hasta, ntree = 200)
    pred_qrf_b  <- predict(modelo_qrf, newdata = x_b, what = taus)

    # -- MODELO PRINCIPAL: QBART --
    # predecir_qbart_grid() ajusta 3 modelos QBART (uno por tau) sobre la
    # ventana x_hasta/y_hasta y predice el bloque x_b de una sola vez por tau.
    # Ya viene con el anti-cruce aplicado (sort por fila) dentro del
    # adaptador; se re-aplica aqui solo para dejar explicito -por pedido- el
    # mismo tratamiento que a rq (es un sort() sobre datos ya ordenados: no
    # cambia nada, pero documenta la intencion).
    pred_qbart_b <- predecir_qbart_grid(x_hasta, y_hasta, x_b, taus, hp_qbart)
    pred_qbart_b <- t(apply(pred_qbart_b, 1, sort))

    pred_fina$random_walk[idx_bloque, ]  <- matrix(cuantiles_rw, nrow = length(idx_bloque),
                                                    ncol = n_tau, byrow = TRUE)
    pred_fina$quantile_reg[idx_bloque, ] <- pred_rq_b
    pred_fina$quantile_rf[idx_bloque, ]  <- pred_qrf_b
    pred_fina$qbart[idx_bloque, ]        <- pred_qbart_b

    transcurrido_min <- as.numeric(difftime(Sys.time(), t_inicio_eval, units = "mins"))
    n_hechos  <- b - bloque_inicial + 1
    restante_min <- transcurrido_min / n_hechos * (n_bloques - b)
    cat(sprintf("[bloque %2d/%2d] origenes %4d-%4d  (ventana hasta %s)  ->  %.1f min transcurridos, ~%.1f min restantes\n",
                b, n_bloques, min(idx_bloque), max(idx_bloque),
                as.character(df_hasta$fecha[nrow(df_hasta)]),
                transcurrido_min, restante_min))

    # --- checkpoint tras cada bloque ---
    saveRDS(list(pred_fina = pred_fina, ultimo_bloque = b, config = config_actual),
            archivo_checkpoint)
  }
}

cat("\nEvaluacion recursiva QBART completada. Tiempo total:",
    round(as.numeric(difftime(Sys.time(), t_inicio_eval, units = "mins")), 1), "min\n\n")


# ============================================================================
# 9. RESULTADOS DE LA EVALUACION
# ============================================================================
# CAMBIO: reemplaza a la §8 del script original. Se elimina el qwCRPS (no
# aplica a 3 tau puntuales) y se AGREGA cobertura empirica -pedida
# explicitamente-: fraccion de retornos reales por debajo del cuantil
# estimado; debe acercarse a tau (0.05, 0.50, 0.95) si el modelo esta
# bien calibrado.

## ---- 9.1 Quantile Score por metodo y cuantil ----
tabla_quantile_score <- data.frame(
  tau          = taus,
  random_walk  = qs_metodo(pred_fina$random_walk,  y_real, taus),
  quantile_reg = qs_metodo(pred_fina$quantile_reg, y_real, taus),
  quantile_rf  = qs_metodo(pred_fina$quantile_rf,  y_real, taus),
  qbart        = qs_metodo(pred_fina$qbart,        y_real, taus)
)

cat("--- Quantile Score por metodo y cuantil (menor = mejor) ---\n")
print(round(tabla_quantile_score, 5), row.names = FALSE)

write.csv(tabla_quantile_score,
          file.path(ruta_out, paste0("qbart_quantile_score_comparacion", sufijo_out, ".csv")),
          row.names = FALSE)

## ---- 9.2 Cobertura empirica ----
cobertura <- function(pred_matrix, y_real, taus) {
  sapply(seq_along(taus), function(i) mean(y_real <= pred_matrix[, i]))
}

tabla_cobertura <- data.frame(
  tau          = taus,
  random_walk  = cobertura(pred_fina$random_walk,  y_real, taus),
  quantile_reg = cobertura(pred_fina$quantile_reg, y_real, taus),
  quantile_rf  = cobertura(pred_fina$quantile_rf,  y_real, taus),
  qbart        = cobertura(pred_fina$qbart,        y_real, taus)
)

cat("\n--- Cobertura empirica por metodo y cuantil (debe acercarse a 'tau') ---\n")
print(round(tabla_cobertura, 4), row.names = FALSE)

write.csv(tabla_cobertura,
          file.path(ruta_out, paste0("qbart_cobertura_empirica", sufijo_out, ".csv")),
          row.names = FALSE)

## ---- 9.3 Desempenio relativo al Random Walk ----
ratios <- tabla_quantile_score
ratios[, -1] <- sweep(ratios[, -1], 1, ratios$random_walk, "/")
cat("\n--- Quantile Score RELATIVO al Random Walk (<1 => mejor que el RW) ---\n")
print(round(ratios, 4), row.names = FALSE)


# ============================================================================
# 10. DETERMINANTES: IMPORTANCIA POR CUANTIL (NATIVA EN QBART)
# ============================================================================
# CAMBIO: reemplaza a las §9 y §10 completas del script original (el
# desglose media/varianza + la dependencia parcial bayesiana con
# reconstruir_draws()). QBART da la importancia por cuantil DIRECTAMENTE:
# como cada tau es su propio arbol-suma, la frecuencia de uso de una variable
# en los arboles de tau=0.05 vs. tau=0.50 vs. tau=0.95 YA ES la respuesta a
# "¿los determinantes de las colas difieren de los del centro?", sin
# necesidad de dependencia parcial simulada.

idx_det <- 1:corte_val
x_det   <- as.matrix(datos_modelo[idx_det, predictores])
y_det   <- datos_modelo$y_objetivo[idx_det]

cat("\n=== Determinantes QBART (importancia por cuantil), ", length(idx_det),
    "obs. ===\n")

tabla_importancia_qbart <- importancia_qbart_por_cuantil(x_det, y_det, taus, hp_qbart)

cat("\n--- Importancia por cuantil (frecuencia de uso en los arboles de cada tau) ---\n")
print(head(tabla_importancia_qbart[order(-tabla_importancia_qbart[[paste0("tau_", taus[1])]]), ], 15),
      row.names = FALSE, digits = 3)

write.csv(tabla_importancia_qbart,
          file.path(ruta_out, paste0("qbart_importancia_por_cuantil", sufijo_out, ".csv")),
          row.names = FALSE)

## ---- 10.1 Heatmap de importancia por cuantil ----
cols_tau <- paste0("tau_", taus)
heat <- tabla_importancia_qbart %>%
  select(variable, all_of(cols_tau)) %>%
  pivot_longer(-variable, names_to = "cuantil", values_to = "importancia") %>%
  group_by(variable) %>%
  mutate(importancia_rel = importancia / max(importancia)) %>%
  ungroup()

g_heat <- ggplot(heat, aes(x = cuantil, y = reorder(variable, importancia_rel),
                           fill = importancia_rel)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = "QBART: importancia de cada predictor POR CUANTIL",
       subtitle = "Frecuencia de uso en los arboles, normalizada por variable",
       x = "Cuantil", y = NULL, fill = "Importancia\nrelativa") +
  theme_minimal()

ggsave(file.path(ruta_figs, paste0("qbart_importancia_por_cuantil_heatmap", sufijo_out, ".png")),
       g_heat, width = 8, height = 9)


# ============================================================================
# 11. RESUMEN FINAL
# ============================================================================
cat("\n============================================================\n")
cat("RESUMEN QBART (MODO:", toupper(MODO), ")\n")
cat("  Hiperparametros: kappa =", hp_qbart$kappa, "| m =", hp_qbart$m,
    "| burn =", hp_qbart$burn, "| ndpost =", hp_qbart$ndpost, "\n")
cat("  paso_reest =", paso_reest, "| Origenes evaluados:", n_or, "\n")
cat("  Mejor metodo por pinball promedio:",
    names(tabla_quantile_score)[-1][which.min(colMeans(tabla_quantile_score[-1]))], "\n")
cat("\nArchivos guardados en", ruta_out, ":\n")
cat("    - qbart_tuneo_grilla", sufijo_out, ".csv\n", sep = "")
cat("    - qbart_quantile_score_comparacion", sufijo_out, ".csv\n", sep = "")
cat("    - qbart_cobertura_empirica", sufijo_out, ".csv\n", sep = "")
cat("    - qbart_importancia_por_cuantil", sufijo_out, ".csv\n", sep = "")
cat("    - checkpoint_ventana.rds (progreso de la ventana expansiva)\n")
cat("Figura guardada en", ruta_figs, ":\n")
cat("    - qbart_importancia_por_cuantil_heatmap", sufijo_out, ".png\n", sep = "")
if (MODO != "completa") {
  cat("\n*** Esto fue MODO =", paste0("'", MODO, "'"), ". Cambia MODO <- \"completa\" en la",
      "§6 para la corrida completa (~60h estimadas; ver cabecera del script). ***\n")
}
cat("============================================================\n")
