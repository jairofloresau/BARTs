# ============================================================================
# TESIS: Pronostico e Identificacion de Determinantes de los Retornos
#        Extremos del Tipo de Cambio en Peru: Un Enfoque BART
#
# SCRIPT 02: Riesgo de cola con BART - distribucion completa y determinantes
# ============================================================================
#
# Que hace este script, en orden:
#
#   1. Carga "base_final_diaria" (construida en "Base de datos.R") y arma:
#        - la variable dependiente: retorno acumulado a h dias, r_(t+h)
#        - el set de predictores X_t con TODAS las variables disponibles en
#          la base (no solo las del Cuadro 2 del avance): ademas del riesgo
#          sistemico global, condiciones financieras EE.UU., fundamentales
#          macro Peru, materias primas e incertidumbre politica, se agregan
#          la prima por plazo (term premia 1/5/10 anios), el spread high
#          yield, plata, petroleo (WTI), indices bursatiles (S&P 500, MXEF)
#          y el resto de monedas frente al dolar (DXY y los cruces
#          individuales). BART no impone una forma funcional ni penaliza por
#          incluir variables adicionales de la misma forma que un modelo
#          lineal, por lo que ampliar el set de predictores permite que sea
#          el propio algoritmo -no la teoria previa- el que revele cuales
#          determinantes importan en cada cuantil.
#   2. Separa la muestra en entrenamiento y evaluacion fuera de muestra
#      (split cronologico, sin barajar, porque es una serie de tiempo)
#   3. Estima 3 metodos de referencia (benchmarks) del avance:
#        - Caminata aleatoria (Random Walk)
#        - Regresion cuantilica lineal (Koenker y Bassett, 1978)
#        - Quantile Random Forest (Meinshausen, 2006)
#   4. Estima BART (Chipman, George & McCulloch, 2010; paquete `BART`) y
#      reconstruye la distribucion predictiva posterior COMPLETA sumando
#      ruido con la sigma posterior a cada extraccion MCMC de la media
#      condicional, siguiendo a Clark, Huber, Koop, Marcellino & Pfarrhofer
#      (2023, "Tail forecasting with multivariate BART")
#   5. Extrae de esa distribucion posterior los cuantiles de interes
#      tau = {0.01, 0.05, 0.10, 0.90, 0.95, 0.99} definidos en el avance
#   6. Evalua el desempeno fuera de muestra de los 4 metodos con:
#        - Quantile Score (funcion pinball) por cuantil
#        - Quantile-weighted CRPS (ajuste ponderado a toda la distribucion)
#   7. Identifica los determinantes principales de cada cuantil:
#        a) importancia global de variables (frecuencia de uso en los
#           arboles, `varcount` de BART)
#        b) efecto marginal por cuantil (dependencia parcial bayesiana):
#           como se desplaza cada Q_tau cuando se mueve un predictor a lo
#           largo de su rango empirico, manteniendo los demas en su mediana
#
# NOTA METODOLOGICA: BART no necesita un modelo distinto por cuantil (a
# diferencia de la regresion cuantilica lineal y el QRF). Se estima UN solo
# modelo sobre la media condicional; todos los cuantiles y todos los efectos
# marginales por cuantil se derivan de la misma distribucion predictiva
# posterior. Esa es precisamente la ventaja que argumenta el avance frente
# a la regresion cuantilica lineal.
# ============================================================================


# ============================================================================
# 1. SETUP
# ============================================================================

paquetes <- c(
  "dplyr", "tidyr", "purrr", "ggplot2",   # manipulacion y graficos
  "BART",                                  # Chipman, George & McCulloch (2010)
  "quantreg",                              # Koenker y Bassett (1978)
  "quantregForest"                         # Meinshausen (2006)
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

# --- Horizonte de pronostico (en dias) y cuantiles de interes (avance, ec. 2) ---
h_pronostico <- 1
taus <- c(0.01, 0.05, 0.10, 0.90, 0.95, 0.99)

# --- Hiperparametros MCMC de BART (defaults de Chipman et al., 2010) ---
n_arboles <- 200   # m = 200, el default del paper fundacional
n_burn    <- 1000
n_post    <- 2000


# ============================================================================
# 2. CARGAR BASE Y CONSTRUIR VARIABLE DEPENDIENTE + PREDICTORES
# ============================================================================

base <- readRDS(file.path(ruta_datos, "base_final_diaria.rds")) %>%
  arrange(fecha)

## ---- 2.1 Variable dependiente: retorno acumulado a h dias, r_(t+h) ----
# r_(t+h) = 100 * (log S_(t+h) - log S_t) = suma de los retornos diarios
# individuales entre t+1 y t+h. Para h = 1 esto es simplemente retorno_tc
# adelantado un periodo.
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

## ---- 2.2 Predictores X_t (TODAS las variables de la base) ----
# Las variables ya acotadas -tasas, spreads y primas por plazo, que se
# interpretan directamente como un porcentaje o una diferencia de tasas- se
# dejan en NIVEL. Las variables de precio/indice sin cota superior (metales,
# petroleo, plata, indices bursatiles, monedas frente al dolar) se
# transforman a variacion logaritmica diaria para aproximarlas a
# estacionariedad, igual que se hace con retorno_tc para el propio tipo de
# cambio. El diferencial de tasas BCRP-Fed se construye explicitamente,
# tal como lo describe el Cuadro 2 del avance.
## NOTA: "tasa_peru" no se incluye como nivel individual porque, junto con
## "tasa_usa", forma una combinacion lineal EXACTA con "diff_tasas" (definida
## mas abajo como tasa_peru - tasa_usa). Incluir los 3 a la vez deja la
## matriz de diseno singular en la regresion cuantilica lineal (rq). BART y
## el QRF no tendrian ese problema, pero se excluye "tasa_peru" en los tres
## metodos para que todos usen exactamente el mismo set de predictores.
variables_nivel <- c(
  "vix", "baa10y",                     # riesgo sistemico global
  "tasa_usa", "nfci", "yield_to_worst", # condiciones financieras EE.UU.
  "tp_1y", "tp_5y", "tp_10y",           # prima por plazo (term premia ACM)
  "embig_peru"                          # incertidumbre politica Peru
)

variables_precio <- c(
  "rin",                                       # fundamentales macro Peru
  "oro", "cobre", "zinc", "plata", "wti",       # materias primas
  "spx", "mxef",                                # indices bursatiles
  "dxy", "jpy", "chf", "eur", "gbp",             # monedas (economias avanzadas)
  "clp", "cop", "mxn", "brl", "cny"              # monedas (America Latina + China)
)

base <- base %>%
  mutate(across(all_of(variables_precio), ~ c(NA, diff(log(.x))),
                .names = "d_{.col}")) %>%
  mutate(
    diff_tasas = tasa_peru - tasa_usa,
    # interv_bcrp es una variable de FLUJO con NA genuinos (ver "Base de
    # datos.R", seccion 5). Para poder usarla como predictor en BART -que
    # no admite NA en X- se asume que un dia sin reporte de intervencion
    # equivale a intervencion neta nula ese dia.
    interv_bcrp = ifelse(is.na(interv_bcrp), 0, interv_bcrp)
  )

predictores <- c(
  variables_nivel,
  "diff_tasas", "interv_bcrp",
  paste0("d_", variables_precio),
  "retorno_tc"   # retorno del dia t (momentum/rezago de la propia serie)
)

datos_modelo <- base %>%
  select(fecha, y_objetivo, all_of(predictores)) %>%
  na.omit()

cat("Observaciones disponibles para estimar:", nrow(datos_modelo), "\n")
cat("Rango de fechas:", as.character(min(datos_modelo$fecha)), "a",
    as.character(max(datos_modelo$fecha)), "\n")


# ============================================================================
# 3. SPLIT ENTRENAMIENTO / EVALUACION (cronologico, 80/20)
# ============================================================================
# En series de tiempo NO se separa aleatoriamente: el conjunto de evaluacion
# debe ser posterior en el tiempo al de entrenamiento.
n_total   <- nrow(datos_modelo)
corte     <- floor(0.8 * n_total)

train <- datos_modelo[1:corte, ]
test  <- datos_modelo[(corte + 1):n_total, ]

x_train <- as.matrix(train[, predictores])
y_train <- train$y_objetivo
x_test  <- as.matrix(test[, predictores])
y_test  <- test$y_objetivo

cat("Entrenamiento:", nrow(train), "obs. | Evaluacion:", nrow(test), "obs.\n")


# ============================================================================
# 4. BENCHMARK 1: CAMINATA ALEATORIA (Random Walk)
# ============================================================================
# Bajo random walk sin deriva, la mejor prediccion puntual es 0 (no hay
# cambio esperado en el retorno). Los cuantiles se estiman a partir de los
# cuantiles EMPIRICOS de los retornos historicos observados hasta t (regla
# de la raiz cuadrada del tiempo para escalar la dispersion a h periodos,
# seccion 5.1 del avance), y se mantienen fijos para todo el periodo de
# evaluacion (el propio random walk no usa covariables).
cuantiles_rw <- quantile(train$y_objetivo, probs = taus, na.rm = TRUE) *
  sqrt(h_pronostico)

pred_rw <- matrix(rep(cuantiles_rw, each = nrow(test)),
                   nrow = nrow(test), ncol = length(taus),
                   dimnames = list(NULL, paste0("tau_", taus)))


# ============================================================================
# 5. BENCHMARK 2: REGRESION CUANTILICA LINEAL (Koenker y Bassett, 1978)
# ============================================================================
# Un modelo lineal POR CADA cuantil tau (a diferencia de BART, que estima
# un unico modelo para toda la distribucion).
formula_rq <- as.formula(paste("y_objetivo ~", paste(predictores, collapse = " + ")))

modelos_rq <- map(taus, function(tau) rq(formula_rq, tau = tau, data = train))
names(modelos_rq) <- paste0("tau_", taus)

pred_rq <- sapply(modelos_rq, function(m) predict(m, newdata = test))
colnames(pred_rq) <- paste0("tau_", taus)


# ============================================================================
# 6. BENCHMARK 3: QUANTILE RANDOM FOREST (Meinshausen, 2006)
# ============================================================================
# Un unico bosque estima la distribucion condicional completa; los cuantiles
# se extraen de esa distribucion (a diferencia de rq, que necesita un ajuste
# por cuantil, pero igual que BART).
modelo_qrf <- quantregForest(x = x_train, y = y_train, ntree = 500)

pred_qrf <- predict(modelo_qrf, newdata = x_test, what = taus)
colnames(pred_qrf) <- paste0("tau_", taus)


# ============================================================================
# 7. MODELO PRINCIPAL: BART (Chipman, George & McCulloch, 2010)
# ============================================================================
# wbart() estima el modelo de suma de arboles Y = sum_j g(x; T_j, M_j) + eps,
# con eps ~ N(0, sigma^2). Guarda x.test para que devuelva de una vez las
# extracciones posteriores de la media condicional en el set de evaluacion.
bart_fit <- wbart(
  x.train = x_train, y.train = y_train, x.test = x_test,
  ntree = n_arboles, nskip = n_burn, ndpost = n_post
)

## ---- 7.1 Reconstruir la distribucion predictiva posterior COMPLETA ----
# yhat.test son extracciones de f(x) = E(y|x) para cada draw MCMC. Para
# obtener la distribucion PREDICTIVA (no solo la de la media condicional)
# hay que sumarle a cada extraccion un ruido ~ N(0, sigma_draw^2), donde
# sigma_draw es la extraccion de sigma correspondiente a esa misma iteracion
# (Clark et al., 2023, seccion de extraccion de cuantiles desde la
# distribucion predictiva posterior).
sigma_post <- tail(bart_fit$sigma, n_post)  # sigma post-burn-in, alineada con yhat.test

draws_predictivos_test <- bart_fit$yhat.test +
  matrix(rnorm(n_post * nrow(test), mean = 0, sd = sigma_post),
         nrow = n_post, ncol = nrow(test))

## ---- 7.2 Cuantiles condicionales Q_tau(x) para cada observacion de test ----
pred_bart <- apply(draws_predictivos_test, 2, quantile, probs = taus)
pred_bart <- t(pred_bart)
colnames(pred_bart) <- paste0("tau_", taus)


# ============================================================================
# 8. EVALUACION FUERA DE MUESTRA: QUANTILE SCORE Y QW-CRPS
# ============================================================================

## ---- 8.1 Quantile Score (funcion de perdida pinball, ec. en seccion 5.2) ----
pinball <- function(y_real, y_pred, tau) {
  u <- y_real - y_pred
  mean(u * (tau - as.numeric(u < 0)))
}

qs_metodo <- function(pred_matrix, y_real, taus) {
  sapply(seq_along(taus), function(i) pinball(y_real, pred_matrix[, i], taus[i]))
}

qs_rw   <- qs_metodo(pred_rw,   y_test, taus)
qs_rq   <- qs_metodo(pred_rq,   y_test, taus)
qs_qrf  <- qs_metodo(pred_qrf,  y_test, taus)
qs_bart <- qs_metodo(pred_bart, y_test, taus)

tabla_quantile_score <- data.frame(
  tau               = taus,
  random_walk       = qs_rw,
  quantile_reg      = qs_rq,
  quantile_rf       = qs_qrf,
  bart              = qs_bart
)

## ---- 8.2 Quantile-weighted CRPS: promedio del QS sobre todos los tau ----
# Aproxima el CRPS ponderando la perdida pinball en una grilla fina de
# cuantiles (aqui, los 6 tau de interes del avance como aproximacion).
qw_crps <- colMeans(tabla_quantile_score[, -1])

cat("\n--- Quantile Score por metodo y cuantil ---\n")
print(round(tabla_quantile_score, 5))
cat("\n--- Quantile-weighted CRPS (promedio sobre todos los tau) ---\n")
print(round(qw_crps, 5))

write.csv(tabla_quantile_score,
          file.path(ruta_out, "quantile_score_comparacion.csv"), row.names = FALSE)


# ============================================================================
# 9. DETERMINANTES PRINCIPALES: IMPORTANCIA GLOBAL (BART)
# ============================================================================
# Frecuencia con la que cada predictor es usado como regla de division en
# los arboles, promediada sobre todas las extracciones MCMC posteriores al
# burn-in (Chipman et al., 2010, seccion de "MCMC backfitting bayesiano").
importancia_global <- colMeans(bart_fit$varcount)
importancia_global <- sort(importancia_global, decreasing = TRUE)

tabla_importancia <- data.frame(
  variable   = names(importancia_global),
  importancia = as.numeric(importancia_global)
)

cat("\n--- Importancia global de variables (BART, frecuencia de uso en arboles) ---\n")
print(tabla_importancia)

write.csv(tabla_importancia,
          file.path(ruta_out, "bart_importancia_global.csv"), row.names = FALSE)

grafico_importancia <- ggplot(tabla_importancia,
                               aes(x = reorder(variable, importancia), y = importancia)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Importancia global de predictores (BART)",
       x = NULL, y = "Frecuencia promedio de uso en los arboles") +
  theme_minimal()

ggsave(file.path(ruta_figs, "bart_importancia_global.png"),
       grafico_importancia, width = 7, height = 5)


# ============================================================================
# 10. DETERMINANTES POR CUANTIL: DEPENDENCIA PARCIAL BAYESIANA
# ============================================================================
# Para cada predictor, se construye una grilla con su rango empirico
# (percentiles 5, 25, 50, 75, 95 de la muestra de entrenamiento) manteniendo
# TODOS los demas predictores fijos en su mediana. Se predice con el MISMO
# modelo BART ya estimado (sin volver a entrenar) y se reconstruye, igual
# que en la seccion 7, la distribucion predictiva posterior en cada punto
# de la grilla. Esto revela como se desplaza cada Q_tau al mover un
# predictor -el "mapa completo" de la seccion 5.2 del avance, pero
# obtenido de un unico modelo no lineal en vez de un modelo lineal por tau.
medianas_train <- apply(x_train, 2, median)

calcular_dependencia_parcial <- function(variable, n_grilla = 5) {
  grilla <- quantile(x_train[, variable], probs = c(0.05, 0.25, 0.50, 0.75, 0.95))

  x_grilla <- matrix(rep(medianas_train, each = length(grilla)),
                      nrow = length(grilla), dimnames = list(NULL, predictores))
  x_grilla[, variable] <- grilla

  # predict.wbart reutiliza los arboles ya ajustados (no reestima el modelo)
  pred <- predict(bart_fit, newdata = x_grilla)
  sigma_pred <- tail(bart_fit$sigma, nrow(pred))

  draws <- pred + matrix(rnorm(length(pred), mean = 0, sd = sigma_pred),
                          nrow = nrow(pred), ncol = ncol(pred))

  cuantiles_grilla <- t(apply(draws, 2, quantile, probs = taus))
  colnames(cuantiles_grilla) <- paste0("tau_", taus)

  data.frame(variable = variable, valor_predictor = grilla, cuantiles_grilla,
             row.names = NULL)
}

dependencia_parcial <- map_dfr(predictores, calcular_dependencia_parcial)

cat("\n--- Dependencia parcial por cuantil (primeras filas) ---\n")
print(head(dependencia_parcial, 10))

write.csv(dependencia_parcial,
          file.path(ruta_out, "bart_dependencia_parcial_por_cuantil.csv"),
          row.names = FALSE)

## ---- 10.1 Grafico: efecto marginal por cuantil, para cada predictor ----
dependencia_larga <- dependencia_parcial %>%
  pivot_longer(cols = starts_with("tau_"), names_to = "cuantil", values_to = "valor_y")

grafico_dependencia <- ggplot(dependencia_larga,
                               aes(x = valor_predictor, y = valor_y, color = cuantil)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variable, scales = "free_x") +
  labs(title = "Efecto marginal de cada predictor por cuantil (BART)",
       x = "Valor del predictor", y = "Cuantil condicional del retorno (r_t+h)",
       color = "Cuantil") +
  theme_minimal()

ggsave(file.path(ruta_figs, "bart_dependencia_parcial_por_cuantil.png"),
       grafico_dependencia, width = 18, height = 16, limitsize = FALSE)


# ============================================================================
# 11. RESUMEN FINAL
# ============================================================================
cat("\n============================================================\n")
cat("Resultados guardados en:\n")
cat(" -", file.path(ruta_out, "quantile_score_comparacion.csv"), "\n")
cat(" -", file.path(ruta_out, "bart_importancia_global.csv"), "\n")
cat(" -", file.path(ruta_out, "bart_dependencia_parcial_por_cuantil.csv"), "\n")
cat(" -", file.path(ruta_figs, "bart_importancia_global.png"), "\n")
cat(" -", file.path(ruta_figs, "bart_dependencia_parcial_por_cuantil.png"), "\n")
cat("============================================================\n")
