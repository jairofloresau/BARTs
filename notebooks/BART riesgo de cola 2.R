# ============================================================================
# TESIS: Pronostico e Identificacion de Determinantes de los Retornos
#        Extremos del Tipo de Cambio en Peru: Un Enfoque BART
#
# SCRIPT 02: Riesgo de cola con BART - distribucion completa y determinantes
#            [VERSION CORREGIDA]
# ============================================================================
#
# QUE HACE ESTE SCRIPT, EN ORDEN:
#
#   1. Carga "base_final_diaria" (construida en "Base de datos.R") y arma la
#      variable dependiente r_(t+h) y el set completo de predictores X_t.
#   2. Parte la muestra en TRES bloques cronologicos: entrenamiento (60%),
#      validacion (20%, para tunear hiperparametros) y evaluacion final
#      (20%, que no se toca hasta tener el modelo elegido).
#   3. TUNEA los hiperparametros de BART (k, ntree, sparse) por busqueda en
#      grilla sobre el bloque de validacion, seleccionando por qwCRPS.
#   4. Estima BART HETEROSCEDASTICO en dos etapas: un BART para la media
#      f(x) y otro para la log-varianza log sigma^2(x). De ahi reconstruye la
#      distribucion predictiva posterior COMPLETA (Clark, Huber, Koop,
#      Marcellino & Pfarrhofer, 2023).
#   5. Evalua fuera de muestra con VENTANA EXPANSIVA (recursiva) los 4
#      metodos sobre EXACTAMENTE los mismos origenes:
#        - Caminata aleatoria (Random Walk)
#        - Regresion cuantilica lineal (Koenker y Bassett, 1978)
#        - Quantile Random Forest (Meinshausen, 2006)
#        - BART heteroscedastico (modelo principal)
#      Metricas: Quantile Score (pinball) por cuantil y qwCRPS PONDERADO POR
#      COLA (Gneiting & Ranjan, 2011).
#   6. Identifica los determinantes con TRES vistas complementarias:
#        a) determinantes de la MEDIA      -> efecto de NIVEL (todos los tau)
#        b) determinantes de la VARIANZA   -> efecto de DISPERSION (las COLAS)
#        c) importancia POR CUANTIL        -> tabla variable x tau
#      mas la dependencia parcial bayesiana promediada sobre la muestra.
#
# ----------------------------------------------------------------------------
# NOTA METODOLOGICA CENTRAL (por que el modelo es HETEROSCEDASTICO)
# ----------------------------------------------------------------------------
# Un BART homoscedastico implica  y | x ~ N( f(x), sigma^2 ), y por lo tanto
#
#       Q_tau(x) = f(x) + sigma * z_tau ,     con z_tau = qnorm(tau)
#
# El unico termino que depende de x -f(x)- es EL MISMO para todos los tau;
# entre cuantiles solo cambia el desplazamiento vertical sigma*z_tau. Es decir:
# con sigma constante, los determinantes de las colas y del centro son
# identicos POR CONSTRUCCION DEL MODELO, no porque lo digan los datos. Eso
# vaciaria de contenido la pregunta central de la tesis ("¿los determinantes
# de los retornos EXTREMOS son los mismos que los del centro?").
#
# Por eso se deja que la VARIANZA dependa de x (en la linea del hBART de
# Clark et al., 2023, y de Pratola et al., 2020):
#
#       y | x ~ N( f(x), sigma^2(x) )
#       Q_tau(x)      = f(x) + sigma(x) * z_tau
#       dQ_tau / dx   = df/dx  +  z_tau * dsigma/dx
#
# El segundo termino escala con z_tau: es ~0 en el centro (z_0.5 = 0) y crece
# en las colas. Entonces:
#   - las variables que mueven la MEDIA f(x)        -> afectan TODOS los cuantiles
#   - las variables que mueven la VARIANZA sigma(x) -> pesan sobre todo en las COLAS
# Recien asi los determinantes PUEDEN diferir por cuantil, y que difieran (o no)
# pasa a ser un HALLAZGO EMPIRICO en vez de un artefacto del modelo.
#
# Ademas, se trata de retornos cambiarios DIARIOS, donde el clustering de
# volatilidad es fuerte y gobierna precisamente las colas. El resultado de
# Clark et al. de que "la heterocedasticidad importa menos cuando BART modela
# bien la media" se obtiene con datos macro trimestrales/mensuales; no es
# trasladable sin mas a frecuencia diaria.
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

# --- Horizonte de pronostico (en dias) y cuantiles de interes ---
h_pronostico <- 1
taus <- c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99)
# tau = 0.50 es la REFERENCIA CENTRAL: sin ella no se puede contrastar
# "colas vs centro", que es el objetivo de identificacion de la tesis.

# --- MCMC de BART: cadena COMPLETA (para el modelo final) ---
n_burn <- 1000
n_post <- 2000

# --- MCMC de BART: cadena REDUCIDA (solo para el tuneo) ---
# Para COMPARAR configuraciones entre si no hace falta la cadena completa: el
# ordenamiento relativo se estabiliza mucho antes que la posterior. La
# configuracion ganadora se reestima despues con la cadena completa.
n_burn_tune <- 250
n_post_tune <- 500

# --- Control del costo computacional de la ventana expansiva ---
# Reestimar BART en cada dia del periodo de evaluacion es inviable. Se
# reestima cada 'paso_reest' observaciones y, entre reestimaciones, se
# REUTILIZA el ultimo modelo aplicandolo a la nueva x_t (los coeficientes no
# cambian, pero la prediccion si, porque cambia el x del dia).
#   paso_reest = 20  -> ~1 mes bursatil. Mas fiel, mas lento.
#   paso_reest = 60  -> ~1 trimestre. Corrida exploratoria.
paso_reest <- 60

# Numero de filas de la muestra usadas para promediar la dependencia parcial
# (ver Seccion 9). Subirlo mejora la aproximacion; bajarlo acelera.
n_muestra_pdp <- 100


# ============================================================================
# 2. CARGAR BASE Y CONSTRUIR VARIABLE DEPENDIENTE + PREDICTORES
# ============================================================================

base <- readRDS(file.path(ruta_datos, "base_final_diaria.rds")) %>%
  arrange(fecha)

# Pasar los retornos de decimales a porcentaje
base <- base %>% mutate(retorno_tc = 100 * retorno_tc)

## ---- 2.1 Variable dependiente: retorno acumulado a h dias, r_(t+h) ----
# r_(t+h) = suma de los retornos diarios entre t+1 y t+h. Para h = 1 esto es
# simplemente retorno_tc adelantado un periodo. Notar que la informacion usada
# para predecir r_(t+h) es la disponible EN t (no hay fuga de informacion).
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

## ---- 2.2 Predictores X_t ----
# Las variables ya acotadas -tasas, spreads y primas por plazo, que se
# interpretan directamente como un porcentaje o una diferencia de tasas- se
# dejan en NIVEL. Las variables de precio/indice sin cota superior (metales,
# petroleo, indices bursatiles, monedas frente al dolar) se transforman a
# variacion logaritmica diaria para aproximarlas a estacionariedad, igual que
# se hace con retorno_tc para el propio tipo de cambio.
#
# NOTA: "tasa_peru" no se incluye como nivel individual porque, junto con
# "tasa_usa", forma una combinacion lineal EXACTA con "diff_tasas". Incluir
# las 3 a la vez deja la matriz de disenio singular en la regresion cuantilica
# lineal (rq). BART y el QRF no tendrian ese problema, pero se excluye en los
# tres metodos para que todos usen exactamente el mismo set de predictores.
variables_nivel <- c(
  "vix", "baa10y",                      # riesgo sistemico global
  "tasa_usa", "nfci", "yield_to_worst", # condiciones financieras EE.UU.
  "tp_1y", "tp_5y", "tp_10y",           # prima por plazo (term premia ACM)
  "embig_peru"                          # riesgo pais / incertidumbre Peru
)

variables_precio <- c(
  "rin",                                  # fundamentales macro Peru
  "oro", "cobre", "zinc", "plata", "wti", # materias primas
  "spx", "mxef",                          # indices bursatiles
  "dxy", "jpy", "chf", "eur", "gbp",      # monedas (economias avanzadas)
  "clp", "cop", "mxn", "brl", "cny"       # monedas (America Latina + China)
)

base <- base %>%
  mutate(across(all_of(variables_precio), ~ 100 * c(NA, diff(log(.x))),
                .names = "d_{.col}")) %>%
  mutate(
    diff_tasas = tasa_peru - tasa_usa,
    # interv_bcrp es una variable de FLUJO con NA genuinos. Un dia sin reporte
    # de intervencion se asume como intervencion neta nula ese dia.
    interv_bcrp = ifelse(is.na(interv_bcrp), 0, interv_bcrp)
  )

predictores <- c(
  variables_nivel,
  "diff_tasas", "interv_bcrp",
  paste0("d_", variables_precio),
  "retorno_tc"   # retorno del dia t (momentum / rezago de la propia serie)
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
# ============================================================================
# En series de tiempo NO se separa aleatoriamente. Ademas se necesitan TRES
# bloques y no dos: el bloque de VALIDACION sirve para tunear hiperparametros,
# y el de TEST queda intacto para la evaluacion final. Si se tunea mirando el
# test, el desempenio fuera de muestra queda contaminado (look-ahead).

n_total   <- nrow(datos_modelo)
corte_tr  <- floor(0.60 * n_total)   # fin del entrenamiento
corte_val <- floor(0.80 * n_total)   # fin de la validacion = inicio del test

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


# ============================================================================
# 4. METRICAS DE EVALUACION
# ============================================================================

## ---- 4.1 Quantile Score (funcion de perdida pinball) ----
pinball <- function(y_real, y_pred, tau) {
  u <- y_real - y_pred
  mean(u * (tau - as.numeric(u < 0)))
}

qs_metodo <- function(pred_matrix, y_real, taus) {
  sapply(seq_along(taus), function(i) pinball(y_real, pred_matrix[, i], taus[i]))
}

## ---- 4.2 qwCRPS: CRPS PONDERADO POR COLA (Gneiting & Ranjan, 2011) ----
# El CRPS admite una representacion como integral de la perdida pinball sobre
# tau. Ponderando esa integral se enfoca la evaluacion en una region concreta:
#
#   qwCRPS = sum_tau [ w(tau) * QS_tau ] / sum_tau w(tau)
#
#   w(tau) = 1              -> CRPS "plano" (toda la distribucion)
#   w(tau) = (1 - tau)^2    -> COLA IZQUIERDA (retornos muy negativos)
#   w(tau) = tau^2          -> COLA DERECHA   (retornos muy positivos =
#                              depreciaciones bruscas del sol, dado que el TC
#                              esta en S/. por US$)
#   w(tau) = (2*tau - 1)^2  -> AMBAS COLAS simultaneamente
#
# IMPORTANTE: se integra sobre una grilla FINA de tau (no solo los 7 tau que
# se reportan), que es lo correcto para aproximar la integral del CRPS. La
# version anterior del script promediaba el QS sobre los 7 tau sin pesos: eso
# NO es un CRPS ponderado, es un promedio simple.
taus_crps <- seq(0.01, 0.99, by = 0.01)

qwcrps <- function(pred_q_fina, y_real, taus_grid = taus_crps,
                   peso = c("ambas_colas", "izquierda", "derecha", "plano")) {
  peso <- match.arg(peso)
  
  w <- switch(peso,
              plano       = rep(1, length(taus_grid)),
              izquierda   = (1 - taus_grid)^2,
              derecha     = taus_grid^2,
              ambas_colas = (2 * taus_grid - 1)^2
  )
  
  qs <- vapply(seq_along(taus_grid),
               function(i) pinball(y_real, pred_q_fina[, i], taus_grid[i]),
               numeric(1))
  
  sum(w * qs) / sum(w)
}

# Metrica de SELECCION del tuneo y de comparacion principal entre metodos.
# "ambas_colas" es la coherente con el objetivo de la tesis (retornos EXTREMOS
# en ambos sentidos). Si el interes fuera solo la depreciacion brusca del sol,
# se cambia a "derecha".
peso_seleccion <- "ambas_colas"


# ============================================================================
# 5. BART HETEROSCEDASTICO EN DOS ETAPAS
# ============================================================================
# Etapa 1: BART sobre la media          -> f(x)  y residuos  e = y - f_hat(x)
# Etapa 2: BART sobre log(e^2 + c)      -> log sigma^2(x)  ->  sigma(x)
# Reconstruccion: y*_k(x) = f_k(x) + sigma(x) * N(0,1)
#
# ADVERTENCIAS METODOLOGICAS (declararlas en la tesis):
#  - Es la version IMPLEMENTABLE con el paquete `BART`: estima media y varianza
#    en DOS ETAPAS (estilo GLS factible), no con un muestreador conjunto. La
#    version conjunta (paquete `rbart`, Pratola et al. 2020) o con volatilidad
#    estocastica quedan como chequeos de ROBUSTEZ.
#  - sigma(x) entra como PLUG-IN (media posterior de la funcion de log-varianza).
#    La incertidumbre de la MEDIA si se propaga via los draws MCMC; la de la
#    varianza no. Esto subestima levemente la incertidumbre total.
#  - Modelar log(e^2) y volver con exp() introduce un sesgo pequenio de escala
#    (E[log X] != log E[X]). Para IDENTIFICAR que variables mueven la dispersion
#    es suficiente; para calibracion fina de la cola, es una limitacion.
#  - Los hiperparametros tuneados se aplican al BART de la MEDIA. El de la
#    VARIANZA queda en defaults, para no multiplicar por 24 el costo del tuneo.

# --- Helper: evita el caso degenerado de x.test/newdata de UNA sola fila ---
# El paquete BART colapsa dimensiones cuando x.test tiene 1 fila (yhat.test /
# predict() dejan de devolver una matriz ndpost x 1 consistente). Se duplica
# la fila, se predice sobre 2, y despues se recorta al numero de filas
# realmente pedido.
blindar_x <- function(x) {
  x <- as.matrix(x)
  if (nrow(x) == 1) rbind(x, x) else x
}

ajustar_bart_het <- function(x_train, y_train, x_pred,
                             k = 2, ntree = 200, sparse = FALSE,
                             nskip = n_burn, ndpost = n_post) {
  
  prt      <- ndpost + 1L
  n_real   <- nrow(as.matrix(x_pred))
  x_pred_b <- blindar_x(x_pred)
  
  ## --- Etapa 1: media condicional f(x)  [HIPERPARAMETROS TUNEADOS] ---
  invisible(capture.output({
    fit_media <- wbart(
      x.train = x_train, y.train = y_train, x.test = x_pred_b,
      k = k, ntree = ntree, sparse = sparse,
      nskip = nskip, ndpost = ndpost, printevery = prt
    )
  }))
  
  ## --- Etapa 2: log-varianza condicional ---
  # ntree = 50: la funcion de log-varianza es SUAVE, no necesita 200 arboles.
  # Con 200 esta etapa costaba mas que el modelo principal.
  f_hat <- colMeans(matrix(fit_media$yhat.train, nrow = ndpost))
  resid <- y_train - f_hat
  c_est <- 1e-6 * stats::var(resid)
  
  invisible(capture.output({
    fit_var <- wbart(
      x.train = x_train, y.train = log(resid^2 + c_est), x.test = x_pred_b,
      ntree = 50,
      nskip = nskip, ndpost = ndpost, printevery = prt
    )
  }))
  
  ym <- matrix(fit_media$yhat.test, nrow = ndpost)[, 1:n_real, drop = FALSE]
  yv <- matrix(fit_var$yhat.test,   nrow = ndpost)[, 1:n_real, drop = FALSE]
  
  sigma_pred <- sqrt(exp(colMeans(yv)))
  
  ruido <- matrix(rnorm(ndpost * n_real), nrow = ndpost, ncol = n_real)
  ruido <- sweep(ruido, 2, sigma_pred, `*`)
  
  list(
    draws          = ym + ruido,
    sigma_pred     = sigma_pred,
    varcount_media = colMeans(fit_media$varcount),
    varcount_var   = colMeans(fit_var$varcount),
    varprob_media  = if (!is.null(fit_media$varprob)) colMeans(fit_media$varprob) else NULL,
    fit_media      = fit_media,
    fit_var        = fit_var,
    ndpost         = ndpost
  )
}

# Reutiliza un modelo YA ajustado para reconstruir la distribucion predictiva
# en puntos nuevos, SIN volver a estimar. Se usa (a) entre reestimaciones de la
# ventana expansiva y (b) en la dependencia parcial.
reconstruir_draws <- function(modelo, x_nuevo) {
  nd     <- modelo$ndpost
  n_real <- nrow(as.matrix(x_nuevo))
  x_b    <- blindar_x(x_nuevo)
  
  invisible(capture.output({
    pm_raw <- predict(modelo$fit_media, newdata = x_b)
    pv_raw <- predict(modelo$fit_var,   newdata = x_b)
  }))
  
  pm <- matrix(pm_raw, nrow = nd)[, 1:n_real, drop = FALSE]
  pv <- matrix(pv_raw, nrow = nd)[, 1:n_real, drop = FALSE]
  
  sig   <- sqrt(exp(colMeans(pv)))
  ruido <- sweep(matrix(rnorm(nd * n_real), nrow = nd, ncol = n_real), 2, sig, `*`)
  pm + ruido
}

# ============================================================================
# 6. TUNEO DE HIPERPARAMETROS: GRILLA 4 x 3 x 2 = 24 CONFIGURACIONES
# ============================================================================
#   k      in {1, 2, 3, 5}   -> shrinkage de los nodos terminales. LA PERILLA
#                               MAS INFLUYENTE. k alto encoge mas (ajuste mas
#                               suave, menos sobreajuste). Con retornos diarios
#                               -senial debil, mucho ruido- es esperable que
#                               gane un k alto.
#   ntree  in {50, 100, 200} -> numero de arboles (la "m" de Chipman et al.).
#                               Rendimientos decrecientes pasados los ~200.
#   sparse in {FALSE, TRUE}  -> prior de Dirichlet disperso (DART, Linero 2018).
#                               Apaga predictores irrelevantes. Con ~30
#                               predictores (metales, monedas, spreads) tiene
#                               doble valor: mejora el pronostico Y limpia la
#                               identificacion de determinantes.
#
# alpha (base=0.95), beta (power=2), nu (sigdf=3) y q (sigquant=0.90) se dejan
# en los defaults de Chipman et al. (2010): aportan poco frente a k y ntree, y
# tunearlos multiplicaria la grilla sin ganancia clara.

grilla <- expand.grid(
  k      = c(1, 2, 3, 5),
  ntree  = c(50, 100, 200),
  sparse = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

cat("=== TUNEO: evaluando", nrow(grilla), "configuraciones sobre la validacion ===\n\n")

resultados_tuneo <- vector("list", nrow(grilla))

for (i in seq_len(nrow(grilla))) {
  
  cfg <- grilla[i, ]
  t0  <- Sys.time()
  
  ajuste <- ajustar_bart_het(
    x_train = x_train, y_train = y_train, x_pred = x_val,
    k = cfg$k, ntree = cfg$ntree, sparse = cfg$sparse,
    nskip = n_burn_tune, ndpost = n_post_tune
  )
  
  # cuantiles en la grilla FINA (para el qwCRPS)
  Q_fina <- t(apply(ajuste$draws, 2, quantile, probs = taus_crps))
  
  m_sel   <- qwcrps(Q_fina, y_val, peso = peso_seleccion)
  m_izq   <- qwcrps(Q_fina, y_val, peso = "izquierda")
  m_der   <- qwcrps(Q_fina, y_val, peso = "derecha")
  m_plano <- qwcrps(Q_fina, y_val, peso = "plano")
  
  # RMSE solo como REFERENCIA. NO se usa para seleccionar: optimizaria el
  # centro de la distribucion, y esta tesis es sobre las COLAS.
  rmse <- sqrt(mean((y_val - colMeans(ajuste$fit_media$yhat.test))^2))
  
  resultados_tuneo[[i]] <- data.frame(
    k = cfg$k, ntree = cfg$ntree, sparse = cfg$sparse,
    qwcrps_seleccion = m_sel,
    qwcrps_izquierda = m_izq,
    qwcrps_derecha   = m_der,
    crps_plano       = m_plano,
    rmse_referencia  = rmse,
    minutos          = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  )
  
  cat(sprintf("[%2d/%2d] k=%g  ntree=%3d  sparse=%-5s  ->  qwCRPS=%.5f  (%.1f min)\n",
              i, nrow(grilla), cfg$k, cfg$ntree, cfg$sparse, m_sel,
              resultados_tuneo[[i]]$minutos))
  
  rm(ajuste); gc(verbose = FALSE)   # los objetos wbart pesan; liberar memoria
}

tabla_tuneo <- do.call(rbind, resultados_tuneo)
tabla_tuneo <- tabla_tuneo[order(tabla_tuneo$qwcrps_seleccion), ]

cat("\n--- Resultados del tuneo (ordenados por qwCRPS; menor = mejor) ---\n")
print(round(tabla_tuneo, 5), row.names = FALSE)

write.csv(tabla_tuneo, file.path(ruta_out, "bart_tuneo_grilla.csv"), row.names = FALSE)

## ---- 6.1 Configuracion ganadora ----
mejor      <- tabla_tuneo[1, ]
k_opt      <- mejor$k
ntree_opt  <- mejor$ntree
sparse_opt <- mejor$sparse

cat("\n============================================================\n")
cat("CONFIGURACION SELECCIONADA (menor qwCRPS en validacion):\n")
cat("  k      =", k_opt, "\n")
cat("  ntree  =", ntree_opt, "\n")
cat("  sparse =", sparse_opt, "\n")
cat("  qwCRPS =", round(mejor$qwcrps_seleccion, 5), "\n")
cat("============================================================\n")

## ---- 6.2 Comparacion contra el DEFAULT: la referencia honesta ----
# BART es conocido por funcionar bien con sus defaults; Chipman et al. los
# calibraron justamente para eso. Si la ganancia del tuneo frente a
# (k=2, ntree=200, sparse=FALSE) es marginal, ESO ES UN RESULTADO REPORTABLE,
# no un fracaso: valida la practica estandar de la literatura y muestra que el
# desempenio no se inflo con busqueda excesiva de hiperparametros.
fila_default <- subset(tabla_tuneo, k == 2 & ntree == 200 & sparse == FALSE)
if (nrow(fila_default) == 1) {
  ganancia <- 100 * (fila_default$qwcrps_seleccion - mejor$qwcrps_seleccion) /
    fila_default$qwcrps_seleccion
  cat("\nqwCRPS del default (k=2, ntree=200, sparse=FALSE):",
      round(fila_default$qwcrps_seleccion, 5), "\n")
  cat("Mejora del tuneo sobre el default:", round(ganancia, 2), "%\n\n")
}


sd(datos_modelo$y_objetivo)   # debería dar ~0.3 a 0.5

## AGREGADO

# ---- Sustituto de la §6 (nos saltamos el tuneo) ----
k_opt <- 2; ntree_opt <- 200; sparse_opt <- FALSE

# ---- Objeto que normalmente se crea dentro de la §7 ----
formula_rq <- as.formula(paste("y_objetivo ~", paste(predictores, collapse = " + ")))

# Ventana tipica (la del primer origen del test)
df_h <- datos_modelo[1:5508, ]
x_h  <- as.matrix(df_h[, predictores])
y_h  <- df_h$y_objetivo
df_1 <- datos_modelo[5509, , drop = FALSE]
x_1  <- as.matrix(df_1[, predictores])

cat("--- AJUSTE (una vez cada 60 dias) ---\n")
system.time({ m_rq  <- rq(formula_rq, tau = taus_crps, data = df_h) })
system.time({ m_qrf <- quantregForest(x = x_h, y = y_h, ntree = 500) })
system.time({ m_bart <- ajustar_bart_het(x_h, y_h, x_1, k = k_opt,
                                         ntree = ntree_opt, sparse = sparse_opt,
                                         nskip = n_burn, ndpost = n_post) })

cat("--- PREDICT (20 dias, para extrapolar a los 59 intermedios) ---\n")
system.time({ for(i in 1:20) predict(m_rq, newdata = df_1) })
system.time({ for(i in 1:20) predict(m_qrf, newdata = x_1, what = taus_crps) })
system.time({ for(i in 1:20) reconstruir_draws(m_bart, x_1) })

# ============================================================================
# 7. EVALUACION FUERA DE MUESTRA CON VENTANA EXPANSIVA (RECURSIVA POR BLOQUES)
# ============================================================================
# Para cada origen t del periodo de test se pronostica r_(t+h) con la
# informacion disponible hasta t-1. La ventana se EXPANDE al avanzar t.
#
# Por que importa: es el estandar de la literatura de pronostico de colas
# (Clark et al., 2023) y, sobre todo, acumula muchas mas observaciones de cola
# para evaluar tau = 0.01 y 0.99. Con un corte unico, el Quantile Score al 1%
# podria estar determinado por un punado de dias.
#
# LOS 4 METODOS SE EVALUAN SOBRE EXACTAMENTE LOS MISMOS ORIGENES, condicion
# necesaria para que la comparacion sea justa.
#
# OPTIMIZACION POR BLOQUES: entre dos reestimaciones consecutivas el modelo
# (rq, QRF, BART) NO cambia -solo cambia el x_t que se le pasa-. En vez de
# iterar dia por dia y, entre reestimaciones, volver a llamar a predict()
# fila por fila (via reconstruir_draws()), se agrupan los "paso_reest"
# origenes de cada bloque y se predicen TODOS DE UNA SOLA VEZ pasando la
# matriz completa del bloque como x.test/newdata. Es matematicamente
# identico -misma ventana de entrenamiento, mismo x_t por fila-, solo evita
# reabrir el objeto ajustado una vez por fila. Esto reduce drasticamente el
# numero de llamadas a predict() de BART (de ~1 por dia a 1 por bloque).

origenes <- idx_test

# Contenedores: una matriz por metodo, con los cuantiles de la grilla FINA
# (de ahi se extraen tanto los 7 tau reportados como el qwCRPS).
n_or <- length(origenes)
n_tf <- length(taus_crps)

pred_fina <- list(
  random_walk  = matrix(NA_real_, n_or, n_tf),
  quantile_reg = matrix(NA_real_, n_or, n_tf),
  quantile_rf  = matrix(NA_real_, n_or, n_tf),
  bart         = matrix(NA_real_, n_or, n_tf)
)
y_real <- datos_modelo$y_objetivo[origenes]

formula_rq <- as.formula(paste("y_objetivo ~", paste(predictores, collapse = " + ")))

# Bloques de "paso_reest" origenes consecutivos. Como "origenes" es una
# secuencia de indices de fila CONSECUTIVOS (idx_test), agrupar por posicion
# equivale a agrupar por "t - ultimo_ajuste >= paso_reest" de la version dia
# a dia: en cada bloque se reestima UNA vez, con la ventana hasta justo antes
# del primer origen del bloque, y se predicen de un tiron los origenes de ese
# bloque (los ultimos, "adentro" del bloque, se predicen con un modelo
# estimado con menos informacion que si se hubiera reestimado dia a dia; ese
# es exactamente el trade-off que ya asumia la version original con
# reconstruir_draws()).
n_bloques <- ceiling(n_or / paso_reest)

cat("=== EVALUACION RECURSIVA POR BLOQUES:", n_or, "origenes,",
    n_bloques, "bloques de hasta", paso_reest, "obs., reestimando 1 vez por bloque ===\n\n")

t_inicio_eval <- Sys.time()

for (b in seq_len(n_bloques)) {

  idx_bloque <- ((b - 1) * paso_reest + 1):min(b * paso_reest, n_or)
  t_bloque   <- origenes[idx_bloque]
  t_ini      <- min(t_bloque)

  ## ---- Reestimacion, una vez por bloque, con la ventana hasta t_ini - 1 ----
  idx_hasta <- 1:(t_ini - 1)
  df_hasta  <- datos_modelo[idx_hasta, ]
  x_hasta   <- as.matrix(df_hasta[, predictores])
  y_hasta   <- df_hasta$y_objetivo

  df_b <- datos_modelo[t_bloque, , drop = FALSE]
  x_b  <- as.matrix(df_b[, predictores])

  # -- BENCHMARK 1: Random Walk --
  # Bajo random walk sin deriva la mejor prediccion puntual es 0; la
  # distribucion se aproxima con los cuantiles EMPIRICOS de los retornos
  # observados hasta t_ini-1. El RW no usa covariables, asi que sus cuantiles
  # son los mismos para todo x del bloque: solo se actualizan al reestimar.
  #
  # CORRECCION: la version anterior multiplicaba por sqrt(h_pronostico). Eso
  # era un DOBLE ESCALAMIENTO: y_objetivo YA es el retorno acumulado a h
  # dias, asi que sus cuantiles empiricos ya incorporan la dispersion de h
  # periodos. La regla de la raiz del tiempo sirve para escalar cuantiles de
  # 1 dia HACIA h dias, no aqui. Inocuo con h=1 (sqrt(1)=1), pero descalibra
  # el benchmark en cuanto h > 1.
  cuantiles_rw <- quantile(y_hasta, probs = taus_crps, na.rm = TRUE)

  # -- BENCHMARK 2: Regresion cuantilica lineal (Koenker y Bassett, 1978) --
  # Un modelo lineal POR CADA cuantil. Se ajusta sobre la grilla fina para
  # que el qwCRPS sea comparable con el de los demas metodos. Se predicen de
  # una vez todas las filas del bloque (predict() de rq acepta newdata con
  # varias filas y devuelve una matriz n_bloque x n_tau).
  # NOTA: al estimar cada tau por separado, rq puede producir CRUCE DE
  # CUANTILES (Q_0.05 > Q_0.10). Se reordena cada fila con sort() -la
  # solucion de reordenamiento de Chernozhukov et al. (2010)-. Que BART no
  # necesite esta correccion, por venir todos los cuantiles de una unica
  # distribucion posterior, es justamente una de sus ventajas.
  modelo_rq  <- rq(formula_rq, tau = taus_crps, data = df_hasta)
  pred_rq_b  <- predict(modelo_rq, newdata = df_b)
  pred_rq_b  <- t(apply(pred_rq_b, 1, sort))

  # -- BENCHMARK 3: Quantile Random Forest (Meinshausen, 2006) --
  # ntree baja de 500 a 200: con la reestimacion por bloques ya no es el
  # cuello de botella, pero 500 arboles seguia costando ~126s por
  # reestimacion sin mejorar materialmente los cuantiles.
  modelo_qrf  <- quantregForest(x = x_hasta, y = y_hasta, ntree = 200)
  pred_qrf_b  <- predict(modelo_qrf, newdata = x_b, what = taus_crps)

  # -- MODELO PRINCIPAL: BART heteroscedastico con hiperparametros TUNEADOS --
  # Se predicen las hasta "paso_reest" filas del bloque en una sola llamada a
  # wbart (x_pred = x_b), en vez de una llamada a predict() por fila via
  # reconstruir_draws(). Matematicamente identico: misma ventana de
  # entrenamiento para todas las filas del bloque.
  ajuste_bart <- ajustar_bart_het(
    x_train = x_hasta, y_train = y_hasta, x_pred = x_b,
    k = k_opt, ntree = ntree_opt, sparse = sparse_opt,
    nskip = n_burn, ndpost = n_post
  )
  pred_bart_b <- t(apply(ajuste_bart$draws, 2, quantile, probs = taus_crps))

  ## ---- Prediccion de los cuantiles para todos los origenes del bloque ----
  pred_fina$random_walk[idx_bloque, ]  <- matrix(cuantiles_rw, nrow = length(idx_bloque),
                                                  ncol = n_tf, byrow = TRUE)
  pred_fina$quantile_reg[idx_bloque, ] <- pred_rq_b
  pred_fina$quantile_rf[idx_bloque, ]  <- pred_qrf_b
  pred_fina$bart[idx_bloque, ]         <- pred_bart_b

  ## ---- Tiempo transcurrido y estimado restante ----
  transcurrido_min <- as.numeric(difftime(Sys.time(), t_inicio_eval, units = "mins"))
  restante_min     <- transcurrido_min / b * (n_bloques - b)
  cat(sprintf("[bloque %2d/%2d] origenes %4d-%4d  (ventana hasta %s)  ->  %.1f min transcurridos, ~%.1f min restantes\n",
              b, n_bloques, min(idx_bloque), max(idx_bloque),
              as.character(df_hasta$fecha[nrow(df_hasta)]),
              transcurrido_min, restante_min))

  rm(ajuste_bart); gc(verbose = FALSE)   # los objetos wbart pesan; liberar memoria
}

cat("\nEvaluacion recursiva por bloques completada. Tiempo total:",
    round(as.numeric(difftime(Sys.time(), t_inicio_eval, units = "mins")), 1), "min\n\n")


# ============================================================================
# 8. RESULTADOS DE LA EVALUACION
# ============================================================================
# Los 7 tau reportados se extraen de la grilla fina (que ya fue predicha).
idx_taus <- match(taus, taus_crps)

## ---- 8.1 Quantile Score por metodo y cuantil ----
tabla_quantile_score <- data.frame(
  tau          = taus,
  random_walk  = qs_metodo(pred_fina$random_walk[,  idx_taus, drop = FALSE], y_real, taus),
  quantile_reg = qs_metodo(pred_fina$quantile_reg[, idx_taus, drop = FALSE], y_real, taus),
  quantile_rf  = qs_metodo(pred_fina$quantile_rf[,  idx_taus, drop = FALSE], y_real, taus),
  bart         = qs_metodo(pred_fina$bart[,         idx_taus, drop = FALSE], y_real, taus)
)

cat("--- Quantile Score por metodo y cuantil (menor = mejor) ---\n")
print(round(tabla_quantile_score, 5), row.names = FALSE)

## ---- 8.2 qwCRPS por metodo, con los cuatro esquemas de ponderacion ----
tabla_qwcrps <- do.call(rbind, lapply(names(pred_fina), function(m) {
  data.frame(
    metodo           = m,
    qwcrps_colas     = qwcrps(pred_fina[[m]], y_real, peso = "ambas_colas"),
    qwcrps_izquierda = qwcrps(pred_fina[[m]], y_real, peso = "izquierda"),
    qwcrps_derecha   = qwcrps(pred_fina[[m]], y_real, peso = "derecha"),
    crps_plano       = qwcrps(pred_fina[[m]], y_real, peso = "plano")
  )
}))

cat("\n--- qwCRPS por metodo (menor = mejor) ---\n")
print(round_df <- (function(d) { d[-1] <- round(d[-1], 5); d })(tabla_qwcrps),
      row.names = FALSE)

## ---- 8.3 Desempenio RELATIVO al Random Walk ----
# Es la lectura estandar en la literatura cambiaria desde Meese y Rogoff (1983):
# el RW es un benchmark durisimo. Ratio < 1 => el metodo LE GANA al RW.
ratios <- tabla_quantile_score
ratios[, -1] <- sweep(ratios[, -1], 1, ratios$random_walk, "/")

cat("\n--- Quantile Score RELATIVO al Random Walk (<1 => mejor que el RW) ---\n")
print(round(ratios, 4), row.names = FALSE)

write.csv(tabla_quantile_score,
          file.path(ruta_out, "quantile_score_comparacion.csv"), row.names = FALSE)
write.csv(tabla_qwcrps, file.path(ruta_out, "qwcrps_comparacion.csv"), row.names = FALSE)
write.csv(ratios,
          file.path(ruta_out, "quantile_score_relativo_rw.csv"), row.names = FALSE)


# ============================================================================
# 9. DETERMINANTES: MODELO FINAL SOBRE TRAIN + VALIDACION
# ============================================================================
# Para el analisis de determinantes se estima el modelo con los hiperparametros
# ya tuneados sobre TODA la informacion previa al test (train + validacion =
# el 80% inicial). Asi los determinantes se identifican con la muestra mas
# grande posible sin tocar el bloque de evaluacion.

idx_det <- 1:corte_val
x_det   <- as.matrix(datos_modelo[idx_det, predictores])
y_det   <- datos_modelo$y_objetivo[idx_det]

cat("\n=== Estimando modelo final para determinantes (", length(idx_det),
    "obs.) ===\n")

modelo_det <- ajustar_bart_het(
  x_train = x_det, y_train = y_det, x_pred = x_det,
  k = k_opt, ntree = ntree_opt, sparse = sparse_opt,
  nskip = n_burn, ndpost = n_post
)

## ---- 9.1 VISTA (i): determinantes de la MEDIA ----
# Variables que desplazan f(x): efecto de NIVEL. Mueven TODOS los cuantiles
# por igual (suben o bajan la distribucion entera).
# Si sparse_opt = TRUE, ademas del varcount (frecuencia de uso) se dispone de
# varprob: la PROBABILIDAD DE INCLUSION POSTERIOR de cada predictor, que es la
# medida de importancia propia de DART y mas limpia que la simple frecuencia.
imp_media <- sort(modelo_det$varcount_media, decreasing = TRUE)
tabla_imp_media <- data.frame(variable = names(imp_media),
                              importancia_media = as.numeric(imp_media),
                              row.names = NULL)

if (!is.null(modelo_det$varprob_media)) {
  vp <- modelo_det$varprob_media
  tabla_imp_media$prob_inclusion_dart <- as.numeric(vp[tabla_imp_media$variable])
}

cat("\n--- (i) Determinantes de la MEDIA (efecto de NIVEL, todos los tau) ---\n")
print(head(tabla_imp_media, 12), row.names = FALSE)

## ---- 9.2 VISTA (ii): determinantes de la VARIANZA ----
# Variables que mueven sigma(x): efecto de DISPERSION. Ensanchan o estrechan
# la distribucion, es decir PESAN SOBRE TODO EN LAS COLAS. Esta es la vista
# que la version homoscedastica del script NO podia producir.
#
# DIAGNOSTICO CLAVE: si este varcount sale practicamente PLANO (todas las
# variables usadas con frecuencia similar), significa que sigma(x) es casi
# constante -> no se rechaza homocedasticidad -> los determinantes de cola y
# centro coinciden. ESO TAMBIEN ES UN RESULTADO VALIDO Y REPORTABLE, no un
# fracaso: seria evidencia de que los extremos del sol responden a los mismos
# drivers que los dias normales, solo que amplificados.
imp_var <- sort(modelo_det$varcount_var, decreasing = TRUE)
tabla_imp_var <- data.frame(variable = names(imp_var),
                            importancia_varianza = as.numeric(imp_var),
                            row.names = NULL)

cat("\n--- (ii) Determinantes de la VARIANZA (efecto de DISPERSION, las COLAS) ---\n")
print(head(tabla_imp_var, 12), row.names = FALSE)

# Coeficiente de variacion del varcount de la varianza: si es ~0, sigma(x) es
# practicamente constante (ver diagnostico arriba).
cv_var <- sd(modelo_det$varcount_var) / mean(modelo_det$varcount_var)
cat("\n[Diagnostico] Coef. de variacion del varcount de la varianza:",
    round(cv_var, 3),
    "\n  (cerca de 0 => sigma(x) casi constante => homocedasticidad no rechazada)\n")

write.csv(tabla_imp_media, file.path(ruta_out, "bart_importancia_media.csv"),   row.names = FALSE)
write.csv(tabla_imp_var,   file.path(ruta_out, "bart_importancia_varianza.csv"), row.names = FALSE)


# ============================================================================
# 10. DEPENDENCIA PARCIAL BAYESIANA E IMPORTANCIA POR CUANTIL
# ============================================================================
# Para cada predictor se recorre su rango empirico (percentiles 5, 25, 50, 75,
# 95) y se observa como se desplaza cada Q_tau.
#
# CORRECCION respecto de la version anterior: antes se fijaban los demas
# predictores en su MEDIANA. Eso es un corte tipo ICE, no la PDP de Friedman.
# La PDP correcta PROMEDIA la prediccion sobre los valores OBSERVADOS de las
# demas covariables. La diferencia importa justamente cuando hay interacciones,
# que es lo que BART esta disenado para capturar; fijar en la mediana puede
# distorsionar el efecto marginal.
#
# Implementacion: para cada valor v de la grilla se toma una submuestra de
# n_muestra_pdp filas reales, se les reemplaza SOLO la columna del predictor
# por v, y se AGRUPAN (pool) todos los draws predictivos resultantes. Los
# cuantiles de esa mezcla son la distribucion predictiva marginalizada.

set.seed(20222834)
filas_pdp <- sample(seq_len(nrow(x_det)), size = min(n_muestra_pdp, nrow(x_det)))
x_base_pdp <- x_det[filas_pdp, , drop = FALSE]

calcular_dependencia_parcial <- function(variable) {
  
  grilla_v <- quantile(x_det[, variable], probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
  
  res <- lapply(seq_along(grilla_v), function(g) {
    x_mod <- x_base_pdp
    x_mod[, variable] <- grilla_v[g]           # se altera SOLO esa columna
    
    draws <- reconstruir_draws(modelo_det, x_mod)   # ndpost x n_muestra_pdp
    q     <- quantile(as.numeric(draws), probs = taus)  # pool de todos los draws
    
    data.frame(variable = variable,
               valor_predictor = as.numeric(grilla_v[g]),
               t(q), row.names = NULL)
  })
  
  out <- do.call(rbind, res)
  colnames(out)[-(1:2)] <- paste0("tau_", taus)
  out
}

cat("\n=== Calculando dependencia parcial (", length(predictores),
    "predictores) ===\n")
dependencia_parcial <- map_dfr(predictores, calcular_dependencia_parcial)

write.csv(dependencia_parcial,
          file.path(ruta_out, "bart_dependencia_parcial_por_cuantil.csv"),
          row.names = FALSE)

## ---- 10.1 VISTA (iii): IMPORTANCIA POR CUANTIL ----
# Esta es la tabla que responde LA pregunta de la tesis.
# Para cada predictor y cada tau: AMPLITUD (max - min) del desplazamiento de
# Q_tau a lo largo del rango del predictor. Como Q_tau = f(x) + z_tau*sigma(x),
# esta amplitud YA depende de tau (en el modelo homoscedastico habria sido
# identica para todos los tau, por construccion).
#
# COMO LEER LA TABLA:
#   - Variable con importancia ALTA y PAREJA en todos los tau  -> driver de NIVEL
#   - Variable con importancia BAJA en tau=0.50 y ALTA en 0.01/0.99
#         -> driver ESPECIFICO DE LAS COLAS  <== el hallazgo buscado
#   - Ratio importancia(tau=0.01) / importancia(tau=0.50) >> 1 -> idem
importancia_por_cuantil <- dependencia_parcial %>%
  group_by(variable) %>%
  summarise(across(starts_with("tau_"), ~ max(.x) - min(.x)), .groups = "drop") %>%
  as.data.frame()

# Ratio cola/centro: el indicador sintetico del hallazgo
importancia_por_cuantil <- importancia_por_cuantil %>%
  mutate(
    ratio_cola_izq_centro = .data[["tau_0.01"]] / .data[["tau_0.5"]],
    ratio_cola_der_centro = .data[["tau_0.99"]] / .data[["tau_0.5"]]
  ) %>%
  arrange(desc(pmax(ratio_cola_izq_centro, ratio_cola_der_centro)))

cat("\n--- (iii) IMPORTANCIA POR CUANTIL (amplitud del efecto sobre Q_tau) ---\n")
cat("    Ordenada por ratio cola/centro: arriba, los drivers ESPECIFICOS DE COLA\n")
print(head(importancia_por_cuantil, 15), row.names = FALSE, digits = 3)

write.csv(importancia_por_cuantil,
          file.path(ruta_out, "bart_importancia_por_cuantil.csv"), row.names = FALSE)


# ============================================================================
# 11. GRAFICOS
# ============================================================================

## ---- 11.1 Importancia: media vs varianza (nivel vs dispersion) ----
comparacion_imp <- tabla_imp_media %>%
  select(variable, importancia_media) %>%
  left_join(tabla_imp_var, by = "variable") %>%
  pivot_longer(cols = c(importancia_media, importancia_varianza),
               names_to = "componente", values_to = "importancia") %>%
  mutate(componente = recode(componente,
                             importancia_media    = "Media f(x): efecto de NIVEL",
                             importancia_varianza = "Varianza sigma(x): efecto en COLAS"))

g_imp <- ggplot(comparacion_imp,
                aes(x = reorder(variable, importancia), y = importancia,
                    fill = componente)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  labs(title = "Determinantes: efecto de nivel vs. efecto de dispersion (BART)",
       subtitle = "Frecuencia promedio de uso en los arboles",
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(ruta_figs, "bart_importancia_media_vs_varianza.png"),
       g_imp, width = 9, height = 8)

## ---- 11.2 Heatmap de importancia por cuantil (el grafico clave) ----
heat <- importancia_por_cuantil %>%
  select(variable, starts_with("tau_")) %>%
  pivot_longer(-variable, names_to = "cuantil", values_to = "importancia") %>%
  group_by(variable) %>%
  mutate(importancia_rel = importancia / max(importancia)) %>%   # normaliza por fila
  ungroup()

g_heat <- ggplot(heat, aes(x = cuantil, y = reorder(variable, importancia_rel),
                           fill = importancia_rel)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = "Importancia de cada predictor POR CUANTIL",
       subtitle = "Amplitud del efecto sobre Q_tau, normalizada por variable",
       x = "Cuantil", y = NULL, fill = "Importancia\nrelativa") +
  theme_minimal()

ggsave(file.path(ruta_figs, "bart_importancia_por_cuantil_heatmap.png"),
       g_heat, width = 8, height = 9)

## ---- 11.3 Efecto marginal por cuantil, predictor por predictor ----
dependencia_larga <- dependencia_parcial %>%
  pivot_longer(cols = starts_with("tau_"), names_to = "cuantil", values_to = "valor_y")

g_dep <- ggplot(dependencia_larga,
                aes(x = valor_predictor, y = valor_y, color = cuantil)) +
  geom_line() + geom_point(size = 0.8) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(title = "Efecto marginal de cada predictor por cuantil (BART heteroscedastico)",
       subtitle = "Dependencia parcial promediada sobre la muestra",
       x = "Valor del predictor", y = "Cuantil condicional de r_(t+h)",
       color = "Cuantil") +
  theme_minimal()

ggsave(file.path(ruta_figs, "bart_dependencia_parcial_por_cuantil.png"),
       g_dep, width = 18, height = 16, limitsize = FALSE)

## ---- 11.4 Quantile Score relativo al Random Walk ----
g_qs <- ratios %>%
  select(-random_walk) %>%
  pivot_longer(-tau, names_to = "metodo", values_to = "ratio") %>%
  ggplot(aes(x = factor(tau), y = ratio, fill = metodo)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Quantile Score relativo al Random Walk",
       subtitle = "Por debajo de la linea = mejor que el Random Walk",
       x = "Cuantil (tau)", y = "QS / QS Random Walk", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(ruta_figs, "quantile_score_relativo_rw.png"),
       g_qs, width = 9, height = 5)


# ============================================================================
# 12. RESUMEN FINAL
# ============================================================================
cat("\n============================================================\n")
cat("RESUMEN\n")
cat("  Hiperparametros tuneados: k =", k_opt, "| ntree =", ntree_opt,
    "| sparse =", sparse_opt, "\n")
cat("  Origenes evaluados (ventana expansiva):", n_or, "\n")
cat("  Mejor metodo por qwCRPS (ambas colas):",
    tabla_qwcrps$metodo[which.min(tabla_qwcrps$qwcrps_colas)], "\n")
cat("\nArchivos guardados:\n")
cat("  Tablas ->", ruta_out, "\n")
cat("    - bart_tuneo_grilla.csv\n")
cat("    - quantile_score_comparacion.csv\n")
cat("    - quantile_score_relativo_rw.csv\n")
cat("    - qwcrps_comparacion.csv\n")
cat("    - bart_importancia_media.csv\n")
cat("    - bart_importancia_varianza.csv\n")
cat("    - bart_importancia_por_cuantil.csv\n")
cat("    - bart_dependencia_parcial_por_cuantil.csv\n")
cat("  Figuras ->", ruta_figs, "\n")
cat("    - bart_importancia_media_vs_varianza.png\n")
cat("    - bart_importancia_por_cuantil_heatmap.png\n")
cat("    - bart_dependencia_parcial_por_cuantil.png\n")
cat("    - quantile_score_relativo_rw.png\n")
cat("============================================================\n")

# Medir el costo real de cada pieza sobre una ventana típica
df_h <- datos_modelo[1:5508, ]
x_h  <- as.matrix(df_h[, predictores])
y_h  <- df_h$y_objetivo
df_1 <- datos_modelo[5509, , drop = FALSE]
x_1  <- as.matrix(df_1[, predictores])

system.time({ m_rq  <- rq(formula_rq, tau = taus_crps, data = df_h) })
system.time({ for(i in 1:20) predict(m_rq, newdata = df_1) })   # x20 días
system.time({ m_qrf <- quantregForest(x = x_h, y = y_h, ntree = 500) })
system.time({ for(i in 1:20) predict(m_qrf, newdata = x_1, what = taus_crps) })

# ============================================================================
# PENDIENTES / CHEQUEOS DE ROBUSTEZ SUGERIDOS
# ============================================================================
#  1. Test de Diebold-Mariano (o Giacomini-White) sobre las diferencias de
#     Quantile Score, para saber si la ventaja de BART sobre el RW es
#     ESTADISTICAMENTE significativa y no solo numerica.
#  2. Robustez del modelo heteroscedastico: reestimar con `rbart` (Pratola
#     et al., 2020), que muestrea media y varianza CONJUNTAMENTE, y comparar.
#  3. Repetir con h = 5 y h = 20 dias (el codigo ya lo soporta; el fix del
#     doble escalamiento del RW era condicion necesaria para esto).
#  4. Test de cobertura de los cuantiles (Kupiec / Christoffersen): verificar
#     que el 1% de las realizaciones efectivamente cae bajo Q_0.01.
# ============================================================================