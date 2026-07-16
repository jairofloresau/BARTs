############################################################################
##  ADAPTADOR QBART  ->  reemplaza el motor de estimacion de tu SCRIPT 02
##
##  Idea: tu pipeline (tuneo, ventana expansiva, qwCRPS, competidores,
##  determinantes) consume "una matriz de cuantiles por origen". Este archivo
##  produce EXACTAMENTE esa matriz, pero con BayesQArt (cuantil directo) en
##  lugar del BART heteroscedastico location-scale.
##
##  Motor: usa el paquete Rcpp original 'BayesQArt' (rapido) si esta instalado;
##  si no, cae al motor en R puro 'bayesqart.R'. Misma interfaz en ambos casos.
##
##  QBART ajusta UN modelo por cada tau -> es mas caro que tu enfoque de 2 BART.
##  Para el qwCRPS (que integra sobre una grilla fina de 99 tau) NO se ajustan
##  99 modelos: se ajusta una grilla moderada 'taus_q' y se interpola de forma
##  monotona hacia la grilla fina. Declararlo asi en la tesis.
############################################################################

if (!exists("bayesqart"))
  source("bayesqart.R")   # motor en R puro (fallback)

## Hiperparametros de QBART. Mapeo desde tu BART:
##   k (shrinkage)   ~  kappa   (controla la varianza del prior de nodos)
##   ntree           ~  m       (numero de arboles)
##   sparse (DART)   ->  no disponible en QBART aqui (ver nota al final)
hp_qbart_default <- list(
  m = 100L, burn = 1000L, ndpost = 2000L,
  nc = 100L, min_obs = 5L, kappa = 2.0, maxdepth = 3L,
  aa = 3.0, bb = 1.0            # prior IG(aa/2, bb/2) sobre phi
)

## --- Interruptor para IGNORAR el paquete compilado aunque este instalado ---
# Benchmark a escala real (n~6200,p=30,m=100): el motor Rcpp de 2016 resulto
# ~1.7x MAS LENTO por iteracion que bayesqart.R (0.82 vs 0.43 s/iter, y bajar
# 'nc' no lo arreglo -el costo no esta en la busqueda de cortes). Mientras eso
# no se optimice, conviene forzar el motor en R puro aunque BayesQArt este
# instalado. Poner FORZAR_MOTOR_R <- FALSE (antes de source() de este archivo)
# para volver a dejar la eleccion automatica.
if (!exists("FORZAR_MOTOR_R")) FORZAR_MOTOR_R <- FALSE

## --- Ajuste de UN cuantil. Devuelve prediccion en x_pred e importancia. -----
.fit_qbart_una <- function(x_train, y_train, x_pred, tau, hp = hp_qbart_default) {
  x_train <- as.matrix(x_train); x_pred <- as.matrix(x_pred)
  if (is.null(colnames(x_train)))
    colnames(x_train) <- colnames(x_pred) <- paste0("X", seq_len(ncol(x_train)))

  if (!FORZAR_MOTOR_R && requireNamespace("BayesQArt", quietly = TRUE)) {
    ## ---- motor Rcpp original (recomendado para la corrida final) ----
    # CORRECCION: (a) BayesQArt() no esta exportada por el paquete -> hace
    # falta ":::" en vez de "::", si no tira "not an exported object". (b) el
    # paquete llama a un C-callable interno de Rcpp (enterRNGScope) que solo
    # queda registrado si Rcpp esta CARGADO (library), no solo instalado; sin
    # esto tira "function 'enterRNGScope' not provided by package 'Rcpp'".
    library(Rcpp)
    out <- BayesQArt:::BayesQArt(
      y = y_train, X = x_train, Xtest = x_pred, quantile = tau,
      burn = hp$burn, nd = hp$ndpost, m = hp$m, min_obs_node = hp$min_obs,
      aa_parm = hp$aa, bb_parm = hp$bb, nc = hp$nc,
      pbd = 0.4, pb = 0.5, alpha = 0.95, betap = 2.0,
      kappa = hp$kappa, maxdepth = hp$maxdepth)
    pred <- out$pred_test
    imp  <- tabulate(out$vars_used + 1L, nbins = ncol(x_train))  # C++ es 0-indexado
    imp  <- setNames(if (sum(imp) > 0) imp / sum(imp) else imp, colnames(x_train))
  } else {
    ## ---- motor en R puro (bayesqart.R) ----
    out <- bayesqart(
      y = y_train, X = x_train, Xtest = x_pred, tau = tau,
      n_trees = hp$m, burn = hp$burn, n_draws = hp$ndpost, nc = hp$nc,
      min_obs = hp$min_obs, kappa = hp$kappa, maxdepth = hp$maxdepth,
      aa = hp$aa, bb = hp$bb, verbose = 0)
    pred <- out$pred_test
    imp  <- out$var_importance
  }
  list(pred = as.numeric(pred), imp = imp)
}


## ==========================================================================
## FUNCION 1: matriz de cuantiles  (reemplaza a ajustar_bart_het + draws)
## ==========================================================================
## Devuelve una matriz  n_pred x length(taus_q)  con los cuantiles condicionales
## en x_pred. Aplica el reordenamiento de Chernozhukov et al. (2010) por fila
## para evitar cruces (cada tau es un modelo aparte y podria cruzarse).
predecir_qbart_grid <- function(x_train, y_train, x_pred, taus_q,
                                hp = hp_qbart_default) {
  x_pred <- as.matrix(x_pred)
  Q <- vapply(taus_q,
              function(tt) .fit_qbart_una(x_train, y_train, x_pred, tt, hp)$pred,
              numeric(nrow(x_pred)))
  if (is.null(dim(Q))) Q <- matrix(Q, nrow = nrow(x_pred))
  Q <- t(apply(Q, 1, sort))                      # anti-cruce
  colnames(Q) <- paste0("tau_", taus_q)
  Q
}

## Interpola de forma monotona una matriz de cuantiles de 'taus_q' hacia la
## grilla fina 'taus_fine' que necesita tu qwCRPS(). Lineal en tau, extrapola
## plano en los extremos (rule = 2). Como cada fila entra ordenada, el resultado
## se mantiene monotono.
expandir_a_grilla_fina <- function(Q_coarse, taus_q, taus_fine) {
  t(apply(Q_coarse, 1, function(qc)
    stats::approx(x = taus_q, y = qc, xout = taus_fine, rule = 2)$y))
}


## ==========================================================================
## FUNCION 2: determinantes  ->  importancia POR CUANTIL (nativa en QBART)
## ==========================================================================
## Reemplaza al desglose media/varianza de tu §9. Ajusta QBART en cada tau y
## devuelve la tabla  variable x tau  con la frecuencia de uso en las
## divisiones. Comparar columnas responde directamente la pregunta de la tesis.
importancia_qbart_por_cuantil <- function(x, y, taus_q, hp = hp_qbart_default) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  x_dummy <- x[seq_len(min(5, nrow(x))), , drop = FALSE]

  cols <- lapply(taus_q, function(tt) {
    r <- .fit_qbart_una(x, y, x_dummy, tt, hp)
    r$imp[colnames(x)]
  })
  M <- do.call(cbind, cols)
  colnames(M) <- paste0("tau_", taus_q); rownames(M) <- colnames(x)
  tab <- data.frame(variable = rownames(M), M, row.names = NULL,
                    check.names = FALSE)
  ## razon cola/centro: >1 => la variable pesa mas en la cola que en el centro
  centro <- paste0("tau_", taus_q[which.min(abs(taus_q - 0.5))])
  for (col in grep("^tau_", names(tab), value = TRUE))
    if (col != centro)
      tab[[paste0("ratio_", col, "_centro")]] <- tab[[col]] / tab[[centro]]
  tab
}
