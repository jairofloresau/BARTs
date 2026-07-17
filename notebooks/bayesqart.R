############################################################################
##  BayesQArt: motor local EXPERIMENTAL
##  Reimplementacion aproximada en R puro del metodo de Kindo et al. (2016).
##  No usar para resultados de investigacion: el flujo principal exige el
##  kernel C++ oficial. Este archivo se conserva para desarrollo y docencia.
##  "Bayesian quantile additive regression trees", arXiv:1607.02676
##
##  Aplicacion: pronostico e identificacion de determinantes de los retornos
##  extremos del tipo de cambio (estimando cuantiles condicionales tau).
##
##  Estructura del modelo (ec. 5 del paper):
##     y_i = G(x_i; T,M) + v1 * nu_i + v2 * sqrt(phi * nu_i) * z_i
##  con error Laplace asimetrica (ALD) representada como mezcla
##  normal-exponencial, donde:
##     v1  = (1 - 2*tau) / (tau*(1-tau))          [vartheta_1]
##     v2^2 = 2 / (tau*(1-tau))                    [vartheta_2^2]
##  Esquema MCMC (Seccion 2.2): (i) muestrear latentes nu ~ GIG,
##  (ii) actualizar cada arbol por Metropolis-Hastings (GROW/PRUNE/CHANGE/SWAP)
##  y sus nodos terminales ~ Normal, (iii) muestrear la escala phi ~ IG.
############################################################################


## ============================================================================
## 1. UTILIDADES DE MUESTREO
## ============================================================================

## Inversa Gaussiana IG(mu, lambda) por el metodo de Michael-Schucany-Haas (1976).
## Vectorizada en mu (lambda escalar).
.rinvgauss <- function(mu, lambda) {
  n  <- length(mu)
  nu <- rnorm(n); y <- nu * nu
  x  <- mu + (mu * mu * y) / (2 * lambda) -
        (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu * mu * y * y)
  u  <- runif(n)
  ifelse(u <= mu / (mu + x), x, (mu * mu) / x)
}

## Generalized Inverse Gaussian GIG(lambda = 1/2, chi, psi).
## Se usa la propiedad de reciprocidad: si X ~ GIG(1/2, chi, psi) entonces
## 1/X ~ GIG(-1/2, psi, chi) = InvGauss(mu = sqrt(psi/chi), lambda = psi).
## chi puede ser un vector (psi escalar). Este es exactamente el latente nu.
.rgig_half <- function(chi, psi) {
  chi <- pmax(chi, 1e-12)          # evita division por cero si el residuo es 0
  mu  <- sqrt(psi / chi)
  1 / .rinvgauss(mu, lambda = psi)
}


## ============================================================================
## 2. PUNTOS DE CORTE (grilla de cutpoints por variable, como en BART)
## ============================================================================

## Devuelve una lista: para cada columna de X, 'nc' cortes candidatos.
make_cutpoints <- function(X, nc) {
  lapply(seq_len(ncol(X)), function(j) {
    rng <- range(X[, j])
    if (diff(rng) == 0) return(rng[1])          # variable constante
    seq(rng[1], rng[2], length.out = nc + 2L)[-c(1L, nc + 2L)]
  })
}


## ============================================================================
## 3. ESTRUCTURA DE ARBOL  (indexacion tipo heap: hijos de i -> 2i y 2i+1)
## ============================================================================
## Un arbol es una lista con $nodes: lista nombrada por id (character).
## Cada nodo: list(id, leaf, var, cut, mu, depth). Nodo interno: leaf=FALSE con
## var/cut; nodo terminal: leaf=TRUE con mu. La raiz tiene id = 1, depth = 0.

new_tree <- function(mu0 = 0) {
  list(nodes = list("1" = list(id = 1L, leaf = TRUE,
                               var = NA_integer_, cut = NA_real_,
                               mu = mu0, depth = 0L)))
}

.depth   <- function(id) as.integer(floor(log2(id)))
.terminal_ids <- function(tree) {
  ids <- as.integer(names(tree$nodes))
  ids[vapply(tree$nodes, function(nd) nd$leaf, logical(1))]
}
.internal_ids <- function(tree) {
  ids <- as.integer(names(tree$nodes))
  ids[vapply(tree$nodes, function(nd) !nd$leaf, logical(1))]
}

## Asigna cada fila de X a su nodo terminal (descenso vectorizado por sub-indices).
assign_terminal <- function(tree, X) {
  n <- nrow(X); out <- integer(n)
  stack <- list(list(id = 1L, idx = seq_len(n)))
  while (length(stack) > 0) {
    cur <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
    nd  <- tree$nodes[[as.character(cur$id)]]
    if (nd$leaf) { out[cur$idx] <- cur$id; next }
    goleft <- X[cur$idx, nd$var] < nd$cut
    li <- cur$idx[goleft]; ri <- cur$idx[!goleft]
    if (length(li)) stack[[length(stack) + 1L]] <- list(id = 2L * cur$id,      idx = li)
    if (length(ri)) stack[[length(stack) + 1L]] <- list(id = 2L * cur$id + 1L, idx = ri)
  }
  out
}

## Ajuste del arbol: valor mu del nodo terminal de cada observacion.
fit_tree <- function(tree, X) {
  leaf_of <- assign_terminal(tree, X)
  mus <- vapply(tree$nodes, function(nd) if (nd$leaf) nd$mu else NA_real_, numeric(1))
  names(mus) <- names(tree$nodes)
  as.numeric(mus[as.character(leaf_of)])
}


## ============================================================================
## 4. VEROSIMILITUD MARGINAL POR NODO (integrando mu, ec. 15 del paper)
## ============================================================================
## Para un nodo terminal con residuos-pseudo w y latentes nu:
##   A = sum(w / nu),  B = sum(1 / nu)
## log del factor de integracion de mu (los terminos que NO dependen de la
## particion se cancelan en el cociente MH, por eso solo guardamos este):
##   0.5*log( v2sq*phi / (v2sq*phi + s0sq*B) ) + s0sq*A^2 / (2*v2sq*phi*(...))
.node_logg <- function(A, B, phi, v2sq, s0sq) {
  denom <- v2sq * phi + s0sq * B
  0.5 * log((v2sq * phi) / denom) + (s0sq * A * A) / (2 * v2sq * phi * denom)
}

## Estadisticos (A, B, n) de cada nodo terminal, dados w (pseudo-residuo) y nu.
.leaf_stats <- function(tree, X, w, nu) {
  leaf_of <- assign_terminal(tree, X)
  inv_nu  <- 1 / nu
  ids <- .terminal_ids(tree)
  A <- B <- nn <- setNames(numeric(length(ids)), as.character(ids))
  for (id in ids) {
    sel <- leaf_of == id
    k <- as.character(id)
    B[k]  <- sum(inv_nu[sel])
    A[k]  <- sum(w[sel] * inv_nu[sel])
    nn[k] <- sum(sel)
  }
  list(A = A, B = B, n = nn, leaf_of = leaf_of)
}

## Score marginal del arbol completo = suma de .node_logg sobre nodos terminales.
.tree_score <- function(stats, phi, v2sq, s0sq) {
  sum(.node_logg(stats$A, stats$B, phi, v2sq, s0sq))
}


## ============================================================================
## 5. MOVIMIENTOS METROPOLIS-HASTINGS  (GROW, PRUNE, CHANGE, SWAP)
## ============================================================================
## Probabilidad de que un nodo a profundidad d se divida: alpha/(1+d)^beta.
.psplit <- function(d, alpha, beta) alpha / (1 + d)^beta

## Nodos terminales que se pueden dividir (profundidad < maxdepth y con
## suficientes observaciones para respetar min_obs en ambos hijos).
.growable <- function(tree, stats, maxdepth, min_obs) {
  ids <- .terminal_ids(tree)
  ids[vapply(ids, function(id)
    .depth(id) < maxdepth && stats$n[as.character(id)] >= 2L * min_obs,
    logical(1))]
}

## Nodos internos "podables": ambos hijos son terminales.
.prunable <- function(tree) {
  ii <- .internal_ids(tree)
  ii[vapply(ii, function(id) {
    l <- as.character(2L * id); r <- as.character(2L * id + 1L)
    !is.null(tree$nodes[[l]]) && !is.null(tree$nodes[[r]]) &&
      tree$nodes[[l]]$leaf && tree$nodes[[r]]$leaf
  }, logical(1))]
}

## Propone un corte valido (variable, valor) en un nodo terminal 'id' tal que
## ambos hijos reciban al menos min_obs observaciones. Devuelve NULL si no hay.
.propose_split <- function(tree, X, xi, id, min_obs) {
  idx  <- assign_terminal(tree, X) == id
  Xsub <- X[idx, , drop = FALSE]
  p    <- ncol(X)
  vars <- sample(seq_len(p))                       # orden aleatorio de variables
  for (v in vars) {
    cuts <- xi[[v]]
    if (length(cuts) == 0) next
    cuts <- sample(cuts)                            # cortes en orden aleatorio
    for (c in cuts) {
      nl <- sum(Xsub[, v] <  c); nr <- sum(Xsub[, v] >= c)
      if (nl >= min_obs && nr >= min_obs)
        return(list(var = v, cut = c))
    }
  }
  NULL
}

## --- GROW ---------------------------------------------------------------
.move_grow <- function(tree, X, w, nu, xi, pr, stats) {
  gid <- .growable(tree, stats, pr$maxdepth, pr$min_obs)
  if (length(gid) == 0) return(NULL)
  id  <- gid[sample.int(length(gid), 1L)]
  spl <- .propose_split(tree, X, xi, id, pr$min_obs)
  if (is.null(spl)) return(NULL)

  ## arbol propuesto
  tp <- tree
  k  <- as.character(id); d <- tree$nodes[[k]]$depth
  tp$nodes[[k]]$leaf <- FALSE; tp$nodes[[k]]$mu <- NA_real_
  tp$nodes[[k]]$var  <- spl$var; tp$nodes[[k]]$cut <- spl$cut
  tp$nodes[[as.character(2L * id)]]     <- list(id = 2L * id,      leaf = TRUE,
                                    var = NA, cut = NA, mu = 0, depth = d + 1L)
  tp$nodes[[as.character(2L * id + 1L)]] <- list(id = 2L * id + 1L, leaf = TRUE,
                                    var = NA, cut = NA, mu = 0, depth = d + 1L)
  sp <- .leaf_stats(tp, X, w, nu)

  ## cociente de verosimilitud (solo cambia el nodo dividido)
  kp <- as.character(id); kl <- as.character(2L*id); kr <- as.character(2L*id+1L)
  ll <- .node_logg(sp$A[kl], sp$B[kl], pr$phi, pr$v2sq, pr$s0sq) +
        .node_logg(sp$A[kr], sp$B[kr], pr$phi, pr$v2sq, pr$s0sq) -
        .node_logg(stats$A[kp], stats$B[kp], pr$phi, pr$v2sq, pr$s0sq)

  ## cociente de prior de estructura
  ps_d  <- .psplit(d,     pr$alpha, pr$beta)
  ps_d1 <- .psplit(d + 1, pr$alpha, pr$beta)
  lp    <- log(ps_d) + 2 * log(1 - ps_d1) - log(1 - ps_d)

  ## cociente de transicion  (P_prune/P_grow) * (#growables / #podables tras crecer)
  b  <- length(gid)
  w2 <- length(.prunable(tp))
  lt <- log(pr$p_prune / pr$p_grow) + log(b / w2)

  list(tree = tp, logratio = ll + lp + lt, stats = sp)
}

## --- PRUNE --------------------------------------------------------------
.move_prune <- function(tree, X, w, nu, xi, pr, stats) {
  pid <- .prunable(tree)
  if (length(pid) == 0) return(NULL)
  id  <- pid[sample.int(length(pid), 1L)]
  k   <- as.character(id); d <- tree$nodes[[k]]$depth
  kl  <- as.character(2L*id); kr <- as.character(2L*id+1L)

  tp <- tree
  tp$nodes[[k]]$leaf <- TRUE; tp$nodes[[k]]$mu <- 0
  tp$nodes[[k]]$var  <- NA_integer_; tp$nodes[[k]]$cut <- NA_real_
  tp$nodes[[kl]] <- NULL; tp$nodes[[kr]] <- NULL
  sp <- .leaf_stats(tp, X, w, nu)

  ll <- .node_logg(sp$A[k], sp$B[k], pr$phi, pr$v2sq, pr$s0sq) -
        .node_logg(stats$A[kl], stats$B[kl], pr$phi, pr$v2sq, pr$s0sq) -
        .node_logg(stats$A[kr], stats$B[kr], pr$phi, pr$v2sq, pr$s0sq)

  ps_d  <- .psplit(d,     pr$alpha, pr$beta)
  ps_d1 <- .psplit(d + 1, pr$alpha, pr$beta)
  lp    <- log(1 - ps_d) - log(ps_d) - 2 * log(1 - ps_d1)

  w1 <- length(pid)                                   # podables antes
  b2 <- length(.growable(tp, sp, pr$maxdepth, pr$min_obs))  # growables despues
  lt <- log(pr$p_grow / pr$p_prune) + log(w1 / b2)

  list(tree = tp, logratio = ll + lp + lt, stats = sp)
}

## --- CHANGE  (cambia variable/valor de corte en un nodo interno) --------
## Propuesta simetrica => el cociente de prior/transicion es 1; se acepta con
## min(1, exp(delta score)). Se rechaza si crea nodos con < min_obs.
.move_change <- function(tree, X, w, nu, xi, pr, stats) {
  ii <- .internal_ids(tree)
  if (length(ii) == 0) return(NULL)
  id <- ii[sample.int(length(ii), 1L)]
  k  <- as.character(id)

  ## conjunto de observaciones que caen bajo el subarbol de 'id'
  reach <- .subtree_idx(tree, X, id)
  Xsub  <- X[reach, , drop = FALSE]
  spl   <- .propose_split_free(xi, Xsub, ncol(X))    # variable/cut cualquiera
  if (is.null(spl)) return(NULL)

  tp <- tree
  tp$nodes[[k]]$var <- spl$var; tp$nodes[[k]]$cut <- spl$cut
  ## validez: todos los terminales del subarbol deben conservar min_obs
  sp <- .leaf_stats(tp, X, w, nu)
  if (any(sp$n < pr$min_obs)) return(NULL)

  ll <- .tree_score(sp, pr$phi, pr$v2sq, pr$s0sq) -
        .tree_score(stats, pr$phi, pr$v2sq, pr$s0sq)
  list(tree = tp, logratio = ll, stats = sp)
}

## --- SWAP  (intercambia la regla de un nodo interno con la de un hijo interno) -
.move_swap <- function(tree, X, w, nu, xi, pr, stats) {
  ii <- .internal_ids(tree)
  ## padres internos con al menos un hijo interno
  cand <- ii[vapply(ii, function(id) {
    l <- as.character(2L*id); r <- as.character(2L*id+1L)
    (!is.null(tree$nodes[[l]]) && !tree$nodes[[l]]$leaf) ||
    (!is.null(tree$nodes[[r]]) && !tree$nodes[[r]]$leaf)
  }, logical(1))]
  if (length(cand) == 0) return(NULL)
  id <- cand[sample.int(length(cand), 1L)]
  l  <- 2L*id; r <- 2L*id + 1L
  ichild <- c(l, r)[c(!is.null(tree$nodes[[as.character(l)]]) &&
                        !tree$nodes[[as.character(l)]]$leaf,
                      !is.null(tree$nodes[[as.character(r)]]) &&
                        !tree$nodes[[as.character(r)]]$leaf)]
  cid <- ichild[sample.int(length(ichild), 1L)]

  tp <- tree; kp <- as.character(id); kc <- as.character(cid)
  tmp_v <- tp$nodes[[kp]]$var; tmp_c <- tp$nodes[[kp]]$cut
  tp$nodes[[kp]]$var <- tp$nodes[[kc]]$var; tp$nodes[[kp]]$cut <- tp$nodes[[kc]]$cut
  tp$nodes[[kc]]$var <- tmp_v;              tp$nodes[[kc]]$cut <- tmp_c

  sp <- .leaf_stats(tp, X, w, nu)
  if (any(sp$n < pr$min_obs)) return(NULL)
  ll <- .tree_score(sp, pr$phi, pr$v2sq, pr$s0sq) -
        .tree_score(stats, pr$phi, pr$v2sq, pr$s0sq)
  list(tree = tp, logratio = ll, stats = sp)
}

## Indices de X que alcanzan el subarbol con raiz 'id'.
.subtree_idx <- function(tree, X, id) {
  n <- nrow(X); reach <- logical(n)
  stack <- list(list(id = 1L, idx = seq_len(n)))
  while (length(stack) > 0) {
    cur <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
    if (cur$id == id) { reach[cur$idx] <- TRUE; next }
    nd <- tree$nodes[[as.character(cur$id)]]
    if (nd$leaf) next
    goleft <- X[cur$idx, nd$var] < nd$cut
    stack[[length(stack)+1L]] <- list(id = 2L*cur$id,      idx = cur$idx[goleft])
    stack[[length(stack)+1L]] <- list(id = 2L*cur$id + 1L, idx = cur$idx[!goleft])
  }
  which(reach)
}

## Corte candidato libre (para CHANGE): variable y valor uniformes que dividan Xsub.
.propose_split_free <- function(xi, Xsub, p) {
  for (v in sample(seq_len(p))) {
    cuts <- xi[[v]]; if (length(cuts) == 0) next
    c <- cuts[sample.int(length(cuts), 1L)]
    if (sum(Xsub[, v] < c) > 0 && sum(Xsub[, v] >= c) > 0)
      return(list(var = v, cut = c))
  }
  NULL
}


## ============================================================================
## 6. ACTUALIZAR PARAMETROS DE NODOS TERMINALES  (drmu, ec. 19)
## ============================================================================
.draw_mu <- function(tree, stats, pr) {
  for (id in .terminal_ids(tree)) {
    k <- as.character(id)
    A <- stats$A[k]; B <- stats$B[k]
    prec <- 1 / pr$s0sq + B / (pr$v2sq * pr$phi)   # precision posterior
    post_var  <- 1 / prec
    post_mean <- post_var * A / (pr$v2sq * pr$phi)
    tree$nodes[[k]]$mu <- rnorm(1, post_mean, sqrt(post_var))
  }
  tree
}


## ============================================================================
## 7. FUNCION PRINCIPAL
## ============================================================================
#' Ajusta BayesQArt para un cuantil tau y predice sobre datos de prueba.
#'
#' @param y         vector respuesta (p.ej. retornos del tipo de cambio)
#' @param X         matriz de predictores de entrenamiento (n x p)
#' @param Xtest     matriz de predictores de prueba (np x p); puede ser NULL
#' @param tau       cuantil a estimar en (0,1), p.ej. 0.05, 0.50, 0.95
#' @param n_trees   numero de arboles en la suma (m)
#' @param burn      iteraciones de calentamiento
#' @param n_draws   iteraciones post-burn que se promedian
#' @param nc        numero de cortes candidatos por variable
#' @param min_obs   minimo de observaciones por nodo terminal
#' @param alpha,beta  hiperparametros del prior de arbol  (psi1, psi2)
#' @param kappa     controla la varianza del prior de nodos: s0^2 = 1/(2*kappa*sqrt(m))
#' @param maxdepth  profundidad maxima de cada arbol
#' @param aa,bb     hiperparametros del prior IG(aa/2, bb/2) sobre phi
#' @param verbose   imprime progreso cada 'verbose' iteraciones (0 = silencio)
#' @return lista con pred_train, pred_test, var_importance, y metadatos.
bayesqart <- function(y, X, Xtest = NULL, tau = 0.5,
                      n_trees = 50L, burn = 1000L, n_draws = 2000L,
                      nc = 50L, min_obs = 5L,
                      alpha = 0.95, beta = 2.0, kappa = 2.0, maxdepth = 3L,
                      aa = 2.0, bb = 3.0, verbose = 200L) {

  X <- as.matrix(X); y <- as.numeric(y)
  n <- length(y); p <- ncol(X)
  stopifnot(nrow(X) == n, tau > 0, tau < 1)

  ## constantes de la ALD
  v1   <- (1 - 2 * tau) / (tau * (1 - tau))     # vartheta_1
  v2sq <- 2 / (tau * (1 - tau))                 # vartheta_2^2

  ## escalado de la respuesta a [-0.5, 0.5]  (invariante al centrado; ver paper)
  ymin <- min(y); yrange <- max(y) - min(y)
  ys   <- -0.5 + (y - ymin) / yrange

  ## grilla de cortes
  xi <- make_cutpoints(X, nc)

  ## hiperparametros y estado en un contenedor 'pr'
  pr <- list(v1 = v1, v2sq = v2sq, alpha = alpha, beta = beta,
             s0sq = 1 / (2 * kappa * sqrt(n_trees)),
             maxdepth = maxdepth, min_obs = min_obs, phi = 1.0,
             p_grow = 0.25, p_prune = 0.25, p_change = 0.10, p_swap = 0.40)

  ## inicializacion: m arboles-raiz con mu = 0
  trees  <- replicate(n_trees, new_tree(0), simplify = FALSE)
  allfit <- rep(0, n)                            # suma de ajustes de todos los arboles

  ## acumuladores de salida
  pmean  <- rep(0, n)
  ppred  <- if (!is.null(Xtest)) rep(0, nrow(as.matrix(Xtest))) else NULL
  Xtest  <- if (!is.null(Xtest)) as.matrix(Xtest) else NULL
  var_counts <- rep(0, p)                         # importancia = # de cortes por variable

  psi_lat <- (2 + v1 * v1 / v2sq)                # numerador de psi (falta /phi)
  total   <- burn + n_draws

  for (it in seq_len(total)) {
    if (verbose > 0 && it %% verbose == 0)
      cat(sprintf("  [tau=%.2f] iteracion %d/%d\n", tau, it, total))

    ## (i) muestrear latentes nu ~ GIG(1/2, chi, psi)  (una vez por iteracion)
    chi <- (ys - allfit)^2 / (pr$phi * v2sq)
    psi <- psi_lat / pr$phi
    nu  <- .rgig_half(chi, psi)

    ## (ii) actualizar cada arbol por Metropolis-Hastings + sus nodos terminales
    for (j in seq_len(n_trees)) {
      fj      <- fit_tree(trees[[j]], X)
      allfit  <- allfit - fj                     # quitar el arbol j
      w       <- ys - allfit - v1 * nu           # pseudo-residuo (ec. 12)

      stats   <- .leaf_stats(trees[[j]], X, w, nu)

      u <- runif(1)
      prop <- if (u < pr$p_grow) {
                .move_grow(trees[[j]], X, w, nu, xi, pr, stats)
              } else if (u < pr$p_grow + pr$p_prune) {
                .move_prune(trees[[j]], X, w, nu, xi, pr, stats)
              } else if (u < pr$p_grow + pr$p_prune + pr$p_change) {
                .move_change(trees[[j]], X, w, nu, xi, pr, stats)
              } else {
                .move_swap(trees[[j]], X, w, nu, xi, pr, stats)
              }

      if (!is.null(prop) && is.finite(prop$logratio) &&
          log(runif(1)) < prop$logratio) {
        trees[[j]] <- prop$tree
        stats      <- prop$stats
      }

      ## muestrear mu de los nodos terminales
      trees[[j]] <- .draw_mu(trees[[j]], stats, pr)

      allfit <- allfit + fit_tree(trees[[j]], X) # volver a sumar el arbol j
    }

    ## (iii) muestrear la escala phi ~ IG  (ec. 20)
    resid <- ys - allfit - v1 * nu
    ss    <- sum(resid * resid / (2 * v2sq * nu))
    pr$phi <- 1 / rgamma(1, shape = n / 2 + aa / 2, rate = (bb + ss) / 2)

    ## acumular tras el burn-in
    if (it > burn) {
      pmean <- pmean + allfit
      if (!is.null(Xtest)) {
        pt <- rep(0, nrow(Xtest))
        for (j in seq_len(n_trees)) pt <- pt + fit_tree(trees[[j]], Xtest)
        ppred <- ppred + pt
      }
      for (j in seq_len(n_trees))
        for (id in .internal_ids(trees[[j]]))
          var_counts[trees[[j]]$nodes[[as.character(id)]]$var] <-
            var_counts[trees[[j]]$nodes[[as.character(id)]]$var] + 1
    }
  }

  ## promediar y des-escalar a la escala original de y
  pmean <- pmean / n_draws
  unscale <- function(s) ymin + (s + 0.5) * yrange
  out <- list(
    tau        = tau,
    pred_train = unscale(pmean),
    pred_test  = if (!is.null(Xtest)) unscale(ppred / n_draws) else NULL,
    var_importance = setNames(var_counts / sum(var_counts),
                              if (!is.null(colnames(X))) colnames(X)
                              else paste0("X", seq_len(p))),
    phi = pr$phi
  )
  out
}


## ============================================================================
## 8. EVALUACION: check loss / desviacion absoluta media ponderada (MWAD, ec. 21)
## ============================================================================
## Es la funcion de perdida propia de la regresion cuantilica; menor es mejor.
wmad <- function(y_true, y_hat, tau) {
  d <- y_true - y_hat
  mean(tau * abs(d) * (d > 0) + (1 - tau) * abs(d) * (d <= 0))
}
