#' @title Multiscale Geographically Weighted Negative Binomial Regression
#'
#' @description Fits a geographically weighted regression model with different scales for each covariate. Uses the negative binomial distribution as default, but also accepts the normal, Poisson, or logistic distributions. Can fit the global versions of each regression and also the geographically weighted alternatives with only one scale, since they are all particular cases of the multiscale approach.
#'
#' @param data name of the dataset.
#' @param formula regression model formula as in \code{lm}.
#' @param weight name of the variable containing the sample weights, default value is \code{NULL}.
#' @param lat name of the variable containing the latitudes in the dataset.
#' @param long name of the variable containing the longitudes in the dataset.
#' @param globalmin logical value indicating whether to find a global minimum in the optimization process, default value is \code{TRUE}.
#' @param method indicates the method to be used for the bandwidth calculation (\code{adaptive_bsq}, \code{fixed_bsq}, \code{fixed_g}).
#' @param model indicates the model to be used for the regression (\code{gaussian}, \code{poisson}, \code{negbin}, \code{logistic}), default value is\code{"negbin"}.
#' @param mgwr logical value indicating if multiscale should be used (\code{TRUE}, \code{FALSE}), default value is \code{TRUE}.
#' @param bandwidth indicates the criterion to be used for the bandwidth calculation (\code{cv}, \code{aic}), default value is \code{"cv"}.
#' @param offset name of the variable containing the offset values, if null then is set to a vector of zeros, default value is \code{NULL}.
#' @param distancekm logical value indicating whether to calculate the distances in km, default value is \code{FALSE}.
#' @param int integer indicating the number of iterations, default value is \code{50}.
#' @param h integer indicating a predetermined bandwidth value, default value is \code{NULL}.
#'
#' @return A list that contains:
#'
#' \itemize{
#' \item \code{general_bandwidth} - General bandwidth value.
#' \item \code{band} - Bandwidth values for each covariate.
#' \item \code{measures} - Goodness of fit statistics.
#' \item \code{ENP} - Effective number of parameters.
#' \item \code{mgwr_param_estimates} - MGWR parameter estimates.
#' \item \code{qntls_mgwr_param_estimates} - Quantiles of MGWR parameter estimates.
#' \item \code{descript_stats_mgwr_param_estimates} - Descriptive statistics of MGWR parameter estimates.
#' \item \code{p_values} - P-values for the t tests on parameter significance.
#' \item \code{t_critical} - Critical values for the t tests on parameter significance.
#' \item \code{mgwr_se} - MGWR standard errors.
#' \item \code{qntls_mgwr_se} - Quantiles of MGWR standard errors.
#' \item \code{descript_stats_se} - Descriptive statistics of MGWR standard errors.
#' \item \code{global_param_estimates} - Parameter estimates for the global model.
#' \item \code{t_test_dfs} - Denominator degrees of freedom for the t tests.
#' \item \code{global_measures} - Goodness of fit statistics for the global model.
#' }
#'
#' @examples
#' ## Data
#'
#'
#' data(georgia)
#'
#' for (var in c("PctFB", "PctBlack")){
#'   georgia[, var] <- as.data.frame(scale(georgia[, var]))
#' }
#'
#'
#' ## Model
#'
#' mod <- mgwnbr(data=georgia, formula=PctBach~PctBlack+PctFB,
#'  lat="Y", long="X", globalmin=FALSE, method="adaptive_bsq", bandwidth="cv",
#'   model="gaussian", mgwr=FALSE, h=136)
#'
#' ## Bandwidths
#' mod$general_bandwidth
#'
#' ## Goodness of fit measures
#' mod$measures
#'
#' @importFrom sp spDistsN1
#'
#' @importFrom stats model.extract model.matrix pnorm pt qt quantile dist
#'
#'
#' @export

# Pure local version of gwr for alpha optimization
gwr_local <- function(H, y, x, fi, alpha_val = 0.5, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta, verbose = FALSE) {
  nvar <- ncol(x)
  bim <- rep(0, nvar*N)
  yhatm <- rep(0, N)
  sm <- matrix(0, N, N)
  sm3 <- matrix(0, N, nvar)
  alphai <- matrix(0, N, 3)
  # No assign or get to parent/global, all local!
  for (i in 1:N){
    for (j in 1:N){
      seqi <- rep(i, N)
      distan <- cbind(seqi, sequ, sp::spDistsN1(COORD, COORD[i,]))
      if (distancekm){
        distan[,3] <- distan[,3]*111
      }
    }
    u <- nrow(distan)
    w <- rep(0, u)
    if (method=="fixed_g"){
      for (jj in 1:u){
        w[jj] <- exp(-(distan[jj,3]/H)^2)
      }
    }
    else if (method=="fixed_bsq"){
      for (jj in 1:u){
        w[jj] <- (1-(distan[jj,3]/H)^2)^2
      }
    }
    else if (method=="adaptive_bsq"){
      distan <- distan[order(distan[, 3]), ]
      distan <- cbind(distan, 1:nrow(distan))
      w <- matrix(0, N, 2)
      hn <- distan[H,3]
      for (jj in 1:N){
        if (distan[jj,4]<=H){
          w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
        } else{
          w[jj,1] <- 0
        }
        w[jj,2] <- distan[jj,2]
      }
      w <- w[order(w[, 2]), ]
      w <- w[,1]
    }
    else if (method == "adaptive_bsq_smr") {
      beta_val <- 1 - alpha_val
      # --- w1 ---
      dist_spatial_sorted <- distan[order(distan[, 3]), ]
      dist_spatial_sorted <- cbind(dist_spatial_sorted, 1:nrow(dist_spatial_sorted))
      w1_matrix <- matrix(0, N, 2)
      hn_spatial <- dist_spatial_sorted[H, 3]
      if (hn_spatial == 0) hn_spatial <- 1e-6 
      for (jj in 1:N) {
        if (dist_spatial_sorted[jj, 4] <= H) {
          w1_matrix[jj, 1] <- (1 - (dist_spatial_sorted[jj, 3] / hn_spatial)^2)^2
        } else {
          w1_matrix[jj, 1] <- 0
        }
        w1_matrix[jj, 2] <- dist_spatial_sorted[jj, 2] 
      }
      w1_matrix <- w1_matrix[order(w1_matrix[, 2]), ]
      w1 <- w1_matrix[, 1]
      
      # --- w2 ---
      sim_dists <- sqrt(rowSums(sweep(x, 2, x[i, ], "-")^2))
      dist_sim <- cbind(sequ, sim_dists)
      dist_sim_sorted <- dist_sim[order(dist_sim[, 2]), ]
      dist_sim_sorted <- cbind(dist_sim_sorted, 1:nrow(dist_sim_sorted))
      w2_matrix <- matrix(0, N, 2)
      hn_sim <- dist_sim_sorted[H, 2]
      if (hn_sim == 0) hn_sim <- 1e-6 
      for (jj in 1:N) {
        if (dist_sim_sorted[jj, 3] <= H) {
          w2_matrix[jj, 1] <- (1 - (dist_sim_sorted[jj, 2] / hn_sim)^2)^2
        } else {
          w2_matrix[jj, 1] <- 0
        }
        w2_matrix[jj, 2] <- dist_sim_sorted[jj, 1]
      }
      w2_matrix <- w2_matrix[order(w2_matrix[, 2]), ]
      w2 <- w2_matrix[, 1]
      w <- alpha_val * w1 + beta_val * w2
    }
    # Verbose weight matrix printing (silent if verbose=FALSE)
    if (verbose && i <= 10) { 
      cat("\n[w-check] i =", i, "H =", H, "\n")
      print(w)
      cat("---------------------------------------------------\n\n")
    }
    ## MODEL SELECTION ##
    if (model=="gaussian"){
      if (det(t(x)%*%(w*x*wt))==0){
        b <- rep(0, nvar)
      }
      else{
        b <- MASS::ginv(t(x)%*%(w*x*wt))%*%t(x)%*%(w*y*wt)
      }
      uj <- x%*%b
      # update sm matrix
      if (nvar==nvarg){
        if (det(t(x)%*%(w*x*wt))==0){
          sm[i, ] <- rep(0, N)
        } else {
          sm[i, ] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt))
        }
      }
    }
    else if (model=="poisson" | model=="negbin"){
      uj <- rep(mean(y), N)
      par <- parg
      nj <- log(uj)
      ddpar <- 1
      cont <- 1
      while (abs(ddpar)>0.000001 & cont<100){
        dpar <- 1
        parold <- par
        cont1 <- 1
        cont3 <- 1
        if (model=="poisson"){
          alpha <- E^-6
          par <- 1/alpha
        }
        else{
          if (par<=E^-5 & i>1){
            par=1/alphai[i-1,2]
          }
          while (abs(dpar)>0.000001 & cont1<200){
            par <- ifelse(par<E^-10, E^-10, par)
            g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
            hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
            hess <- ifelse(hess==0, E^-23, hess)
            par0 <- par
            par <- par0-MASS::ginv(hess)*g
            if (cont1>50 & par>E^5){
              dpar <- 0.0001
              cont3 <- cont3+1
              if (cont3==1){
                par <- 2
              }
              else if (cont3==2){
                par <- E^5
              }
              else if (cont3==3){
                par <- 0.0001
              }
            }
            else{
              dpar <- par-par0
            }
            cont1 <- cont1+1
            if (par>E^6){
              par <- E^6
              dpar <- 0
            }
            if (par<=E^-5){
              par <- E^-3
              dpar <- 0
            }
          }
          alpha <- as.vector(1/par)
        }
        # Initialize loop counters and variables for tracking convergence and deviance, ensuring they start with safe finite values
        cont2  <- 0L            # iteration counter for inner loop
        olddev <- 0             # previous deviance value, start at 0 to avoid Inf-Inf issues
        dev    <- 0             # current deviance value
        ddev   <- Inf           # change in deviance, initialized large
        
        # Small numerical constant to prevent division by zero and log(0) issues
        .eps <- 1e-10
        # Ensure cont2 is a finite scalar before entering the loop
        if (!is.finite(cont2)) cont2 <- 0L
        
        # Iterative loop with NA-safe condition to avoid unexpected logical NA errors
        while (isTRUE(is.finite(ddev)) && abs(ddev) > 1e-6 && cont2 < 100L) {
          # Clamp uj to avoid exponential overflow in later steps
          uj <- pmin(uj, exp(100))
          
          # Compute ai using your original formulation, then ensure all values are positive and finite
          ai <- as.vector((uj / (1 + alpha * uj)) + (y - uj) * (alpha * uj / (1 + 2 * alpha * uj + alpha^2 * uj * uj)))
          ai <- pmax(ai, exp(-5))
          
          # Update zj vector, incorporating the offsets and current estimates
          zj <- nj + (y - uj) / (ai * (1 + alpha * uj)) - yhat_beta + fi
          
          # Weighted cross-product matrix for solving b; safe determinant calculation
          XtWA <- t(x) %*% (w * ai * x * wt)
          detval <- tryCatch(det(XtWA), error = function(e) NA_real_)
          
          # If the determinant is zero or NA, set coefficients to zero to avoid singular matrix inversion
          if (!is.finite(detval) || detval == 0) {
            b <- rep(0, nvar)
          } else {
            b <- MASS::ginv(XtWA) %*% t(x) %*% (w * ai * wt * zj)
          }
          
          # Update the linear predictor nj, clamping extreme values to avoid overflow
          nj <- x %*% b + yhat_beta - fi
          nj <- ifelse(nj > exp(2), exp(2), nj)
          
          # Update uj from nj, and clamp extremely small values to prevent underflow
          uj <- as.vector(exp(nj))
          uj <- pmax(uj, exp(-150))
          
          # Store previous deviance, then compute the new deviance depending on the model type
          olddev <- dev
          tt <- y / uj
          tt[!is.finite(tt)] <- .eps
          
          if (model == "poisson") {
            dev <- 2 * sum(y * log(pmax(tt, .eps)) - (y - uj))
          } else {
            # Negative binomial branch with safe ratio calculation
            ratio <- (1 + alpha * y) / (1 + alpha * uj)
            dev <- 2 * sum(y * log(pmax(tt, .eps)) - (y + 1 / alpha) * log(pmax(ratio, .eps)))
          }
          
          # Check for non-finite deviance; break loop early if it occurs
          if (!is.finite(dev)) {
            warning("Deviance became non-finite; breaking optimization step.")
            break
          }
          
          # Update change in deviance and increment loop counter
          ddev  <- dev - olddev
          cont2 <- cont2 + 1L
        }
        
        # Maintain your original outer loop counters and parameter differences
        cont  <- cont + 1
        ddpar <- par - parold
      }
      yhatm[i] <- uj[i]
      alphai[i, 2] <- alpha
      # update sm matrix
      if (nvar==nvarg){
        if (det(t(x)%*%(w*ai*x*wt))==0){
          sm[i, ] <- rep(0, N)
          sm3[i,] <- rep(0, nvar) 
        } else {
          C <- MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
          sm[i, ] <- (x[i,] %*% C)
          sm3[i, ] <- t(diag(C %*% diag(1/ai) %*% t(C)))
        }
      }
    }
    m1 <- (i-1)*nvar+1
    m2 <- m1+(nvar-1)
    bim[m1:m2] <- b
    yhatm[i] <- uj[i]
  }
  beta <- matrix(bim, N, byrow=T)
  yhbeta <- cbind(yhatm, beta)
  return(list(
    yhbeta = yhbeta,
    sm = sm,
    alphai = alphai,
    sm3 = sm3
  ))
}
# Golden Section Search for optimal alpha


# filepath: R/mgwr_and_alpha.R
# ----------------------------------------------------------------------------
# This file defines:
#   1) mgwr_local (Gaussian MGWR via backfitting)
#   2) find_optimal_alpha_gss (Golden-Section Search for alpha)
#
# Order: mgwr_local first, then alpha optimizer — so you can source this file
# after your existing `gwr_local` and call both.
# ----------------------------------------------------------------------------

# ------------------------------
# 1) mgwr_local (Gaussian only)
# ------------------------------
mgwr_local <- function(
    H,
    y,
    x,
    fi = NULL,
    alpha_val = 0.5,
    method = c("adaptive_bsq", "fixed_g", "fixed_bsq", "adaptive_bsq_smr"),
    model = c("gaussian"),
    N = NULL,
    nvarg = NULL,
    wt = NULL,
    E = exp(1),
    COORD,
    sequ = NULL,
    distancekm = FALSE,
    Offset = NULL,
    parg = NULL,
    yhat_beta = 0,
    max_iter = 50,
    tol = 1e-5,
    verbose = FALSE
) {
  method <- match.arg(method)
  model  <- match.arg(model)
  if (model != "gaussian") stop("mgwr_local() supports Gaussian only.")
  
  y <- as.numeric(y)
  x <- as.matrix(x)
  if (!is.null(fi)) y <- y + as.numeric(fi)
  if (is.null(N)) N <- nrow(x)
  if (is.null(nvarg)) nvarg <- ncol(x)
  stopifnot(length(y) == N, ncol(x) == nvarg, nrow(x) == N)
  
  if (is.null(wt)) wt <- rep(1, N) else wt <- as.numeric(wt)
  if (length(wt) != N) wt <- rep(wt, length.out = N)
  
  COORD <- as.matrix(COORD)
  stopifnot(nrow(COORD) == N, ncol(COORD) >= 2)
  if (is.null(sequ)) sequ <- seq_len(N)
  
  if (length(H) == 1L) H <- rep(H, nvarg)
  if (length(H) != nvarg) stop("H must be scalar or length equal to ncol(x)")
  
  # ---- Distance helper ----
  dist_i_to_all <- function(i) {
    d <- sqrt(rowSums((t(t(COORD) - COORD[i, ]))^2))
    if (distancekm) d <- d * 111
    cbind(rep.int(i, N), sequ, d)
  }
  
  # ---- Kernels ----
  kernel_fixed_g   <- function(d, h) exp(-(d / h)^2)
  kernel_fixed_bsq <- function(d, h) (pmax(0, 1 - (d / h)^2))^2
  kernel_adaptive_bsq <- function(d_sorted, h_rank) {
    hn <- d_sorted[pmin(h_rank, length(d_sorted))]
    if (hn <= 0) hn <- 1e-6
    (pmax(0, 1 - (d_sorted / hn)^2))^2
  }
  
  feature_sim_dists <- function(i) {
    sqrt(rowSums((t(t(x) - x[i, ]))^2))
  }
  
  # ---- Allocate ----
  beta <- matrix(0, nrow = N, ncol = nvarg)
  colnames(beta) <- colnames(x)
  yhat <- numeric(N)
  
  sm  <- matrix(NA_real_, N, N)     # placeholders for compatibility
  sm3 <- matrix(NA_real_, N, nvarg)
  alphai <- matrix(NA_real_, N, 3)
  
  # ---- Build weights per variable ----
  W_list <- vector("list", nvarg)
  for (k in seq_len(nvarg)) {
    Wk <- matrix(0, nrow = N, ncol = N)
    for (i in seq_len(N)) {
      dmat <- dist_i_to_all(i)
      d <- dmat[, 3]
      if (method == "fixed_g") {
        w <- kernel_fixed_g(d, H[k])
      } else if (method == "fixed_bsq") {
        w <- kernel_fixed_bsq(d, H[k])
      } else if (method == "adaptive_bsq") {
        ord <- order(d)
        wtmp <- numeric(N)
        wtmp[ord] <- kernel_adaptive_bsq(d[ord], as.integer(H[k]))
        w <- wtmp
      } else if (method == "adaptive_bsq_smr") {
        # spatial component
        ord_s <- order(d)
        w1 <- numeric(N)
        w1[ord_s] <- kernel_adaptive_bsq(d[ord_s], as.integer(H[k]))
        # similarity component
        simd <- feature_sim_dists(i)
        ord_f <- order(simd)
        w2 <- numeric(N)
        w2[ord_f] <- kernel_adaptive_bsq(simd[ord_f], as.integer(H[k]))
        w <- alpha_val * w1 + (1 - alpha_val) * w2
      }
      Wk[i, ] <- w * wt
    }
    W_list[[k]] <- Wk
    if (verbose) message(sprintf("weights built for var %d/%d", k, nvarg))
  }
  
  # ---- Backfitting loop ----
  iter <- 0L
  repeat {
    iter <- iter + 1L
    beta_old <- beta
    for (k in seq_len(nvarg)) {
      xk <- x[, k]
      partial <- y - rowSums(x[, -k, drop = FALSE] * beta[, -k, drop = FALSE])
      for (i in seq_len(N)) {
        w <- W_list[[k]][i, ]
        if (all(w == 0)) {
          beta[i, k] <- 0
        } else {
          num <- sum(w * xk * partial)
          den <- sum(w * xk * xk)
          beta[i, k] <- if (den > 0) num / den else 0
        }
      }
    }
    delta <- max(abs(beta - beta_old))
    if (verbose) message(sprintf("iter %d: max|Δβ|=%.3e", iter, delta))
    if (delta < tol || iter >= max_iter) break
  }
  
  for (i in seq_len(N)) yhat[i] <- sum(x[i, ] * beta[i, ])
  
  yhbeta <- cbind(yhat, beta)
  colnames(yhbeta) <- c("yhat", paste0("b_", seq_len(nvarg)))
  
  list(
    yhbeta = yhbeta,
    beta   = beta,
    sm     = sm,
    alphai = alphai,
    sm3    = sm3,
    iters  = iter,
    converged = iter < max_iter
  )
}

# -------------------------------------------------------------
# 2) find_optimal_alpha_gss (GWR: Gaussian/Poisson/NegBin; MGWR: Gaussian)
# -------------------------------------------------------------
find_optimal_alpha_gss <- function(
    H, Y, X, finb, N, nvarg, Offset, method, model,
    wt, E, COORD, sequ, distancekm, parg, yhat_beta,
    tol = 0.01, verbose = TRUE,
    engine = c("gwr", "mgwr"),
    criterion = c("auto", "AICc", "CV")
) {
  engine <- match.arg(engine)
  criterion <- match.arg(criterion)
  
  if (engine == "mgwr" && model != "gaussian") {
    stop("engine='mgwr' supports model='gaussian' only. Use engine='gwr' for Poisson/NegBin.")
  }
  
  if (method != "adaptive_bsq_smr" && verbose) {
    message("[warn] method != 'adaptive_bsq_smr': alpha will have no effect on weights.")
  }
  
  # ---- helpers ----
  safe_trace <- function(S) {
    if (is.null(S) || length(S) == 0L || anyNA(S)) return(NA_real_)
    sum(diag(S))
  }
  rss_of <- function(y, yhat) sum((y - yhat)^2)
  
  gaussian_aicc <- function(y, yhat, S) {
    N <- length(y)
    rss <- rss_of(y, yhat)
    edf <- safe_trace(S)
    if (is.na(edf)) return(NA_real_)
    aic <- N * log(rss / N) + N * log(2 * pi) + N + 2 * (edf + 1)
    aic + 2 * (edf + 1) * (edf + 2) / pmax(N - edf - 2, 1e-6)
  }
  
  negbin_aicc <- function(Y, mu, alphai, nvarg, N, S) {
    v1 <- safe_trace(S); if (is.na(v1)) v1 <- 0
    alpha <- alphai[, 2]
    ll <- sum(Y * log(alpha * mu) - (Y + 1/alpha) * log(1 + alpha * mu) +
                lgamma(Y + 1/alpha) - lgamma(1/alpha) - lgamma(Y + 1))
    AIC <- 2 * (v1 + v1 / nvarg) - 2 * ll
    AIC + 2 * (v1 + v1 / nvarg) * (v1 + v1 / nvarg + 1) / pmax(N - (v1 + v1 / nvarg) - 1, 1e-6)
  }
  
  poisson_aicc <- function(Y, mu, nvarg, N, S) {
    v1 <- safe_trace(S); if (is.na(v1)) v1 <- 0
    ll <- sum(Y * log(pmax(mu, 1e-12)) - mu - lgamma(Y + 1))
    AIC <- 2 * (v1 + v1 / nvarg) - 2 * ll
    AIC + 2 * (v1 + v1 / nvarg) * (v1 + v1 / nvarg + 1) / pmax(N - (v1 + v1 / nvarg) - 1, 1e-6)
  }
  
  eval_alpha <- function(alpha_val) {
    if (engine == "gwr") {
      # Uses user's gwr_local (external) with same signature as provided earlier
      res <- gwr_local(H, Y, X, finb, alpha_val, method, model, N, nvarg, wt, E,
                       COORD, sequ, distancekm, Offset, parg, yhat_beta, verbose = FALSE)
      sm <- res$sm
      yhbeta <- res$yhbeta
      
      if (model == "gaussian") {
        yhat <- yhbeta[, 1]
        crit <- gaussian_aicc(Y, yhat, sm)
        if (is.na(crit)) crit <- rss_of(Y, yhat)
      } else if (model == "poisson") {
        local_beta <- yhbeta[, 2:(nvarg + 1), drop = FALSE]
        eta <- rowSums(X * local_beta) + Offset
        mu <- exp(pmin(eta, 50))
        crit <- poisson_aicc(Y, mu, nvarg, N, sm)
      } else if (model == "negbin") {
        local_beta <- yhbeta[, 2:(nvarg + 1), drop = FALSE]
        eta <- rowSums(X * local_beta) + Offset
        mu <- exp(pmin(eta, 50))
        crit <- negbin_aicc(Y, mu, res$alphai, nvarg, N, sm)
      } else {
        stop("Unsupported model for engine='gwr'.")
      }
      if (verbose) cat("alpha=", sprintf("%.4f", alpha_val), " AICc=", sprintf("%.4f", crit), "\n", sep = "")
      return(crit)
    } else {
      # engine == "mgwr" (Gaussian only)
      res <- mgwr_local(H = H, y = Y, x = X, fi = finb, alpha_val = alpha_val,
                        method = method, model = "gaussian", N = N, nvarg = nvarg,
                        wt = wt, E = E, COORD = COORD, sequ = sequ, distancekm = distancekm,
                        Offset = Offset, parg = parg, yhat_beta = yhat_beta, verbose = FALSE)
      yhat <- res$yhbeta[, 1]
      use <- if (criterion == "auto") "RSS" else criterion
      if (use == "AICc" && !all(is.na(res$sm))) {
        crit <- gaussian_aicc(Y, yhat, res$sm)
      } else if (use == "CV" && !all(is.na(res$sm))) {
        N <- length(Y); rss <- rss_of(Y, yhat); edf <- safe_trace(res$sm)
        crit <- N * rss / pmax((N - edf)^2, 1e-6)
      } else {
        crit <- rss_of(Y, yhat)
      }
      if (verbose) cat("alpha=", sprintf("%.4f", alpha_val), " crit=", sprintf("%.4f", crit), "\n", sep = "")
      return(crit)
    }
  }
  
  # ---- Golden Section Search on [0, 1] ----
  phi <- (1 + sqrt(5)) / 2
  a <- 0.0; b <- 1.0
  fa <- eval_alpha(a)
  fb <- eval_alpha(b)
  step <- 1L
  if (verbose) cat("[alpha-opt] step ", step, ": a=", sprintf("%.4f", a), " f(a)=", sprintf("%.4f", fa),
                   " | b=", sprintf("%.4f", b), " f(b)=", sprintf("%.4f", fb), "\n", sep = "")
  
  while (abs(b - a) > tol) {
    step <- step + 1L
    if (fa < fb) {
      b_new <- a + (b - a) / phi
      fb <- eval_alpha(b_new)
      if (verbose) cat("[alpha-opt] step ", step, ": a=", sprintf("%.4f", a), " f(a)=", sprintf("%.4f", fa),
                       " | b=", sprintf("%.4f", b_new), " f(b)=", sprintf("%.4f", fb), "\n", sep = "")
      b <- b_new
    } else {
      a_new <- b - (b - a) / phi
      fa <- eval_alpha(a_new)
      if (verbose) cat("[alpha-opt] step ", step, ": a=", sprintf("%.4f", a_new), " f(a)=", sprintf("%.4f", fa),
                       " | b=", sprintf("%.4f", b), " f(b)=", sprintf("%.4f", fb), "\n", sep = "")
      a <- a_new
    }
  }
  optimal_alpha <- (a + b) / 2
  if (verbose) cat("Optimal alpha:", sprintf("%.4f", optimal_alpha), "\n")
  optimal_alpha
}



# filepath: R/mgwr_main.R
# ----------------------------------------------------------------------------
# mgwr(): Main driver for (M)GWR with optional spatial–attribute blending.
# Aligns with `gwr_local`, `mgwr_local`, and `find_optimal_alpha_gss` previously defined.
# ----------------------------------------------------------------------------

#' Main MGWR/GWR fitting wrapper with optional alpha optimization
#'
#' @param Y numeric vector (response)
#' @param X numeric matrix (design)
#' @param COORD numeric matrix of coordinates (N x 2+)
#' @param model one of "gaussian", "poisson", "negbin"
#' @param method one of "fixed_g", "fixed_bsq", "adaptive_bsq", "adaptive_bsq_smr"
#' @param mgwr logical; TRUE for MGWR (variable-specific bandwidths)
#' @param h optional scalar bandwidth; if NULL, uses GSS() per context
#' @param wt weights vector length N
#' @param E numeric base (kept for compatibility)
#' @param Offset numeric vector length N (link-scale for GLMs)
#' @param sequ integer index 1..N
#' @param distancekm logical; TRUE if coords in degrees
#' @param parg starting value for NB dispersion
#' @param int max MGWR outer iterations
#' @param verbose logical
#' @return list with fields: header, output, band, optimal_alpha (when used),
#'         yhbeta (fitted y and betas), sm, alphai, misc diagnostics
mgwr <- function(
    Y, X, COORD,
    model = c("gaussian", "poisson", "negbin"),
    method = c("fixed_g", "fixed_bsq", "adaptive_bsq", "adaptive_bsq_smr"),
    mgwr = TRUE,
    h = NULL,
    wt = NULL,
    E = exp(1),
    Offset = 0,
    sequ = NULL,
    distancekm = FALSE,
    parg = NULL,
    int = 20,
    verbose = TRUE
) {
  model  <- match.arg(model)
  method <- match.arg(method)
  
  # ---- prep ----
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  N <- nrow(X); nvarg <- ncol(X)
  if (is.null(wt)) wt <- rep(1, N) else wt <- rep(as.numeric(wt), length.out = N)
  if (length(Offset) == 1L) Offset <- rep(Offset, N)
  if (is.null(sequ)) sequ <- seq_len(N)
  COORD <- as.matrix(COORD)
  
  header <- character(0)
  output <- list()
  
  # thin wrappers for external bandwidth/GWR engines expected in environment
  stopifnot(exists("GSS"), exists("gwr"), exists("gwr_local"))
  
  # ------------------------------
  # Case 1: Standard GWR (single bandwidth)
  # ------------------------------
  if (!mgwr) {
    finb <- rep(0, N)
    yhat_beta <- Offset
    hh <- if (!is.null(h)) h else GSS(Y, X, finb)
    
    header <- append(header, "General Bandwidth")
    output$general_bandwidth <- hh
    
    if (method == "adaptive_bsq_smr") {
      optimal_alpha <- find_optimal_alpha_gss(
        H = hh, Y = Y, X = X, finb = finb, N = N, nvarg = nvarg, Offset = Offset,
        method = method, model = model, wt = wt, E = E, COORD = COORD, sequ = sequ,
        distancekm = distancekm, parg = parg, yhat_beta = yhat_beta,
        engine = "gwr", verbose = verbose
      )
      res <- gwr_local(
        H = hh, y = Y, x = X, fi = finb, alpha_val = optimal_alpha,
        method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta, verbose = verbose
      )
      yhat_beta <- res$yhbeta
      sm <- res$sm
      output$optimal_alpha <- optimal_alpha
    } else {
      yhat_beta <- gwr(hh, Y, X, finb)
      sm <- NULL
    }
    
    beta <- yhat_beta[, 2:(nvarg + 1), drop = FALSE]
    Fi <- X * beta
    
    return(list(
      header = header,
      output = output,
      yhbeta = yhat_beta,
      sm = sm
    ))
  }
  
  # ------------------------------
  # Case 2: MGWR (per-variable bandwidths + optional alpha)
  # ------------------------------
  finb <- rep(0, N)
  yhat_beta <- Offset
  hh <- if (!is.null(h)) h else GSS(Y, X, finb)
  header <- append(header, "General Bandwidth")
  output$general_bandwidth <- hh
  
  # start with a global GWR to get residuals
  yhat_beta <- gwr(hh, Y, X, finb)
  error <- Y - yhat_beta[, 1]
  beta  <- yhat_beta[, 2:(nvarg + 1), drop = FALSE]
  Fi    <- X * beta
  
  mband <- rep(hh, nvarg)
  socf <- 1
  INT <- 1
  mband_socf <- rbind(c(mband, socf))
  
  while (socf > 0.001 && INT < int) {
    fi_old <- Fi
    diffi <- 0
    fi2 <- 0
    
    for (i in seq_len(nvarg)) {
      ferror <- error + Fi[, i]
      mband[i] <- if (!is.null(h)) h else GSS(ferror, as.matrix(X[, i]), finb)
      
      # Optional alpha per-variable during backfitting
      if (method == "adaptive_bsq_smr") {
        alpha_i <- find_optimal_alpha_gss(
          H = mband[i], Y = ferror, X = as.matrix(X[, i]), finb = finb,
          N = N, nvarg = 1, Offset = Offset, method = method, model = model,
          wt = wt, E = E, COORD = COORD, sequ = sequ, distancekm = distancekm,
          parg = parg, yhat_beta = yhat_beta, engine = "gwr", verbose = FALSE
        )
      } else {
        alpha_i <- NULL
      }
      
      yhat_i <- gwr(mband[i], ferror, as.matrix(X[, i]), finb, alpha_val = alpha_i)
      beta[, i] <- yhat_i[, 2]
      Fi[, i] <- X[, i] * beta[, i]
      error <- Y - rowSums(Fi)
      
      diffi <- diffi + mean((Fi[, i] - fi_old[, i])^2)
      fi2 <- fi2 + Fi[, i]
    }
    
    socf <- sqrt(diffi / sum(fi2^2))
    INT <- INT + 1
    mband_socf <- rbind(mband_socf, c(mband, socf))
    if (verbose) message(sprintf("MGWR iter %d: socf=%.4g", INT - 1L, socf))
  }
  
  mband_socf <- mband_socf[-1, , drop = FALSE]
  band <- as.data.frame(mband_socf)
  names(band) <- c(colnames(X), "socf")
  rownames(band) <- NULL
  header <- append(header, "Bandwidth")
  output$band <- band
  
  # Final alpha over selected per-variable bandwidths
  if (method == "adaptive_bsq_smr") {
    optimal_alpha <- find_optimal_alpha_gss(
      H = mband, Y = Y, X = X, finb = finb, N = N, nvarg = nvarg, Offset = Offset,
      method = method, model = model, wt = wt, E = E, COORD = COORD, sequ = sequ,
      distancekm = distancekm, parg = parg, yhat_beta = yhat_beta,
      engine = if (model == "gaussian") "mgwr" else "gwr", verbose = verbose
    )
    output$optimal_alpha <- optimal_alpha
    
    # Apply final alpha using appropriate engine
    if (model == "gaussian") {
      res_final <- mgwr_local(
        H = mband, y = Y, x = X, fi = finb, alpha_val = optimal_alpha,
        method = method, model = "gaussian", N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta, verbose = verbose
      )
    } else {
      res_final <- gwr_local(
        H = mband, y = Y, x = X, fi = finb, alpha_val = optimal_alpha,
        method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
        COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
        parg = parg, yhat_beta = yhat_beta, verbose = verbose
      )
    }
    yhbeta <- res_final$yhbeta
    sm     <- res_final$sm
    alphai <- res_final$alphai
  } else {
    # No alpha blending; finalize with variable-specific bandwidths via gwr_local
    res_final <- gwr_local(
      H = mband, y = Y, x = X, fi = finb, alpha_val = NULL,
      method = method, model = model, N = N, nvarg = nvarg, wt = wt, E = E,
      COORD = COORD, sequ = sequ, distancekm = distancekm, Offset = Offset,
      parg = parg, yhat_beta = yhat_beta, verbose = verbose
    )
    yhbeta <- res_final$yhbeta
    sm     <- res_final$sm
    alphai <- res_final$alphai
  }
  
  v1 <- sum(diag(sm))
  if (model=='gaussian'){
    yhat <- apply(Fi, 1, sum)
    res <- Y-yhat
    rsqr1 <- t(res*wt)%*%res
    ym <- t(Y*wt)%*%Y
    rsqr2 <- ym-(sum(Y*wt)^2)/sum(wt)
    rsqr <- 1-rsqr1/rsqr2
    rsqradj <- 1-((N-1)/(N-v1))*(1-rsqr)
    sigma2 <- as.vector(N*rsqr1/((N-v1)*sum(wt)))
    root_mse <- sqrt(sigma2)
  }
  for (jj in 1:nvarg){
    m1 <- (jj-1)*N+1
    m2 <- m1+(N-1)
    ENP[jj] <- sum(diag(mrj[,m1:m2]))
    if (!mgwr){
      ENP[jj] <- sum(diag(sm))
    }
    if (model=='gaussian'){
      if (!mgwr){
        stdbm[,jj] <- sqrt(sigma2*sm3[,jj])
      } else{
        stdbm[,jj] <- sqrt(diag(sigma2*Cm[,m1:m2]%*%t(Cm[,m1:m2])))
      }
    } else{ #else if (model=='poisson' | model=='negbin' | model=='logistic'){
      if (mgwr){
        stdbm[,jj] <- sqrt(diag(Cm[,m1:m2]%*%diag(1/mAi[,jj])%*%t(Cm[,m1:m2])))
      } else{
        stdbm[,jj] <- sqrt(sm3[,jj])
      }
    }
  }
  if (model=='gaussian'){
    ll <- -N*log(rsqr1/N)/2-N*log(2*acos(-1))/2-sum((Y-yhat)*(Y-yhat))/(2*(rsqr1/N))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    stats_measures <- c(sigma2, root_mse, round(rsqr, 4),
                        round(rsqradj, 4), ll, AIC, AICc)
    names(stats_measures) <- c("sigma2e", "root_mse",
                               "R_square", "Adj_R_square",
                               "full_Log_likelihood",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  } else if (model=='poisson'){
    yhat <- exp(apply(Fi, 1, sum)+Offset)
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    dev <- 2*sum(Y*log(tt)-(Y-yhat))
    ll <- sum(-yhat+Y*log(yhat)-lgamma(Y+1))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(Y*log(tt2)-(Y-mean(Y)))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-v1))*(1-pctdev)
    stats_measures <- c(dev, ll, pctdev, adjpctdev, AIC,
                        AICc)
    names(stats_measures) <- c("deviance",
                               "full_Log_likelihood",
                               "pctdev", "adjpctdev",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  } else if (model=='negbin'){
    yhat <- exp(apply(Fi, 1, sum)+Offset)
    rmse_val <- sqrt(mean((Y - yhat)^2)) 
    
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    dev <- 2*sum(Y*log(tt)-(Y+1/alphai[,2])*log((1+alphai[,2]*Y)/(1+alphai[,2]*yhat)))
    ll <- sum(Y*log(alphai[,2]*yhat)-(Y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(Y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(Y+1))
    AIC <- 2*(v1+v1/nvarg)-2*ll
    AICc <- AIC+2*(v1+v1/nvarg)*(v1+v1/nvarg+1)/(N-(v1+v1/nvarg)-1)
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum(Y*log(tt2)-(Y+1/alphai[,2])*log((1+alphai[,2]*Y)/(1+alphai[,2]*mean(Y))))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-(v1+v1/nvarg)))*(1-pctdev)
    stats_measures <- c(rmse_val, dev, ll, pctdev, adjpctdev, AIC,
                        AICc)
    names(stats_measures) <- c("RMSE", "deviance",
                               "full_Log_likelihood",
                               "pctdev", "adjpctdev",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  } else{ #else if (model=='logistic'){
    yhat <- exp(apply(Fi, 1, sum))/(1+exp(apply(Fi, 1, sum)))
    tt <- Y/yhat
    tt <- ifelse(tt==0, E^-10, tt)
    yhat2 <- ifelse(yhat==1, 0.99999, yhat)
    tt2 <- (1-Y)/(1-yhat2)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    dev <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    lyhat2 <- 1-yhat
    lyhat2 <- ifelse(lyhat2==0, E^-10, lyhat2)
    ll <- sum(Y*log(yhat)+(1-Y)*log(lyhat2))
    AIC <- 2*v1-2*ll
    AICc <- AIC+2*(v1*(v1+1)/(N-v1-1))
    tt <- Y/mean(Y)
    tt <- ifelse(tt==0, E^-10, tt)
    tt2 <- (1-Y)/(1-mean(Y))
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnull <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    pctdev <- 1-dev/devnull
    adjpctdev <- 1-((N-1)/(N-v1))*(1-pctdev)
    stats_measures <- c(dev, ll, pctdev, adjpctdev, AIC,
                        AICc)
    names(stats_measures) <- c("deviance",
                               "full_Log_likelihood",
                               "pctdev", "adjpctdev",
                               "AIC", "AICc")
    header <- append(header, "Measures")
    output <- append(output, list(stats_measures))
    names(output)[length(output)] <- "measures"
  }
  
  ENP[nvarg+1] <- sum(diag(sm))
  ENP[nvarg+2] <- sum(diag(Sm2))
  if (model=='negbin'){
    ENP <- c(ENP, (v1/nvarg))
    names(ENP) <- c('Intercept', XVAR, 'MGWR', 'GWR', 'alpha')
  } else{
    names(ENP) <- c('Intercept', XVAR, 'MGWR', 'GWR')
  }
  header <- append(header, "ENP")
  output <- append(output, list(ENP))
  names(output)[length(output)] <- "ENP"
  dff <- N-v1
  tstat <- beta/stdbm
  probt <- 2*(1-pt(abs(tstat), dff))
  malpha <- ENP
  malpha[1:nvarg] <- 0.05/ENP[1:nvarg]
  malpha[nvarg+1] <- 0.05*(nvarg/v1)
  malpha[nvarg+2] <- 0.05*(nvarg/sum(diag(Sm2)))
  if (!mgwr){
    malpha[1:nvarg] <- 0.05*(nvarg/v1)
  }
  if (model=='negbin'){
    malpha[nvarg+3] <- 0.05*(nvarg/v1)
  }
  t_critical <- abs(qt(malpha/2,dff))
  beta2 <- beta
  if (model=='negbin'){
    alpha <- alphai[,2]
    beta2 <- cbind(beta, alpha)
  }
  qntl <- apply(beta2, 2, quantile, c(0.25, 0.5, 0.75))
  IQR <- (qntl[3,]-qntl[1,])
  qntl <- rbind(round(qntl, 6), IQR=round(IQR, 6))
  descriptb <- rbind(apply(beta2, 2, mean), apply(beta2, 2, min), apply(beta2, 2, max))
  rownames(descriptb) <- c('Mean', 'Min', 'Max')
  
  final_estimates_df <- cbind(y_observed = Y, y_predicted = yhat, beta2)
  
  if (model=='negbin'){ #release 2
    colnames(final_estimates_df) <- c('y_observed', 'y_predicted', 'Intercept', XVAR, 'alpha')#release 2
  } else{ #release 2
    colnames(final_estimates_df) <- c('y_observed', 'y_predicted', 'Intercept', XVAR) #release 2
  } #release 2
  output <- append(output, list(as.data.frame(final_estimates_df))) #release 2
  names(output)[length(output)] <- "mgwr_param_estimates" #release 2
  if (model=='negbin'){
    colnames(qntl) <- c('Intercept', XVAR, 'alpha')
  } else{
    colnames(qntl) <- c('Intercept', XVAR)
  }
  header <- append(header, "Quantiles of MGWR Parameter Estimates")
  output <- append(output, list(qntl))
  names(output)[length(output)] <- "qntls_mgwr_param_estimates"
  if (model=='negbin'){
    colnames(descriptb) <- c('Intercept', XVAR, 'alpha')
  } else{
    colnames(descriptb) <- c('Intercept', XVAR)
  }
  header <- append(header, "Descriptive Statistics")
  output <- append(output, list(descriptb))
  names(output)[length(output)] <- "descript_stats_mgwr_param_estimates"
  stdbeta <- stdbm
  stdbeta2 <- stdbeta
  if (model=='negbin'){
    stdalpha <- alphai[,3]
    stdbeta2 <- cbind(stdbeta, stdalpha)
  }
  qntls <- apply(stdbeta2, 2, quantile, c(0.25, 0.5, 0.75))
  IQR <- (qntls[3,]-qntls[1,])
  qntls <- rbind(round(qntls, 6), IQR=round(IQR, 6))
  descripts <- rbind(apply(stdbeta2, 2, mean), apply(stdbeta2, 2, min), apply(stdbeta2, 2, max))
  rownames(descripts) <- c('Mean', 'Min', 'Max')
  header <- append(header, "alpha-level=0.05")
  output <- append(output, list(malpha))
  names(output)[length(output)] <- "p_values"
  t_critical <- round(t_critical, 2)
  header <- append(header, "t-Critical")
  output <- append(output, list(t_critical))
  names(output)[length(output)] <- "t_critical"
  if (model=='negbin'){ #release 2
    colnames(stdbeta2) <- c('Intercept', XVAR, 'alpha')#release 2
  } else{ #release 2
    colnames(stdbeta2) <- c('Intercept', XVAR) #release 2
  } #release 2
  output <- append(output, list(as.data.frame(stdbeta2))) #release 2
  names(output)[length(output)] <- "mgwr_se" #release 2
  if (model=='negbin'){
    colnames(qntls) <- c('Intercept', XVAR, 'alpha')
  } else{
    colnames(qntls) <- c('Intercept', XVAR)
  }
  header <- append(header, "Quantiles of MGWR Standard Errors")
  output <- append(output, list(qntls))
  names(output)[length(output)] <- "qntls_mgwr_se"
  if (model=='negbin'){
    colnames(descripts) <- c('Intercept', XVAR, 'alpha')
  } else{
    colnames(descripts) <- c('Intercept', XVAR)
  }
  header <- append(header, "Descriptive Statistics of Standard Errors")
  output <- append(output, list(descripts))
  names(output)[length(output)] <- "descripts_stats_se"
  #### global estimates ####
  if (model=='gaussian'){
    bg <- MASS::ginv(t(X)%*%(X*wt))%*%t(X)%*%(Y*wt)
    s2g <- as.vector(t((Y-X%*%bg)*wt)%*%(Y-X%*%bg)/(N-nrow(bg)))
    varg <- diag(MASS::ginv(t(X)%*%(X*wt))*s2g)
  }
  if (is.null(weight)){
    vargd <- varg
    dfg <- N-nrow(bg)
  }
  stdg <- matrix(sqrt(vargd))
  if (model=='negbin'){
    bg <- rbind(bg, alphag)
    stdg <- rbind(stdg, sealphag)
    dfg <- dfg-1
  }
  tg <- bg/stdg
  probtg <- 2*(1-pt(abs(tg), dfg))
  bg_stdg_tg_probtg <- cbind(bg, stdg, tg, probtg)
  if (model=='negbin'){
    rownames(bg_stdg_tg_probtg) <- c('Intercept', XVAR, 'alpha')
  } else{
    rownames(bg_stdg_tg_probtg) <- c('Intercept', XVAR)
  }
  colnames(bg_stdg_tg_probtg) <- c("Par. Est.", "Std Error", "t Value", "Pr > |t|")
  header <- append(header, "Global Parameter Estimates")
  output <- append(output, list(bg_stdg_tg_probtg))
  names(output)[length(output)] <- "global_param_estimates"
  header <- append(header, "NOTE: The denominator degrees of freedom for the t tests is...")
  output <- append(output, list(dfg))
  names(output)[length(output)] <- "t_test_dfs"
  if (model=='gaussian'){
    resg <- (Y-X%*%bg)
    rsqr1g <- t(resg*wt)%*%resg
    ymg <- t(Y*wt)%*%Y
    rsqr2g <- ymg-(sum(Y*wt)^2)/sum(wt)
    rsqrg <- 1-rsqr1g/rsqr2g
    rsqradjg <- 1-((N-1)/(N-nrow(bg)))%*%(1-rsqrg)
    sigma2g <- N*rsqr1g/((N-nrow(bg))*sum(wt))
    root_mseg <- sqrt(sigma2g)
    ll <- -N*log(rsqr1g/N)/2-N*log(2*acos(-1))/2-sum(resg*resg)/(2*(rsqr1g/N))
    AIC <- -2*ll+2*nrow(bg)
    AICc <- -2*ll+2*nrow(bg)*(N/(N-nrow(bg)-1))
    global_measures <- c(sigma2g, root_mseg, round(c(rsqrg, rsqradjg), 4), ll, AIC, AICc)
    names(global_measures) <- c('sigma2e', 'root_mse', "R_square", "Adj_R_square", 'full_Log_likelihood', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  } else if (model=='poisson'){
    yhatg <- exp(X%*%bg+Offset)
    ll <- sum(-yhatg+Y*log(yhatg)-lgamma(Y+1))
    AIC <- -2*ll+2*nvarg
    AICc <- -2*ll+2*nvarg*(N/(N-nvarg-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(Y*log(tt2)-(Y-mean(Y)))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    global_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
    names(global_measures) <- c('deviance', 'full_Log_likelihood', 'pctdevg',
                                'adjpctdevg', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  } else if (model=='negbin'){
    yhatg <- exp(X%*%bg[1:(nrow(bg)-1)]+Offset)
    ll <- sum(Y*log(alphag*yhatg)-(Y+1/alphag)*log(1+alphag*yhatg)+lgamma(Y+1/alphag)-lgamma(1/alphag)-lgamma(Y+1))
    AIC <- -2*ll+2*(nvarg+1)
    AICc <- -2*ll+2*(nvarg+1)*(N/(N-(nvarg+1)-1))
    tt2 <- Y/mean(Y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(Y*log(tt2)-(Y+1/alphag)*log((1+alphag*Y)/(1+alphag*mean(Y))))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    global_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
    names(global_measures) <- c('deviance', 'full_Log_likelihood', 'pctdevg',
                                'adjpctdevg', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  } else{ #else if (model=='logistic'){
    yhatg <- exp(X%*%bg)/(1+exp(X%*%bg))
    lyhat2 <- 1-yhatg
    lyhat2 <- ifelse(lyhat2==0, E^-10, lyhat2)
    ll <- sum(Y*log(yhatg)+(1-Y)*log(lyhat2))
    AIC <- -2*ll+2*nvarg
    AICc <- -2*ll+2*nvarg*(N/(N-nvarg-1))
    tt <- Y/mean(Y)
    tt <- ifelse(tt==0, E^-10, tt)
    tt2 <- (1-Y)/(1-mean(Y))
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvarg))*(1-pctdevg)
    global_measures <- c(devg, ll, pctdevg, adjpctdevg, AIC, AICc)
    names(global_measures) <- c('deviance', 'full_Log_likelihood', 'pctdevg',
                                'adjpctdevg', 'AIC', 'AICc')
    header <- append(header, "Global Measures")
    output <- append(output, list(global_measures))
    names(output)[length(output)] <- "global_measures"
  }
  bistdt <- cbind(COORD, beta, stdbm, tstat, probt)
  colname1 <- c("Intercept", XVAR)
  parameters2 <- as.data.frame(bistdt)
  names(parameters2) <- c('x', 'y', colname1, paste('std_', colname1, sep=''), paste('tstat_', colname1, sep=''), paste('probt_', colname1, sep=''))
  sig <- matrix("not significant at 90%", N, nvarg)
  for (i in 1:N){
    for (j in 1:nvarg){
      if (probt[i,j]<0.01/ENP[j]){
        sig[i,j] <- "significant at 99%"
      } else if (probt[i,j]<0.05/ENP[j]){
        sig[i,j] <- "significant at 95%"
      } else if (probt[i,j]<0.1/ENP[j]){
        sig[i,j] <- "significant at 90%"
      } else{
        sig[i,j] <- "not significant at 90%"
      }
    }
  }
  sig_parameters2 <- as.data.frame(sig)
  names(sig_parameters2) <- c(paste('sig_', colname1, sep=''))
  if (model=='negbin'){
    atstat <- alphai[,2]/alphai[,3]
    aprobtstat <- 2*(1-pnorm(abs(atstat)))
    siga <- rep("not significant at 90%", N)
    for (i in 1:N){
      if (aprobtstat[i]<0.01*(nvarg/v1)){
        siga[i] <- "significant at 99%"
      } else if (aprobtstat[i]<0.05*(nvarg/v1)){
        siga[i] <- "significant at 95%"
      } else if (aprobtstat[i]<0.1*(nvarg/v1)){
        siga[i] <- "significant at 90%"
      } else{
        siga[i] <- "not significant at 90%"
      }
    }
    alphai <- cbind(alphai, atstat, aprobtstat)
    Alpha <- as.data.frame(alphai)
    names(Alpha) <- c("id", "alpha", "std", "tstat", "probt")
    sig_alpha <- as.data.frame(siga)
    names(sig_alpha) <- "sig_alpha"
  }
  ###################################
  min_bandwidth <- as.data.frame(t(mband))
  if (!mgwr){
    names(min_bandwidth) <- 'Intercept'
  } else{
    names(min_bandwidth) <- colname1
  }
  parameters2 <- cbind(parameters2, sig_parameters2)
  if (model=='negbin'){
    Alpha <- cbind(Alpha, sig_alpha)
  }
  # i <- 1
  # for (element in output){
  #   cat(header[i], "\n")
  #   print(element)
  #   i <- i+1
  # }
  message("NOTE: The denominator degrees of freedom for the t tests is ", dfg, ".")
  invisible(output)
  
  return(list(
    header = header,
    output = output,
    yhbeta = yhbeta,
    sm = sm,
    alphai = alphai,
    mband = mband,
    converged = (nrow(band) == 0L) || band$socf[nrow(band)] <= 0.001
  ))
}
