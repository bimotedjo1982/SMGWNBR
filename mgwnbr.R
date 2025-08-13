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
        dev <- 0
        ddev <- 1
        cont2 <- 0
        while (abs(ddev)>0.000001 & cont2<100){
          uj <- ifelse(uj>E^100, E^100, uj)
          ai <- as.vector((uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj)))
          ai <- ifelse(ai<=0, E^-5, ai)
          zj <- nj+(y-uj)/(ai*(1+alpha*uj))-yhat_beta+fi
          if (det(t(x)%*%(w*ai*x*wt))==0){
            b <- rep(0, nvar)
          }
          else{
            b <- MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
          }
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- as.vector(exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          if (model=="poisson"){
            dev <- 2*sum(y*log(tt)-(y-uj))
          }
          else{
            dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
          }
          cont2 <- cont2+1
        }
        cont <- cont+1
        ddpar <- par-parold
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

find_optimal_alpha_gss <- function(
    H, Y, X, finb, N, nvarg, Offset, method, model, wt, E, COORD, sequ, distancekm, parg, yhat_beta, tol = 0.01, verbose = TRUE
) {
  cat("Searching optimal alpha for spatial and attribute similarity weight matrix....\n")
  phi <- (1 + sqrt(5)) / 2
  a <- 0.0
  b <- 1.0
  
  # Evaluate endpoints
  fa <- {
    aicc <- eval_aicc <- function(alpha_val) {
      res <- gwr_local(H, Y, X, finb, alpha_val, method, model, N, nvarg, wt, E, COORD, sequ, distancekm, Offset, parg, yhat_beta)
      local_yhbeta <- res$yhbeta
      local_beta <- local_yhbeta[, 2:(nvarg + 1)]
      Fi <- X * local_beta
      sm <- res$sm
      alphai <- res$alphai
      v1 <- sum(diag(sm))
      yhat <- exp(apply(Fi, 1, sum) + Offset)
      ll <- sum(Y * log(alphai[,2] * yhat) - (Y + 1/alphai[,2]) * log(1 + alphai[,2] * yhat) +
                  lgamma(Y + 1/alphai[,2]) - lgamma(1/alphai[,2]) - lgamma(Y + 1))
      AIC <- 2 * (v1 + v1/nvarg) - 2 * ll
      AICc <- AIC + 2 * (v1 + v1/nvarg) * (v1 + v1/nvarg + 1) / (N - (v1 + v1/nvarg) - 1)
      if (verbose) cat("Tested alpha =", round(alpha_val, 4), "AICc =", round(AICc, 4), "\n")
      return(AICc)
    }
    aicc(a)
  }
  fb <- {
    aicc <- eval_aicc
    aicc(b)
  }
  step <- 1
  if (verbose) cat("[alpha-opt] Step", step, "a=", round(a,4), "AICc_a=", round(fa,4), "b=", round(b,4), "AICc_b=", round(fb,4), "\n")
  while (abs(b - a) > tol) {
    step <- step + 1
    if (fa < fb) {
      # Update upper bound
      b_new <- a + (b - a)/phi
      fb <- eval_aicc(b_new)
      if (verbose) cat("[alpha-opt] Step", step, "a=", round(a,4), "AICc_a=", round(fa,4), "b=", round(b_new,4), "AICc_b=", round(fb,4), "\n")
      b <- b_new
    } else {
      # Update lower bound
      a_new <- b - (b - a)/phi
      fa <- eval_aicc(a_new)
      if (verbose) cat("[alpha-opt] Step", step, "a=", round(a_new,4), "AICc_a=", round(fa,4), "b=", round(b,4), "AICc_b=", round(fb,4), "\n")
      a <- a_new
    }
  }
  optimal_alpha <- (a + b) / 2
  cat("Optimal alpha found:", round(optimal_alpha, 4), "\n")
  return(optimal_alpha)
}

mgwnbr <- function(data, formula, weight=NULL, lat, long,
                   globalmin=TRUE, method, model="negbin",
                   mgwr=TRUE, bandwidth="cv", offset=NULL,
                   distancekm=FALSE, int=50, h=NULL){
  output <- list()
  header <- c()
  yhat_beta <- NULL
  E <- 10
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf)
  mt <- attr(mf, "terms")
  XVAR <- attr(mt, "term.labels")
  Y <- model.extract(mf, "response")
  N <- length(Y)
  X <- model.matrix(mt, mf)
  wt <-rep(1, N)
  if (!is.null(weight)){
    wt <- unlist(data[, weight])
  }
  Offset <- rep(0, N)
  if (!is.null(offset)){
    Offset <- unlist(data[, offset])
  }
  nvarg <- ncol(X)
  yhat <- rep(0, N)
  bi <- matrix(0, nvarg*N, 4)
  alphai <- matrix(0, N, 3)
  s <- rep(0, N)
  mrj <- matrix(0, N, N*nvarg)
  sm <- matrix(0, N, N)
  sm3 <- matrix(0, N, nvarg)
  rj <- matrix(0, N, N)
  Cm <- matrix(0, N, N*nvarg)
  stdbm <- matrix(0, N, nvarg)
  mAi <- matrix(0, N, nvarg)
  ENP <- rep(0, nvarg+2)
  ## global estimates ##
  if (model=="poisson" | model=="negbin"){
    uj <- (Y+mean(Y))/2
    nj <- log(uj)
    parg <- sum((Y-uj)^2/uj)/(N-nvarg)
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
      dpar <- 1
      parold <- parg
      cont1 <- 1
      cont3 <- 1
      if(model=="poisson"){
        alphag <- E^-6
        parg <- 1/alphag
      }
      else{
        if (cont>1){
          parg <- 1/(sum((Y-uj)^2/uj)/(N-nvarg))
        }
        while (abs(dpar)>0.000001 & cont1<200){
          parg <- ifelse(parg<E^-10, E^-10, parg)
          g <- sum(digamma(parg+Y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+Y)/(parg+uj))
          hess <- sum(trigamma(parg+Y)-trigamma(parg)+1/parg-2/(parg+uj)+(Y+parg)/(parg+uj)^2)
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          parg <- par0-MASS::ginv(hess)%*%g
          if (cont1>50 & parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              parg <- 2
            }
            else if (cont3==2){
              parg <- E^5
            }
            else if (cont3==3){
              parg <- 0.0001
            }
          }
          else{
            dpar <- parg-par0
          }
          cont1 <- cont1+1
          if (parg>E^6){
            parg <- E^6
            dpar <- 0
          }
        }
        alphag <- as.vector(1/parg)
      }
      devg <- 0
      ddev <- 1
      cont2 <- 0
      while (abs(ddev)>0.000001 & cont2<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        ai <- as.vector((uj/(1+alphag*uj))+(Y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj)))
        ai <- ifelse(ai<=0, E^-5, ai)
        zj <- nj+(Y-uj)/(ai*(1+alphag*uj))-Offset
        if (det(t(X)%*%(ai*X))==0){
          bg <- rep(0, nvarg)
        }
        else{
          bg <- MASS::ginv(t(X)%*%(ai*X))%*%t(X)%*%(ai*zj)
        }
        nj <- X%*%bg+Offset
        nj <- ifelse(nj>E^2, E^2, nj)
        uj <- as.vector(exp(nj))
        olddev <- devg
        uj <- ifelse(uj<E^-150, E^-150, uj)
        tt <- Y/uj
        tt <- ifelse(tt==0, E^-10, tt)
        if (model=="poisson"){
          devg <- 2*sum(Y*log(tt)-(Y-uj))
        }
        else{
          devg <- 2*sum(Y*log(tt)-(Y+1/alphag)*log((1+alphag*Y)/(1+alphag*uj)))
          sealphag <- sqrt(1/abs(hess))/(parg^2)
        }
        if (cont2>100){
          ddev <- 0.0000001
        }
        else{
          ddev <- devg-olddev
        }
        cont2 <- cont2+1
      }
      ujg <- uj
      yhat <- uj
      cont <- cont+1
      ddpar <- parg-parold
    }
    varg <- diag(MASS::ginv(t(X*wt*ai)%*%X))
  }
  else if (model=="logistic"){
    uj <- (Y+mean(Y))/2
    nj <- log(uj/(1-uj))
    devg <- 0
    ddev <- 1
    cont <- 0
    while (abs(ddev)>0.000001 & cont<100){
      uj <- ifelse(uj>E^100, E^100, uj)
      ai <- as.vector(uj*(1-uj))
      ai <- ifelse(ai<=0, E^-5, ai)
      zj <- nj+(Y-uj)/ai
      if (det(t(X)%*%(wt*ai*X))==0){
        bg <- rep(0, nvarg)
      }
      else{
        bg <- MASS::ginv(t(X)%*%(wt*ai*X))%*%t(X)%*%(wt*ai*zj)
      }
      nj <- X%*%bg
      nj <- ifelse(nj>E^2, E^2, nj)
      uj <- exp(nj)/(1+exp(nj))
      olddev <- devg
      uj <- ifelse(uj<E^-150, E^-150, uj)
      tt <- Y/uj
      tt <- ifelse(tt==0, E^-10, tt)
      uj <- ifelse(uj==1, 0.99999, uj)
      tt2 <- (1-Y)/(1-uj)
      tt2 <- ifelse(tt2==0, E^-10, tt2)
      devg <- 2*sum((Y*log(tt))+(1-Y)*log(tt2))
      ddev <- devg-olddev
      cont <- cont+1
    }
    ujg <- uj
    yhat <- uj
    varg <- diag(MASS::ginv(t(X*wt*ai)%*%X))
  }
  long <- unlist(data[, long])
  lat <- unlist(data[, lat])
  COORD <- matrix(c(long, lat), ncol=2)
  sequ <- 1:N
  cv <- function(H, y, x, fi){
    nvar <- ncol(x)
    for (i in 1:N){
      for (j in 1:N){
        seqi <- rep(i, N)
        dx <- sp::spDistsN1(COORD, COORD[i,])
        distan <- cbind(seqi, sequ, dx)
        if (distancekm){
          distan[,3] <- distan[,3]*111
        }
      }
      u <- nrow(distan)
      w <- rep(0, u)
      for (jj in 1:u){
        w[jj] <- exp(-0.5*(distan[jj,3]/H)^2)
        if (bandwidth=="cv"){
          w[i] <- 0
        }
      }
      if (method=="fixed_bsq"){
        position <- which(distan[,3]>H)
        w[position] <- 0
      }
      else if (method=="adaptive_bsq"){
        distan <- distan[order(distan[, 3]), ]
        distan <- cbind(distan, 1:nrow(distan))
        w <- matrix(0, N, 2)
        hn <- distan[H,3]
        for (jj in 1:N){
          if (distan[jj,4]<=H){
            w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
          }
          else{
            w[jj,1] <- 0
          }
          w[jj,2] <- distan[jj,2]
        }
        if (bandwidth=="cv"){
          w[which(w[,2]==i)] <- 0
        }
        w <- w[order(w[, 2]), ]
        w <- w[ ,1]
      }
      else if (method == "adaptive_bsq_smr") {
        # Define mixing parameters
        alpha_val <- 1
        beta_val <- 1 - alpha_val
        
        # --- w1---
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
        
        sim_dists <- sqrt(rowSums(sweep(x, 2, x[i, ], "-")^2))
        dist_sim <- cbind(sequ, sim_dists)
        
        dist_sim_sorted <- dist_sim[order(dist_sim[, 2]), ]
        dist_sim_sorted <- cbind(dist_sim_sorted, 1:nrow(dist_sim_sorted)) # Add rank
        
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
        
        if (bandwidth == "cv") {
          w[i] <- 0
        }
      }
      
      if (model=="gaussian"){
        if (det(t(x)%*%(w*x*wt))==0){
          b <- rep(0, nvar)
        }
        else{
          b <- MASS::ginv(t(x)%*%(w*x*wt))%*%t(x)%*%(w*y*wt)
        }
        yhat_ <- get("yhat")
        yhat_[i] <- x[i, ]%*%b
        assign("yhat", yhat_, envir=parent.frame())
        yhat[i] <- x[i, ]%*%b
        if (det(t(x)%*%(w*x*wt))==0){
          s_ <- get("s")
          s_[i] <- 0
          assign("s", s_, envir=parent.frame())
          s[i] <- 0
        }
        else{
          s_ <- get("s")
          s_[i] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt))[i]
          assign("s", s_, envir=parent.frame())
          s[i] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt))[i]
        }
        next
      }
      else if (model=="poisson" | model=="negbin"){
        uj <- yhat
        par <- parg
        nj <- log(uj)
        ddpar <- 1
        cont <- 1
        cont3 <- 0
        while (abs(ddpar)>0.000001 & cont<100){
          dpar <- 1
          parold <- par
          cont1 <- 1
          if (model=="poisson"){
            alpha <- E^-6
            par <- 1/alpha
          }
          else{
            if (par<=E^-5 & i>1){
              par <- as.vector(1/alphai[i-1, 2])
            }
            while (abs(dpar)>0.000001 & cont1<200){
              par <- ifelse(par<E^-10, E^-10, par)
              g <- sum(w*wt*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
              hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
              hess <- ifelse(hess==0, E^-23, hess)
              par0 <- par
              par <- as.vector(par0-MASS::ginv(hess)%*%g)
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
            alpha <- 1/par
          }
          dev <- 0
          ddev <- 1
          cont2 <- 1
          while (abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            assign("ai", as.vector((uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj))), envir=parent.frame())
            ai <- as.vector((uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj)))
            assign("ai", ifelse(ai<=0, E^-5, ai), envir=parent.frame())
            ai <- ifelse(ai<=0, E^-5, ai)
            zj <- nj+(y-uj)/(ai*(1+alpha*uj))-yhat_beta+fi
            if (det(t(x)%*%(w*ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
            }
            nj <- x%*%b+yhat_beta-fi
            nj <- ifelse(nj>E^2, E^2, nj)
            uj <- exp(nj)
            olddev <- dev
            uj <- ifelse(uj<E^-150, E^-150, uj)
            tt <- y/uj
            tt <- ifelse(tt==0, E^-10, tt)
            if (model=="poisson"){
              dev <- 2*sum(y*log(tt)-(y-uj))
            }
            else{
              dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
            }
            if (cont2>100){
              ddev <- 0.0000001
            }
            else{
              ddev <- dev-olddev
            }
            cont2 <- cont2+1
          }
          cont <- cont+1
          ddpar <- par-parold
        }
        yhat_ <- get("yhat")
        yhat_[i] <- uj[i]
        assign("yhat", yhat_, envir=parent.frame())
        yhat[i] <- uj[i]
        alphai_ <- get("alphai")
        alphai_[i, 2] <- alpha
        assign("alphai", alphai_, envir=parent.frame())
        alphai[i, 2] <- alpha
        if (det(t(x)%*%(w*ai*x*wt))==0){
          s_ <- get("s")
          s_[i] <- 0
          assign("s", s_, envir=parent.frame())
          s[i] <- 0
        }
        else{
          s_ <- get("s")
          s_[i] <- (x[i, ]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*ai*wt))[i]
          assign("s", s_, envir=parent.frame())
          s[i] <- (x[i, ]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*ai*wt))[i]
        }
        next
      }
      else if (model=="logistic"){
        uj <- yhat
        nj <- log(uj/(1-uj))
        dev <- 0
        ddev <- 1
        cont <- 0
        while (abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100, E^100, uj)
          assign("ai", as.vector(uj*(1-uj)), envir=parent.frame())
          ai <- as.vector(uj*(1-uj))
          assign("ai", ifelse(ai<=0, E^-5, ai), envir=parent.frame())
          ai <- ifelse(ai<=0, E^-5, ai)
          zj <- nj+(y-uj)/ai-yhat_beta+fi
          if (det(t(x)%*%(w*ai*x*wt))==0){
            b <- rep(0, nvar)
          }
          else{
            b <- MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
          }
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- exp(nj)/(1+exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          uj <- ifelse(uj==1, 0.99999, uj)
          tt2 <- (1-y)/(1-uj)
          tt2 <- ifelse(tt2==0,E^-10, tt2)
          dev <- 2*sum((y*log(tt))+(1-y)*log(tt2))
          if (cont>100){
            ddev <-  0.0000001
          }
          else{
            ddev <- dev-olddev
          }
        }
        yhat_ <- get("yhat")
        yhat_[i] <- uj[i]
        assign("yhat", yhat_, envir=parent.frame())
        yhat[i] <- uj[i]
        if (det(t(x)%*%(w*ai*x*wt))==0){
          s_ <- get("s")
          s_[i] <- 0
          assign("s", s_, envir=parent.frame())
          s[i] <- 0
        }
        else{
          s_ <- get("s")
          s_[i] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))[i]
          assign("s", s_, envir=parent.frame())
          s[i] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))[i]
        }
        next
      }
    }
    if (model=="gaussian"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      npar <- sum(s)
      AICc <- 2*N*log(CV/N)+N*log(2*3.14159)+N*(N+npar)/(N-2-npar)
    }
    else if (model=="poisson" | model=="negbin"){
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      if (model=="poisson"){
        ll <- sum(-yhat+y*log(yhat)-lgamma(y+1))
        npar <- sum(s)
      }
      else{
        ll <- sum(y*log(alphai[,2]*yhat)-(y+1/alphai[,2])*log(1+alphai[,2]*yhat)+lgamma(y+1/alphai[,2])-lgamma(1/alphai[,2])-lgamma(y+1))
        npar <- sum(s)+sum(s)/nvar
      }
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(N-npar-1)
    }
    else if (model=="logistic"){
      uj <- ifelse(uj==0, E^-10, uj)
      uj <- ifelse(uj==1, 0.99999, uj)
      CV <- t((y-yhat)*wt)%*%(y-yhat)
      ll <- sum(y*log(uj)-(1-y)*log(1-uj))
      npar <- sum(s)
      AIC <- 2*npar-2*ll
      AICC <- AIC +(2*npar*(npar+1))/(N-npar-1)
    }
    if (bandwidth=="aic"){
      CV <- AICC
    }
    res <- cbind(CV, npar)
    return (res)
  }
  GSS <- function(depy, indepx, fix){
    # DEFINING GOLDEN SECTION SEARCH PARAMETERS #
    if(method=="fixed_g" | method=="fixed_bsq"){
      ax <- 0
      bx <- as.integer(max(dist(COORD))+1)
      if (distancekm){
        bx <- bx*111
      }
    }
    else if (method=="adaptive_bsq" | method=="adaptive_bsq_smr"){
      ax <- 5
      bx <- N
    }
    r <- 0.61803399
    tol <- 0.1
    if (!globalmin){
      lower <- ax
      upper <- bx
      xmin <- matrix(0, 1, 2)
      GMY <- 1
      ax1 <- lower[GMY]
      bx1 <- upper[GMY]
      h0 <- ax1
      h3 <- bx1
      h1 <- bx1-r*(bx1-ax1)
      h2 <- ax1+r*(bx1-ax1)
      res1 <- cv(h1, depy, indepx, fix)
      assign("s", s, envir=parent.frame())
      if (model!="gaussian"){ #release 2
        assign("ai", ai, envir=parent.frame())
      }#release 2
      assign("yhat", yhat, envir=parent.frame())
      assign("alphai", alphai, envir=parent.frame())
      CV1 <- res1[1]
      res2 <- cv(h2,depy,indepx,fix)
      assign("s", s, envir=parent.frame())
      if (model!="gaussian"){ #release 2
        assign("ai", ai, envir=parent.frame())
      } #release 2
      assign("yhat", yhat, envir=parent.frame())
      assign("alphai", alphai, envir=parent.frame())
      CV2 <- res2[1]
      INT <- 1
      while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & INT<200){
        if (CV2<CV1){
          h0 <- h1
          h1 <- h3-r*(h3-h0)
          h2 <- h0+r*(h3-h0)
          CV1 <- CV2
          res2 <- cv(h2,depy,indepx,fix)
          assign("s", s, envir=parent.frame())
          if (model!="gaussian"){ #release 2
            assign("ai", ai, envir=parent.frame())
          } #release 2
          assign("yhat", yhat, envir=parent.frame())
          assign("alphai", alphai, envir=parent.frame())
          CV2 <- res2[1]
        }
        else{
          h3 <- h2
          h1 <- h3-r*(h3-h0)
          h2 <- h0+r*(h3-h0)
          CV2 <- CV1
          res1 <- cv(h1, depy, indepx, fix)
          assign("s", s, envir=parent.frame())
          if (model!="gaussian"){ #release 2
            assign("ai", ai, envir=parent.frame())
          } #release 2
          assign("yhat", yhat, envir=parent.frame())
          assign("alphai", alphai, envir=parent.frame())
          CV1 <- res1[1]
        }
        INT <- INT+1
      }
      if (CV1<CV2){
        golden <- CV1
        xmin[GMY, 1] <- golden
        xmin[GMY, 2] <- h1
        npar <- res1[1]
        if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
          xmin[GMY, 2] <- floor(h1)
          xming <- xmin[GMY, 2]
        }
      }
      else{
        golden <- CV2
        xmin[GMY, 1] <- golden
        xmin[GMY, 2] <- h2
        npar <- res2[1]
        if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
          xmin[GMY, 2] <- floor(h2)
          xming <- xmin[GMY, 2]
        }
      }
      xming <- xmin[GMY, 2]
    }
    else{
      lower <- cbind(ax, (1-r)*bx, r*bx)
      upper <- cbind((1-r)*bx, r*bx, bx)
      xmin <- matrix(0, 3, 2)
      for (GMY in 1:3){
        ax1 <- lower[GMY]
        bx1 <- upper[GMY]
        h0 <- ax1
        h3 <- bx1
        h1 <- bx1-r*(bx1-ax1)
        h2 <- ax1+r*(bx1-ax1)
        res1 <- cv(h1, depy, indepx, fix)
        assign("s", s, envir=parent.frame())
        if (model!="gaussian"){ #release 2
          assign("ai", ai, envir=parent.frame())
        } #release 2
        assign("yhat", yhat, envir=parent.frame())
        assign("alphai", alphai, envir=parent.frame())
        CV1 <- res1[1]
        res2 <- cv(h2,depy,indepx,fix)
        assign("s", s, envir=parent.frame())
        if (model!="gaussian"){ #release 2
          assign("ai", ai, envir=parent.frame())
        } #release 2
        assign("yhat", yhat, envir=parent.frame())
        assign("alphai", alphai, envir=parent.frame())
        CV2 <- res2[1]
        INT <- 1
        while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & INT<200){
          if (CV2<CV1){
            h0 <- h1
            h1 <- h3-r*(h3-h0)
            h2 <- h0+r*(h3-h0)
            CV1 <- CV2
            res2 <- cv(h2,depy,indepx,fix)
            assign("s", s, envir=parent.frame())
            if (model!="gaussian"){ #release 2
              assign("ai", ai, envir=parent.frame())
            } #release 2
            assign("yhat", yhat, envir=parent.frame())
            assign("alphai", alphai, envir=parent.frame())
            CV2 <- res2[1]
          }
          else{
            h3 <- h2
            h1 <- h3-r*(h3-h0)
            h2 <- h0+r*(h3-h0)
            CV2 <- CV1
            res1 <- cv(h1, depy, indepx, fix)
            assign("s", s, envir=parent.frame())
            if (model!="gaussian"){ #release 2
              assign("ai", ai, envir=parent.frame())
            } #release 2
            assign("yhat", yhat, envir=parent.frame())
            assign("alphai", alphai, envir=parent.frame())
            CV1 <- res1[1]
          }
          INT <- INT+1
        }
        if (CV1<CV2){
          golden <- CV1
          xmin[GMY,1] <- golden
          xmin[GMY,2] <- h1
          npar <- res1[1]
          if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
            xmin[GMY,2] <- floor(h1)
            xming <- xmin[GMY,2]
          }
        }
        else{
          golden <- CV2
          xmin[GMY,1] <- golden
          xmin[GMY,2] <- h2
          npar <- res2[1]
          if (method=="adaptive_bsq"| method=="adaptive_bsq_smr"){
            xmin[GMY,2] <- floor(h2)
            xming <- xmin[GMY,2]
          }
        }
        xming <- xmin[GMY,2]
      }
    }
    if (globalmin){
      xming <- xmin[which(xmin[,1]==min(xmin[,1])),2]
    }
    bandwidth <- xming
    return (bandwidth)
  }
  
  gwr <- function(H, y, x, fi, alpha_val = NULL) {
    # alpha_val dapat dipassing dari luar, default ke 0.5 jika NULL
    nvar <- ncol(x)
    bim <- rep(0, nvar*N)
    yhatm <- rep(0, N)
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
        # PATCH: gunakan alpha_val dari parameter, default 0.5 jika NULL
        if (is.null(alpha_val)) alpha_val <- 0.5
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
        dist_sim_sorted <- cbind(dist_sim_sorted, 1:nrow(dist_sim_sorted)) # Add rank
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
      ## MODEL SELECTION ##
      if (model=="gaussian"){
        if (det(t(x)%*%(w*x*wt))==0){
          b <- rep(0, nvar)
        }
        else{
          b <- MASS::ginv(t(x)%*%(w*x*wt))%*%t(x)%*%(w*y*wt)
        }
        uj <- x%*%b
        if (nvar==nvarg){
          if (det(t(x)%*%(w*x*wt))==0){
            sm_ <- get("sm")
            sm_[i, ] <- rep(0, N)
            assign("sm", sm_, envir=parent.frame())
            mrj_ <- get("mrj")
            mrj_[i, ] <- matrix(0, N*nvar)
            assign("mrj", mrj_, envir=parent.frame())
          }
          else{
            ej <- diag(nvar)
            sm_ <- get("sm")
            sm_[i, ] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt))
            assign("sm", sm_, envir=parent.frame())
            sm3_ <- get("sm3")
            sm3_[i, ] <- t(diag((MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt))%*%t(MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt))))
            assign("sm3", sm3_, envir=parent.frame())
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj_ <- get("mrj")
              mrj_[i, m1:m2] <- (x[i,jj]*ej[jj,])%*%MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt)
              assign("mrj", mrj_, envir=parent.frame())
            }
          }
        }
        else{
          if (det(t(x)%*%(w*x*wt))==0){
            rj_ <- get("rj")
            rj_[i, ] <- rep(0, N)
            assign("rj", rj_, envir=parent.frame())
          }
          else{
            rj_ <- get("rj")
            rj_[i, ] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*x*wt))%*%t(x*w*wt))
            assign("rj", rj_, envir=parent.frame())
          }
        }
      }
      else if (model=="poisson" | model=="negbin"){
        uj <- yhat
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
          else{ #model=='NEGBIN'
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
          dev <- 0
          ddev <- 1
          cont2 <- 0
          while (abs(ddev)>0.000001 & cont2<100){
            uj <- ifelse(uj>E^100, E^100, uj)
            assign("ai", as.vector((uj/(1+alpha*uj))+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj*uj))), envir=parent.frame())
            assign("ai", ifelse(ai<=0, E^-5, ai), envir=parent.frame())
            zj <- nj+(y-uj)/(ai*(1+alpha*uj))-yhat_beta+fi
            if (det(t(x)%*%(w*ai*x*wt))==0){
              b <- rep(0, nvar)
            }
            else{
              b <- MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*wt*zj)
            }
            nj <- x%*%b+yhat_beta-fi
            nj <- ifelse(nj>E^2, E^2, nj)
            uj <- as.vector(exp(nj))
            olddev <- dev
            uj <- ifelse(uj<E^-150, E^-150, uj)
            tt <- y/uj
            tt <- ifelse(tt==0, E^-10, tt)
            if (model=="poisson"){
              dev <- 2*sum(y*log(tt)-(y-uj))
            }
            else{ #model=="NEGBIN"
              dev <- 2*sum(y*log(tt)-(y+1/alpha)*log((1+alpha*y)/(1+alpha*uj)))
            }
            cont2 <- cont2+1
          }
          cont <- cont+1
          ddpar <- par-parold
        }
        if (nvar==nvarg){
          if (det(t(x)%*%(w*ai*x*wt))==0){
            sm_ <- get("sm")
            sm_[i, ] <- c(0, N)
            assign("sm", sm_, envir=parent.frame())
            mrj_ <- get("mrj")
            mrj_[i, ] <- rep(0, N*nvar)
            assign("mrj", mrj_, envir=parent.frame())
          }
          else{
            ej <- diag(nvar)
            sm_ <- get("sm")
            sm_[i, ] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
            assign("sm", sm_, envir=parent.frame())
            sm3_ <- get("sm3")
            sm3_[i, ] <- t(diag((MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))%*%diag(1/ai)%*%t(MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))))
            assign("sm3", sm3_, envir=parent.frame())
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj_ <- get("mrj")
              mrj_[i, m1:m2] <- (x[i,jj]*ej[jj,])%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
              assign("mrj", mrj_, envir=parent.frame())
            }
          }
        }
        else{
          if (det(t(x)%*%(w*ai*x*wt))==0){
            rj_ <- get("rj")
            rj_[i, ] <- rep(0, N)
            assign("rj", rj_, envir=parent.frame())
          }
          else{
            rj_ <- get("rj")
            rj_[i, ] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
            assign("rj", rj_, envir=parent.frame())
          }
        }
        if (model=="negbin"){
          hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+exp(yhat_beta))+(y+par)/(par+exp(yhat_beta))^2))
          if (!mgwr){
            hess <- sum(w*wt*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
            hess <- ifelse(hess==0, E^-23, hess)
          }
          sealpha <- sqrt(1/abs(hess))/(par^2)
          alphai_ <- get("alphai")
          alphai_[i, 1] <- i
          assign("alphai", alphai_, envir=parent.frame())
          alphai_ <- get("alphai")
          alphai_[i, 2] <- alpha
          assign("alphai", alphai_, envir=parent.frame())
          alphai_ <- get("alphai")
          alphai_[i, 3] <- sealpha
          assign("alphai", alphai_, envir=parent.frame())
        }
      }
      else{ #else if (model=="logistic"){
        uj <- yhat
        nj <- log(uj/(1-uj))
        dev <- 0
        ddev <- 1
        cont <- 1
        while (abs(ddev)>0.000001 & cont<100){
          cont <- cont+1
          uj <- ifelse(uj>E^100, E^100, uj)
          assign("ai", as.vector(uj*(1-uj)), envir=parent.frame())
          assign("ai", ifelse(ai<=0, E^-5, ai), envir=parent.frame())
          zj <- nj+(y-uj)/ai-yhat_beta+fi
          if (det(t(x)%*%(w*ai*x*wt))==0){
            b <- rep(0, nvar)
          }
          else{
            b <- MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x)%*%(w*ai*zj*wt)
          }
          nj <- x%*%b+yhat_beta-fi
          nj <- ifelse(nj>E^2, E^2, nj)
          uj <- exp(nj)/(1+exp(nj))
          olddev <- dev
          uj <- ifelse(uj<E^-150, E^-150, uj)
          tt <- y/uj
          tt <- ifelse(tt==0, E^-10, tt)
          uj <- ifelse(uj==1, 0.99999, uj)
          tt2 <- (1-y)/(1-uj)
          tt2 <- ifelse(tt2==0, E^-10, tt2)
          dev <- 2*sum((y*log(tt))+(1-y)*log(tt2))
          if (cont>100){
            ddev <- 0.0000001
          }
          else{
            ddev <- dev-olddev
          }
        }
        if (nvar==nvarg){
          if (det(t(x)%*%(w*ai*x*wt))==0){
            sm_ <- get("sm")
            sm_[i, ] <- rep(0, N)
            assign("sm", sm_, envir=parent.frame())
            mrj_ <- get("mrj")
            mrj_[i, ] <- matrix(0, N*nvar)
            assign("mrj", mrj_, envir=parent.frame())
          }
          else{
            ej <- diag(nvar)
            sm_ <- get("sm")
            sm_[i, ] <- x[i,]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
            assign("sm", sm_, envir=parent.frame())
            sm3_ <- get("sm3")
            sm3_[i, ] <- t(diag((MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))%*%diag(1/ai)%*%t(MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))))
            assign("sm3", sm3_, envir=parent.frame())
            for (jj in 1:nvar){
              m1 <- (jj-1)*N+1
              m2 <- m1+(N-1)
              mrj_ <- get("mrj")
              mrj_[i, m1:m2] <- (x[i,jj]*ej[jj,])%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai)
              assign("mrj", mrj_, envir=parent.frame())
            }
          }
        }
        else{
          if (det(t(x)%*%(w*ai*x*wt))==0){
            rj_ <- get("rj")
            rj_[i, ] <- rep(0, N)
            assign("rj", rj_, envir=parent.frame())
          }
          else{
            rj_ <- get("rj")
            rj_[i, ] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*ai*x*wt))%*%t(x*w*wt*ai))
            assign("rj", rj_, envir=parent.frame())
          }
        }
      }
      m1 <- (i-1)*nvar+1
      m2 <- m1+(nvar-1)
      bim[m1:m2] <- b
      yhatm[i] <- uj[i]
      yhat_ <- get("yhat")
      yhat_[i] <- uj[i]
      assign("yhat", yhat_, envir=parent.frame())
    }
    beta <- matrix(bim, N, byrow=T)
    yhbeta <- cbind(yhatm, beta)
    return (yhbeta)
  }
  
  if (!mgwr) {
  finb <- rep(0, N)
  yhat_beta <- Offset
  if (!is.null(h)) {
    hh <- h
  } else {
    hh <- GSS(Y, X, finb)
  }
  header <- append(header, "General Bandwidth")
  output <- append(output, hh)
  names(output) <- "general_bandwidth"
  
  yhat_beta <- gwr(hh, Y, X, finb)
  error <- Y - yhat_beta[, 1]
  beta <- yhat_beta[, 2:(nvarg + 1)]
  Fi <- X * beta
  
  mband <- rep(hh, nvarg)
  socf <- 1
  INT <- 1
  mband_socf <- c(mband, socf)
  
  while (socf > 0.001 && INT < int) {
    fi_old <- Fi
    diffi <- 0
    fi2 <- 0
    
    for (i in 1:nvarg) {
      ferror <- error + Fi[, i]
      if (!is.null(h)) {
        mband[i] <- h
      } else {
        mband[i] <- GSS(ferror, as.matrix(X[, i]), finb)
      }
      
      # Re-evaluate alpha for each bandwidth using GSS
      if (method == "adaptive_bsq_smr") {
        optimal_alpha <- find_optimal_alpha_gss(
          mband[i], ferror, as.matrix(X[, i]), finb, N, 1, Offset,
          method, model, wt, E, COORD, sequ, distancekm, parg, yhat_beta
        )
      } else {
        optimal_alpha <- NULL  # Use default or other method
      }
      
      yhat_beta <- gwr(mband[i], ferror, as.matrix(X[, i]), finb, alpha_val = optimal_alpha)
      beta[, i] <- yhat_beta[, 2]
      Fi[, i] <- X[, i] * beta[, i]
    }
    
    diffi <- sum((Fi - fi_old)^2)
    fi2 <- sum(Fi^2)
    socf <- sqrt(diffi / fi2)
    INT <- INT + 1
  }
  
  if (INT >= int) {
    band <- as.data.frame(mband_socf)
  } else {
    band <- as.data.frame(mband_socf)
  }
  names(band) <- c("Intercept", XVAR, "socf")
  rownames(band) <- NULL
  header <- append(header, "Bandwidth")
  output <- append(output, list(band))
  names(output)[length(output)] <- "band"
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
      }
      else{
        stdbm[,jj] <- sqrt(diag(sigma2*Cm[,m1:m2]%*%t(Cm[,m1:m2])))
      }
    }
    else{ #else if (model=='poisson' | model=='negbin' | model=='logistic'){
      if (mgwr){
        stdbm[,jj] <- sqrt(diag(Cm[,m1:m2]%*%diag(1/mAi[,jj])%*%t(Cm[,m1:m2])))
      }
      else{
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
  }
  else if (model=='poisson'){
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
  }
  else if (model=='negbin'){
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
  }
  else{ #else if (model=='logistic'){
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
  }
  else{
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
  } #release 2
  else{ #release 2
    colnames(final_estimates_df) <- c('y_observed', 'y_predicted', 'Intercept', XVAR) #release 2
  } #release 2
  output <- append(output, list(as.data.frame(final_estimates_df))) #release 2
  names(output)[length(output)] <- "mgwr_param_estimates" #release 2
  if (model=='negbin'){
    colnames(qntl) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    colnames(qntl) <- c('Intercept', XVAR)
  }
  header <- append(header, "Quantiles of MGWR Parameter Estimates")
  output <- append(output, list(qntl))
  names(output)[length(output)] <- "qntls_mgwr_param_estimates"
  if (model=='negbin'){
    colnames(descriptb) <- c('Intercept', XVAR, 'alpha')
  }
  else{
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
  } #release 2
  else{ #release 2
    colnames(stdbeta2) <- c('Intercept', XVAR) #release 2
  } #release 2
  output <- append(output, list(as.data.frame(stdbeta2))) #release 2
  names(output)[length(output)] <- "mgwr_se" #release 2
  if (model=='negbin'){
    colnames(qntls) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    colnames(qntls) <- c('Intercept', XVAR)
  }
  header <- append(header, "Quantiles of MGWR Standard Errors")
  output <- append(output, list(qntls))
  names(output)[length(output)] <- "qntls_mgwr_se"
  if (model=='negbin'){
    colnames(descripts) <- c('Intercept', XVAR, 'alpha')
  }
  else{
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
  }
  else{
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
  }
  else if (model=='poisson'){
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
  }
  else if (model=='negbin'){
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
  }
  else{ #else if (model=='logistic'){
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
      }
      else if (probt[i,j]<0.05/ENP[j]){
        sig[i,j] <- "significant at 95%"
      }
      else if (probt[i,j]<0.1/ENP[j]){
        sig[i,j] <- "significant at 90%"
      }
      else{
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
      }
      else if (aprobtstat[i]<0.05*(nvarg/v1)){
        siga[i] <- "significant at 95%"
      }
      else if (aprobtstat[i]<0.1*(nvarg/v1)){
        siga[i] <- "significant at 90%"
      }
      else{
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
  }
  else{
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
}

