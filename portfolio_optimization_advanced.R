# --- Package Load ---

(function(packages){
  for (pkg in packages) {
    library(pkg, character.only = TRUE)
    cat("Package loaded:", pkg, "\n")
  }
}) (c("quantmod", "ggplot2", "xts", "quadprog", "Rsolnp"))




# --- Return Series ---

get_returns <- (function() {
  symbols <- c("AAPL", "MSFT", "NVDA", "AMZN")
  
  # 1. Download the data and extract the closing prices
  prices_xts <- do.call(
    merge,
    lapply(symbols, function(sym) {
      Cl(getSymbols(sym, from = "2020-01-01", src = "yahoo", auto.assign = FALSE))
    })
  )
  colnames(prices_xts) <- symbols
  
  # 2. Calculate the log returns
  returns_xts <- diff(log(prices_xts))
  colnames(returns_xts) <- paste0(symbols, "_RET")
  
  # Remove the first NA row
  returns_xts <- na.omit(returns_xts)
})()



# --- Markowitz Optimization ---

(function(data, rf = 0) {
  
  # Expected returns and covariance
  Sigma <- cov(data)
  mu <- colMeans(data)
  
  # Risk-free rate (adjusted to scale)
  rf
  
  
  n <- ncol(Sigma)

  ## Minimum Variance Portfolio
  Dmat <- 2 * Sigma
  dvec <- rep(0, n)
  Amat <- cbind(rep(1, n), diag(n))
  bvec <- c(1, rep(0, n))
  
  sol_min <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  w_min <- pmax(sol_min$solution, 0)
  w_min <- w_min / sum(w_min)
  ret_min <- sum(w_min * mu)
  risk_min <- sqrt(crossprod(w_min, crossprod(Sigma, w_min)))
  
  ## Maximum Sharpe Portfolio
  sharpe_fn <- function(w) {
    ret <- sum(w * mu)
    risk <- sqrt(crossprod(w, crossprod(Sigma, w)))
    -(ret - rf) / risk
  }
  
  eq_fn <- function(w) sum(w) - 1
  w0 <- rep(1/n, n)
  
  sol_sharpe <- solnp(pars = w0,
                      fun = sharpe_fn,
                      eqfun = eq_fn, eqB = 0,
                      LB = rep(0, n), UB = rep(1, n))

  w_sharpe <- sol_sharpe$pars
  ret_sharpe <- sum(w_sharpe * mu)
  risk_sharpe <- sqrt(crossprod(w_sharpe, crossprod(Sigma, w_sharpe)))
  sharpe_val <- (ret_sharpe - rf) / risk_sharpe
  
  l <- length(colnames(data))
  
  ## Results
  
  results_df <- data.frame(
    Type   = c("Max Sharpe", "Min Risk"),
    Sharpe = c(round(sharpe_val, l), NA),
    Return = c(round(ret_sharpe, l), round(ret_min, l)),
    Risk   = c(round(risk_sharpe, l), round(risk_min, l))
  )
  
  weighted_df <- rbind(
    Max_Sharpe = round(w_sharpe, l),
    Min_Risk   = round(w_min, l)
  )
  colnames(weighted_df) <- colnames(data)
  weighted_df <- as.data.frame(weighted_df)
  
  ## Output
  return(list(
    results_df = results_df,
    weighted_df = weighted_df
  ))
  }) (get_returns, r = 0)





# Risk-Return Space

(function(data, rf = 0, n_points) {
  mu <- colMeans(data)
  Sigma <- cov(data)
  n <- ncol(Sigma)
  
  # Generate random portfolio weights
  W <- matrix(runif(n_points * n), ncol = n)
  W <- W / rowSums(W)
  
  # Compute returns and risks vectorized
  rets <- W %*% mu
  risks <- sqrt(rowSums((W %*% Sigma) * W))
  sharpes <- (rets - rf) / risks
  
  # Output
  df <- data.frame(Return = rets, Risk = risks, Sharpe = sharpes)
  
  # Efficient Frontier
  grid <- seq(min(df$Risk), max(df$Risk), length.out = 300)
  
  # Plot
  plot(df$Risk, df$Return,
       col = rgb(0, 0, 1, 0.3), pch = 16,
       xlab = "Risk (σ)", ylab = "Beklenen Getiri (μ)",
       main = "Risk-Getiri Uzayı (Vektörize Matris Yöntemi)")
  grid()
  
}) (get_returns, 0, 50000)




















