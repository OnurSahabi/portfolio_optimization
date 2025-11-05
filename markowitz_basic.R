# --- Package Load ---

(function(packages){
for (pkg in packages) {
  library(pkg, character.only = TRUE)
  cat("Package loaded:", pkg, "\n")
  }
  }) (c("quantmod", "ggplot2", "xts"))




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

(function(data, rf = 0, step = 0.01) {
  
  # Risk-free rate (adjusted to scale)
  rf
  
  # --- Expected returns and covariance ---
  mu    <- colMeans(data)
  Sigma <- cov(data)
  n     <- ncol(data)
  grid  <- seq(0, 1, by = step)
  
  
  # --- Generate all possible weight combinations dynamically ---
  W <- expand.grid(rep(list(grid), n))
  W <- W[rowSums(W) == 1, , drop = FALSE]
  W <- as.matrix(W)
  
  # --- Containers for results ---
  best_sharpe <- -Inf
  best_sharpe_ret <- best_sharpe_risk <- NA
  best_sharpe_w <- NULL
  
  min_risk <- Inf
  min_risk_ret <- min_risk_val <- NA
  min_risk_w <- NULL
  
  # --- Single loop for all combinations ---
  for (i in seq_len(nrow(W))) {
    w <- W[i, ]
    ret  <- sum(w * mu)
    risk <- sqrt(as.numeric(t(w) %*% Sigma %*% w))
    sharpe <- (ret - rf) / risk
    
    # Maximum Sharpe portfolio
    if (sharpe > best_sharpe) {
      best_sharpe <- sharpe
      best_sharpe_ret <- ret
      best_sharpe_risk <- risk
      best_sharpe_w <- w
    }
    
    # Minimum risk portfolio
    if (risk < min_risk) {
      min_risk <- risk
      min_risk_ret <- ret
      min_risk_val <- risk
      min_risk_w <- w
    }
  }
  
  l <- length(colnames(data))
  
  # --- Results ---
  results_df <- data.frame(
    Type   = c("Max Sharpe", "Min Risk"),
    Sharpe = c(best_sharpe, NA),
    Return = c(best_sharpe_ret, min_risk_ret),
    Risk   = c(best_sharpe_risk, min_risk_val)
  )
  
  weighted_df <- rbind(
    Max_Sharpe = round(best_sharpe_w, l),
    Min_Risk   = round(min_risk_w, l)
  )
  colnames(weighted_df) <- colnames(data)
  
  list(results_df = results_df, weights = as.data.frame(weighted_df))
  
})(get_returns, r = 0, step = 0.01)





# --- Graph ---

(function(data, step = 0.01) {
  
# Step size
    grid <- seq(0, 1, by = step)
    mu    <- colMeans(data)
    Sigma <- cov(data)
    n     <- length(mu)
    
    # Combination Matrix
    W <- expand.grid(rep(list(grid), n - 1))
    W$w_last <- 1 - rowSums(W)
    W <- W[W$w_last >= 0, ]
    W <- as.matrix(W)
    
    # Risk-Return Vectors
    port_rets  <- as.vector(W %*% mu)
    port_risks <- sqrt(rowSums((W %*% Sigma) * W))
    
    # Graph
    plot(port_risks, port_rets,
         col = rgb(0, 0, 1, 0.3), pch = 16,
         xlab = "Risk (σ)", 
         ylab = "Expected Return (μ)",
         main = "Portfolio Risk-Return Space")
    grid()
  })(get_returns, 0.01)