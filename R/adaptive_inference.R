.adaptive_pvalue <- function(
    obs_stat,
    perm_fun,
    alternative = c("greater", "less", "two.sided"),
    eps = 0.1,
    min_perm = 200,
    max_perm = 1e5,
    batch = 200,
    seed = NULL,
    verbose = FALSE
) {
  alternative <- match.arg(alternative)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Counters
  k <- 0L   # number of extreme permutations
  B <- 0L   # total permutations

  # Optional trace for diagnostics / plotting later
  trace <- list(B = integer(), p_hat = numeric(), RSE = numeric())

  repeat {

    # --- draw a batch of permutations ---
    Tb <- replicate(batch, perm_fun())

    B <- B + length(Tb)

    # --- update extreme count ---
    if (alternative == "greater") {
      k <- k + sum(Tb >= obs_stat)
    } else if (alternative == "less") {
      k <- k + sum(Tb <= obs_stat)
    } else {  # two.sided
      k <- k + sum(Tb >= abs(obs_stat))
    }

    # --- p-value estimate with +1 correction ---
    p_hat <- (k + 1) / (B + 1)

    # --- Monte Carlo precision ---
    SE  <- sqrt(p_hat * (1 - p_hat) / B)
    RSE <- SE / p_hat

    if (verbose) {
      message(
        sprintf(
          "B = %d | p̂ = %.4g | RSE = %.3f",
          B, p_hat, RSE
        )
      )
    }

    # Store trace (cheap, optional)
    trace$B     <- c(trace$B, B)
    trace$p_hat <- c(trace$p_hat, p_hat)
    trace$RSE   <- c(trace$RSE, RSE)

    # --- stopping rules ---
    if (B >= min_perm && RSE < eps) {
      break
    }

    if (B >= max_perm) {
      break
    }
  }

  list(
    p        = p_hat,
    k        = k,
    B        = B,
    RSE      = RSE,
    converged = (RSE < eps),
    trace    = trace
  )
}
