.adaptive_decision <- function(obs_stat,
                               perm_fun,
                               alternative = c("greater", "less", "two.sided"),
                               alpha = 0.05,
                               B_min = 100,
                               B_max = 20000,
                               batch_size = 50,
                               seed = NULL) {

  alternative <- match.arg(alternative)
  if (!is.null(seed)) set.seed(seed)

  B <- 0L
  k <- 0L
  trace <- list()

  repeat {

    # run next batch
    perm_stats <- replicate(batch_size, perm_fun())
    B <- B + batch_size

    # update extreme count
    k_batch <- if (alternative == "greater") {
      sum(perm_stats >= obs_stat)
    } else if (alternative == "less") {
      sum(perm_stats <= obs_stat)
    } else {
      sum(perm_stats >= abs(obs_stat))
    }

    k <- k + k_batch

    # don't evaluate CI too early
    if (B >= B_min) {

      ci <- binom::binom.confint(
        x = k,
        n = B,
        methods = "wilson"
      )

      decision <- .decision_from_ci(ci_low = ci$lower, ci_high = ci$upper)

      trace[[length(trace) + 1]] <- list(
        B = B,
        k = k,
        p_hat = (k + 1) / (B + 1),
        ci_low = ci$lower,
        ci_high = ci$upper,
        decision = decision
      )

      # stop if decision resolved
      if (decision != "borderline") break
    }

    if (B >= B_max) break
  }

  # final values
  p_hat <- (k + 1) / (B + 1)
  mc_se <- sqrt(p_hat * (1 - p_hat) / B)

  ci <- binom::binom.confint(
    x = k,
    n = B,
    methods = "wilson"
  )

  decision_alpha <- .decision_from_ci(ci_low = ci$lower, ci_high = ci$upper)

  list(
    p = p_hat,
    B = B,
    k = k,
    mc_se = mc_se,
    p_ci_low = ci$lower,
    p_ci_high = ci$upper,
    decision_alpha = decision_alpha,
    converged = (decision_alpha != "borderline"),
    trace = trace
  )
}

