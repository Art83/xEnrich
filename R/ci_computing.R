.permutation_inference <- function(k_extreme, n_perm, alpha = 0.05) {

  p_hat <- (k_extreme + 1) / (n_perm + 1)
  mc_se <- sqrt(p_hat * (1 - p_hat) / n_perm)

  ci <- binom::binom.confint(
    x = k_extreme,
    n = n_perm,
    methods = "wilson"
  )

  decision_alpha <- if (ci$upper < alpha) {
    "sig"
  } else if (ci$lower > alpha) {
    "nonsig"
  } else {
    "borderline"
  }

  list(
    k_extreme      = k_extreme,
    n_perm         = n_perm,
    p_hat          = p_hat,
    p_ci_low       = ci$lower,
    p_ci_high      = ci$upper,
    mc_se          = mc_se,
    decision_alpha = decision_alpha,
    stability_flag = decision_alpha == "borderline"
  )
}
