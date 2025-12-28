.decision_from_ci <- function(ci_low, ci_high, alpha = 0.05, tol = 1e-4) {
  if (ci_high <= (alpha - tol)) return("sig")
  if (ci_low  >= (alpha + tol)) return("nonsig")
  "borderline"
}
