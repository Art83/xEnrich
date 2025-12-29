.ite_overlap_bits <- function(N, n, m, k) {
  # N: universe size
  # n: |list|
  # m: |pathway|
  # k: overlap

  # probabilities
  p11 <- k / N
  p10 <- (n - k) / N
  p01 <- (m - k) / N
  p00 <- (N - n - m + k) / N

  p1. <- n / N
  p0. <- 1 - p1.
  p.1 <- m / N
  p.0 <- 1 - p.1

  term <- function(p, pa, pb) {
    if (p <= 0) return(0)
    p * log2(p / (pa * pb))
  }

  sum(
    term(p11, p1., p.1),
    term(p10, p1., p.0),
    term(p01, p0., p.1),
    term(p00, p0., p.0)
  )
}
