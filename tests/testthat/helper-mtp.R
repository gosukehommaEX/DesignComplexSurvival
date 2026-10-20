# Mean vector and covariance matrix of the statistics of two endpoints whose
# score processes have independent bivariate normal increments between the
# analyses, with correlation rho (the canonical joint distribution); the
# same definition as canonical() in inst/validation/python/ref_mtp.py.
# tau1, tau2: information fractions at all analyses; a1, a2: analyses at
# which each endpoint is tested.
mtp_canonical <- function(tau1, tau2, rho, theta1, theta2, a1, a2) {
  cinc <- rho * sqrt(diff(c(0, tau1)) * diff(c(0, tau2)))
  t1 <- tau1[a1]
  t2 <- tau2[a2]
  s11 <- sqrt(outer(t1, t1, pmin) / outer(t1, t1, pmax))
  s22 <- sqrt(outer(t2, t2, pmin) / outer(t2, t2, pmax))
  s12 <- matrix(0, length(a1), length(a2))
  for (p in seq_along(a1)) {
    for (q in seq_along(a2)) {
      s12[p, q] <- sum(cinc[seq_len(min(a1[p], a2[q]))]) / sqrt(t1[p] * t2[q])
    }
  }
  list(mean = c(theta1 * sqrt(t1), theta2 * sqrt(t2)),
       sigma = rbind(cbind(s11, s12), cbind(t(s12), s22)))
}

# Numeric vector from a ";"-separated field of the expected-value files
# (NULL for an empty field)
mtp_num <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(NULL)
  }
  as.numeric(strsplit(x, ";", fixed = TRUE)[[1L]])
}

# Rejection probabilities of gs_prob_mtp() named as in gs_prob_mtp.csv
mtp_quantities <- function(res) {
  lk <- res$looks
  e1 <- lk[lk$endpoint == 1, ]
  e2 <- lk[lk$endpoint == 2, ]
  c(res$reject,
    stats::setNames(e1$prob_upper, paste0("p1_", e1$look)),
    stats::setNames(e1$prob_lower, paste0("pf_", e1$look)),
    stats::setNames(e2$prob_upper, paste0("p2_", e2$look)))
}
