# =====================================================================
#  Mortality modelling helper functions
#  Model: deaths D_x ~ Poisson(E_x * mu_x), with E_x = central exposure
# =====================================================================

# ---- Force of mortality laws ----------------------------------------
# Parameters are on the log scale so optim() can search freely
# (keeps A, B > 0 automatically).

mu_gompertz <- function(par, x) exp(par[1] + par[2] * x)                    # mu = B c^x
mu_makeham  <- function(par, x) exp(par[1]) + exp(par[2] + par[3] * x)      # mu = A + B c^x

# ---- Poisson negative log-likelihood ---------------------------------
# log L = sum( D log(E mu) - E mu - log(D!) ); the last term has no
# parameters, so it is dropped for optimisation and added back for AIC.
poisson_nll <- function(par, mu_fun, x, D, E) {
  mu <- mu_fun(par, x)
  -sum(D * log(E * mu) - E * mu)
}

# ---- Gradient of the negative log-likelihood (the score, with a minus) --
# d(-log L)/d(theta) = -sum( (D / mu - E) * d(mu)/d(theta) )
# Supplying it lets BFGS converge exactly; the parameters are highly
# correlated, and numerical gradients alone stop short of the optimum.
dmu_gompertz <- function(par, x) { g <- exp(par[1] + par[2] * x); cbind(g, g * x) }
dmu_makeham  <- function(par, x) { g <- exp(par[2] + par[3] * x); cbind(exp(par[1]), g, g * x) }

poisson_grad <- function(par, mu_fun, x, D, E) {
  mu   <- mu_fun(par, x)
  dmu  <- if (length(par) == 2) dmu_gompertz(par, x) else dmu_makeham(par, x)
  -colSums((D / mu - E) * dmu)
}

# ---- Fit a law by maximum likelihood ----------------------------------
fit_law <- function(mu_fun, start, x, D, E) {
  opt <- optim(start, poisson_nll, gr = poisson_grad, mu_fun = mu_fun, x = x, D = D, E = E,
               method = "BFGS", hessian = TRUE,
               control = list(maxit = 10000, reltol = 1e-15))
  if (opt$convergence != 0) warning("optim did not converge")
  loglik <- -opt$value - sum(lfactorial(D))
  list(par = opt$par,
       se  = sqrt(diag(solve(opt$hessian))),   # inverse observed information
       loglik = loglik,
       aic = -2 * loglik + 2 * length(start),
       mu_fun = mu_fun)
}

# ---- Integrated hazard over [x, x+1] ----------------------------------
# p_x = exp(-integral of mu from x to x+1). Exact for these laws:
#   integral of B c^t dt from x to x+1 = B c^x (c - 1) / ln c
H_makeham <- function(par, x) {
  A <- if (length(par) == 3) exp(par[1]) else 0
  a <- par[length(par) - 1]; b <- par[length(par)]
  A + exp(a + b * x) * (exp(b) - 1) / b
}

# ---- Life table --------------------------------------------------------
# ages: integer ages from start age up to limiting age omega (q_omega = 1)
life_table <- function(par, ages, radix = 100000, mort_scale = 1) {
  q <- 1 - exp(-mort_scale * H_makeham(par, ages))
  q[length(q)] <- 1                               # close the table
  l <- radix * c(1, cumprod(1 - q))[seq_along(ages)]
  d <- l * q
  # curtate expectation e_x = sum_{k>=1} l_{x+k} / l_x ; complete ~ e_x + 0.5
  e_curt <- rev(cumsum(rev(c(l[-1], 0)))) / l
  data.frame(age = ages, qx = q, lx = l, dx = d, ex = e_curt + 0.5)
}

# ---- Life-contingent present values (annual, payments at start/end) ---
kpx <- function(lt, x) lt$lx[lt$age >= x] / lt$lx[lt$age == x]

annuity_due <- function(lt, x, i, n = Inf) {          # a-double-dot_{x(:n)}
  p <- kpx(lt, x); k <- seq_along(p) - 1
  keep <- k < n
  sum((1 + i)^(-k[keep]) * p[keep])
}

term_assurance <- function(lt, x, n, i) {             # A^1_{x:n}, paid end of year of death
  p <- kpx(lt, x); q <- lt$qx[lt$age >= x]; k <- seq_along(p) - 1
  keep <- k < n
  sum((1 + i)^(-(k[keep] + 1)) * p[keep] * q[keep])
}

# ---- Graduation tests --------------------------------------------------
grad_tests <- function(D, E, mu, n_par) {
  ex  <- E * mu
  z   <- (D - ex) / sqrt(ex)
  m   <- length(z)
  chi <- sum(z^2)

  n_pos <- sum(z > 0); n_neg <- sum(z < 0)
  runs  <- rle(z > 0)
  G     <- sum(runs$values)                          # number of positive groups
  EG    <- n_pos * (n_neg + 1) / (n_pos + n_neg)
  VG    <- (n_pos * n_neg)^2 / ((n_pos + n_neg)^2 * (n_pos + n_neg - 1))
  cumz  <- sum(D - ex) / sqrt(sum(ex))

  data.frame(
    test = c("Chi-square (overall fit)",
             "Standardised deviations: |z| > 2",
             "Signs test (positives)",
             "Cumulative deviations",
             "Grouping of signs (positive runs)"),
    statistic = c(sprintf("%.1f on %d df", chi, m - n_par),
                  sprintf("%d of %d (expect ~%.1f)", sum(abs(z) > 2), m, 0.0455 * m),
                  sprintf("%d of %d", n_pos, m),
                  sprintf("z = %.2f", cumz),
                  sprintf("%d groups (expect %.1f)", G, EG)),
    p_value = c(pchisq(chi, m - n_par, lower.tail = FALSE),
                pbinom(sum(abs(z) > 2) - 1, m, 0.0455, lower.tail = FALSE),  # one-sided: too many
                min(1, 2 * pbinom(min(n_pos, n_neg), m, 0.5)),
                2 * pnorm(-abs(cumz)),
                pnorm(G + 0.5, EG, sqrt(VG)))          # one-sided: too few groups
  )
}
