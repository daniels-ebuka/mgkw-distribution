# ============================================================================
# Monte Carlo study of the MLEs of the MGKw distribution (Table 2)
#
# TRUE PARAMETERS: a = 1.001, b = 0.77
#
# For each sample size the script reports the average estimate, the bias
# (average - true value), the measure of spread (MS, the Monte Carlo
# variance) and the average squared error (ASE = Bias^2 + MS). Bias and ASE
# should shrink towards zero as n grows, as expected of a consistent
# estimator.
#
# NOTE ON set.seed PLACEMENT
# The original scripts called set.seed() inside the replication loop, which
# regenerates an identical sample on every iteration and collapses the Monte
# Carlo variance to zero. The seed is set once, before the loops.
# ============================================================================

source(file.path("R", "mgkw-distribution.R"))

set.seed(2024)          # set ONCE, outside every loop

a_true <- 1.001
b_true <- 0.77

n_vec  <- c(60, 80, 100, 200, 400, 600, 800, 1000)
R      <- 10000         # replications (as in the manuscript)

simulate_one_n <- function(n, R, a_true, b_true) {
  est <- matrix(NA_real_, nrow = R, ncol = 2,
                dimnames = list(NULL, c("a", "b")))
  for (r in seq_len(R)) {
    x <- rmgkw(n, a_true, b_true)
    x <- pmin(pmax(x, 1e-10), 1 - 1e-10)
    fit <- tryCatch(
      fit_mgkw(x, start = c(a_true, b_true)),
      error = function(e) NULL
    )
    if (!is.null(fit) && fit$convergence == 0 && fit$value < 1e9)
      est[r, ] <- fit$par
  }
  est  <- est[complete.cases(est), , drop = FALSE]
  bias <- colMeans(est) - c(a_true, b_true)
  ms   <- apply(est, 2, var)
  ase  <- bias^2 + ms
  data.frame(
    n         = n,
    parameter = c("a", "b"),
    Average   = round(colMeans(est), 4),
    Bias      = round(bias,          4),
    MS        = round(ms,            4),
    ASE       = round(ase,           4),
    converged = sprintf("%d/%d", nrow(est), R),
    row.names = NULL
  )
}

cat(sprintf("Simulation: true a = %g, b = %g, R = %d replications\n",
            a_true, b_true, R))
cat(sprintf("Sample sizes: %s\n\n", paste(n_vec, collapse = ", ")))

results <- do.call(rbind, lapply(n_vec, simulate_one_n,
                                 R = R, a_true = a_true, b_true = b_true))
print(results, digits = 4)

dir.create("output", showWarnings = FALSE)
write.csv(results, file.path("output", "Table2_simulation.csv"),
          row.names = FALSE)

cat("\nConsistency check (|Bias| and ASE should shrink as n grows):\n")
for (p in c("a", "b")) {
  sub <- results[results$parameter == p, ]
  cat(sprintf("  |Bias| for %s: %s\n", p,
              paste(sprintf("%.4f", abs(sub$Bias)), collapse = "  ")))
}
