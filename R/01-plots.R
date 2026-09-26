# ============================================================================
# Figures 1-3: PDF, CDF, hazard and survival shapes of the MGKw distribution
#
# The original plotting notebook reassigned the parameters a and b inside the
# plotting function, so every curve was drawn with the same values and the
# function arguments were ignored. Parameters are passed explicitly here.
# ============================================================================

source(file.path("R", "mgkw-distribution.R"))

dir.create("figures", showWarnings = FALSE)

x <- seq(0.001, 0.999, length.out = 1000)

# Parameter sets used in Figure 1 of the manuscript
pars_left  <- list(c(1.2, 3.5), c(2, 2), c(0.5, 0.7), c(5.5, 4.05),
                   c(2.0, 2.5), c(1.02, 0.15), c(2.02, 8.05))
pars_right <- list(c(0.5, 7), c(3, 4), c(0.2, 0.5), c(0.05, 0.7),
                   c(0.002, 0.05), c(5, 2), c(3, 8))

plot_family <- function(pars, fun, ylab, main, ylim = NULL, file) {
  cols <- grDevices::rainbow(length(pars))
  curves <- lapply(pars, function(p) fun(x, p[1], p[2]))
  if (is.null(ylim)) ylim <- c(0, min(5, max(unlist(curves), na.rm = TRUE)))
  grDevices::png(file.path("figures", file), width = 900, height = 700, res = 120)
  plot(x, curves[[1]], type = "l", col = cols[1], lwd = 2, ylim = ylim,
       xlab = "x", ylab = ylab, main = main)
  for (i in seq_along(pars)[-1]) lines(x, curves[[i]], col = cols[i], lwd = 2)
  legend("topright", bty = "n", cex = 0.75, col = cols, lwd = 2,
         legend = vapply(pars, function(p) sprintf("a=%g, b=%g", p[1], p[2]),
                         character(1)))
  invisible(grDevices::dev.off())
}

plot_family(pars_left,  dmgkw, "f(x)", "MGKw density", file = "fig1a_pdf.png")
plot_family(pars_right, dmgkw, "f(x)", "MGKw density", file = "fig1b_pdf.png")
plot_family(pars_left,  pmgkw, "F(x)", "MGKw distribution function",
            ylim = c(0, 1), file = "fig_cdf.png")
plot_family(pars_left,  smgkw, "S(x)", "MGKw survival function",
            ylim = c(0, 1), file = "fig_survival.png")

# Hazard shapes (Figure 2): two panels — bathtub (left) and IDI/increasing (right)
haz_left  <- list(c(0.1, 1.5), c(0.1, 1.3), c(0.3, 1.01), c(0.02, 0.23),
                  c(0.002, 0.5))
haz_right <- list(c(1.3, 1.4), c(7.1, 0.3), c(3.3, 1.01), c(0.61, 2.1),
                  c(2, 0.5))
plot_family(haz_left,  hmgkw, "h(x)", "MGKw hazard rate (bathtub shapes)",
            ylim = c(0, 10), file = "fig2a_hazard.png")
plot_family(haz_right, hmgkw, "h(x)", "MGKw hazard rate (IDI/increasing shapes)",
            ylim = c(0, 10), file = "fig2b_hazard.png")

cat("Figures written to figures/\n")
