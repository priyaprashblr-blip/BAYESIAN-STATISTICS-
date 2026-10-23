# 05_prior_plots.R - shapes of the prior distributions (how sharp or flat they are)
source("00_common.R")
blue <- "#2a78d6"; orange <- "#eb6834"; green <- "#1baf7a"; amber <- "#eda100"; ink <- "#14232b"
dir.create("fig", showWarnings = FALSE)
nums <- list()   # numbers quoted in the report

# normal priors
png("fig/fig07_prior_normal.png", width = 6.7, height = 5.6, units = "in", res = 300, pointsize = 10)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1), mgp = c(2.4, 0.7, 0), cex.main = 1.0)
x <- seq(-300, 300, length.out = 1000)
plot(x, dnorm(x, 0, 100), type = "l", lwd = 3, col = blue, xlab = "Parameter value", ylab = "Density",
     main = "(a) Normal(0, sd = 100) over a wide range")
abline(v = c(-10, 10), lty = 3, col = orange, lwd = 2)
text(25, max(dnorm(x, 0, 100)) * 0.55, "plausible trend\nrange (-10, 10)", col = orange, cex = 0.9, adj = c(0, 0.5))
x <- seq(-10, 10, length.out = 1000)
plot(x, dnorm(x, 3, 2), type = "l", lwd = 3, col = green, ylim = c(0, 0.22), xlab = "Slope beta (years per decade)",
     ylab = "Density", main = "(b) Same range: flat versus informed")
lines(x, dnorm(x, 0, 100), lwd = 3, col = blue); lines(x, dnorm(x, 0, 1000), lwd = 3, col = amber, lty = 2)
legend("topleft", c("Normal(0, sd 100)", "Normal(0, sd 1000)", "Normal(3, sd 2)"), col = c(blue, amber, green),
       lty = c(1, 2, 1), lwd = 3, bty = "n", cex = 0.85)
x <- seq(0, 130, length.out = 1000)
plot(x, dnorm(x, 65, 15), type = "l", lwd = 3, col = green, ylim = c(0, 0.03), xlab = "Baseline alpha (years)", ylab = "Density",
     main = "(c) Baseline level: informed versus flat")
lines(x, dnorm(x, 0, 100), lwd = 3, col = blue)
legend("topright", c("Normal(0, sd 100)", "Normal(65, sd 15)"), col = c(blue, green), lwd = 3, bty = "n", cex = 0.85)
x <- seq(-30, 30, length.out = 2000)
plot(x, dt(x / 5, 3) / 5, type = "l", lwd = 3, col = orange, ylim = c(0, 0.085), xlab = "Parameter value", ylab = "Density",
     main = "(d) Heavy tails: t(3 df, scale 5) vs Normal(0, sd 5)")
lines(x, dnorm(x, 0, 5), lwd = 3, col = blue)
legend("topright", c("t, 3 df, scale 5", "Normal, sd 5"), col = c(orange, blue), lwd = 3, bty = "n", cex = 0.85)
dev.off()
nums$p_abs_lt10_n100  <- diff(pnorm(c(-10, 10), 0, 100))
nums$p_abs_lt10_n1000 <- diff(pnorm(c(-10, 10), 0, 1000))
nums$dens_ratio_n100  <- dnorm(6, 0, 100) / dnorm(0, 0, 100)
nums$p_abs_gt15_t3    <- 2 * pt(-15 / 5, 3)
nums$p_abs_gt15_n5    <- 2 * pnorm(-15, 0, 5)
nums$p_alpha_in_20_90_n100 <- diff(pnorm(c(20, 90), 0, 100))
nums$p_alpha_in_20_90_n65  <- diff(pnorm(c(20, 90), 65, 15))

#  gamma / exponential priors on tau
set.seed(11)
png("fig/fig08_prior_gamma.png", width = 6.7, height = 5.6, units = "in", res = 300, pointsize = 10)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1), mgp = c(2.4, 0.7, 0), cex.main = 1.0)
tau <- seq(0.001, 2, length.out = 2000)
plot(tau, dgamma(tau, 0.01, 0.01), type = "l", lwd = 3, col = blue, ylim = c(0, 3), xlab = expression(tau),
     ylab = "Density", main = "(a) Gamma priors on precision (linear axes)")
lines(tau, dgamma(tau, 0.001, 0.001), lwd = 3, col = amber, lty = 2); lines(tau, dexp(tau, 0.1), lwd = 3, col = green)
legend("topright", c("Gamma(0.01, 0.01)", "Gamma(0.001, 0.001)", "Exponential(0.1)"), col = c(blue, amber, green),
       lty = c(1, 2, 1), lwd = 3, bty = "n", cex = 0.85)
tl <- 10^seq(-8, 2, length.out = 600)
plot(tl, dgamma(tl, 0.01, 0.01), log = "xy", type = "l", lwd = 3, col = blue, xlab = expression(tau ~ "(log scale)"),
     ylab = "Density (log scale)", main = "(b) Same priors on log-log axes", ylim = c(1e-4, 1e8))
lines(tl, dgamma(tl, 0.001, 0.001), lwd = 3, col = amber, lty = 2); lines(tl, dexp(tl, 0.1), lwd = 3, col = green)
x <- seq(0, 40, length.out = 800)
plot(x, dexp(x, 0.1), type = "l", lwd = 3, col = green, xlab = expression(tau), ylab = "Density",
     main = "(c) Exponential(0.1) = Gamma(1, 0.1)")
polygon(c(0, x, 40), c(0, dexp(x, 0.1), 0), col = adjustcolor(green, 0.2), border = NA)
# implied priors on sigma = 1/sqrt(tau)
n <- 400000
sg <- list("Gamma(0.01, 0.01)" = 1 / sqrt(rgamma(n, 0.01, 0.01)), "Gamma(0.001, 0.001)" = 1 / sqrt(rgamma(n, 0.001, 0.001)),
           "Exponential(0.1)" = 1 / sqrt(rexp(n, 0.1)))
cl <- c(blue, amber, green); lt <- c(1, 2, 1)
dl <- lapply(sg, function(z) density(log10(z), from = -3, to = 8, adjust = 1.5))
plot(NA, xlim = c(-3, 8), ylim = c(0, 1.08 * max(sapply(dl, function(z) max(z$y)))), xlab = expression(log[10](sigma)),
     ylab = "Density", main = expression("(d) Implied prior on " * sigma))
for (i in 1:3) lines(dl[[i]]$x, dl[[i]]$y, lwd = 3, col = cl[i], lty = lt[i])
abline(v = log10(c(0.1, 100)), lty = 3, lwd = 2, col = ink)
text(2.3, max(sapply(dl, function(z) max(z$y))) * 0.62, "plausible sigma\nlies between\nthe dotted lines", cex = 0.85, adj = c(0, 0.5))
legend("topright", names(sg), col = cl, lty = lt, lwd = 3, bty = "n", cex = 0.85)
dev.off()
pm <- function(sh, ra) diff(pgamma(c(1e-4, 100), sh, ra))     # P(0.1 < sigma < 100)
nums$gamma01_mean <- 1; nums$gamma01_var <- 100
nums$gamma01_p_tau_lt_1e3 <- pgamma(1e-3, 0.01, 0.01)
nums$p_sigma_plaus_gamma01  <- pm(0.01, 0.01)
nums$p_sigma_plaus_gamma001 <- pm(0.001, 0.001)
nums$p_sigma_plaus_exp01    <- diff(pexp(c(1e-4, 100), 0.1))
nums$sigma_median_gamma01   <- median(sg[[1]])
nums$sigma_q05_gamma01 <- unname(quantile(sg[[1]], 0.05)); nums$sigma_q95_gamma01 <- unname(quantile(sg[[1]], 0.95))
nums$sigma_median_exp01 <- median(sg[[3]])

#  prior for the t degrees of freedom
png("fig/fig09_prior_nu.png", width = 6.7, height = 3.4, units = "in", res = 300, pointsize = 10)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1), mgp = c(2.4, 0.7, 0), cex.main = 1.0)
nu <- seq(1, 150, length.out = 1000)
plot(nu, dexp(nu - 1, 1 / 29), type = "l", lwd = 3, col = orange, xlab = "nu (degrees of freedom)", ylab = "Density",
     main = "(a) Prior on nu: 1 + Exponential(1/29)")
abline(v = c(10, 30), lty = 3, lwd = 2, col = ink)
xx <- seq(-8, 8, length.out = 600)
plot(xx, dnorm(xx), type = "l", lwd = 3, col = blue, xlab = "Standardised error", ylab = "Density",
     main = "(b) Normal versus t(3) errors")
lines(xx, dt(xx, 3), lwd = 3, col = orange)
legend("topright", c("Normal", "t, 3 df"), col = c(blue, orange), lwd = 3, bty = "n", cex = 0.85)
dev.off()
nums$p_nu_lt10 <- pexp(9, 1 / 29); nums$p_nu_gt30 <- 1 - pexp(29, 1 / 29)
nums$prior_nu_mean <- 1 + 29

write.csv(data.frame(quantity = names(nums), value = unlist(nums)), "out/prior_numbers.csv", row.names = FALSE)
print(data.frame(quantity = names(nums), value = signif(unlist(nums), 4)))
