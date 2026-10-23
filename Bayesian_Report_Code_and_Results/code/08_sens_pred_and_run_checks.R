# (a) 2030 predictions under each prior setting, (b) evidence for the choice of
# burn-in, number of chains and thinning (early-chain R-hat, autocorrelation, ESS per kept draw)
source("00_common.R")
source("utilities.R"); openGraph <- function(...) invisible(NULL)
fits <- list(S1 = readRDS("out/fit_main.rds"), S2 = readRDS("out/fit_S2.rds"), S3 = readRDS("out/fit_S3.rds"), S4 = readRDS("out/fit_S4.rds"))
pc <- c("Australia", "India", "Nigeria", "Germany", "United States", "Brazil")
out <- do.call(rbind, lapply(names(fits), function(s) {
  M <- as.matrix(fits[[s]]$samples)
  rows <- c(paste0("ypred[", 1:6, "]"), paste0("yNew[", 1:6, "]"))
  tg <- c(pc, paste("New country,", conts))
  do.call(rbind, lapply(seq_along(rows), function(i) { h <- HDIofMCMC(M[, rows[i]])
    data.frame(setting = s, target = tg[i], mean = mean(M[, rows[i]]), hdi_low = h[1], hdi_high = h[2]) }))
}))
write.csv(out, "out/sens_pred.csv", row.names = FALSE)
w <- reshape(out, idvar = "target", timevar = "setting", direction = "wide")
print(round(w[, c("mean.S1", "mean.S2", "mean.S3", "mean.S4")], 2))
cat("max abs diff in predictive mean vs S1:", max(abs(w[, c("mean.S2", "mean.S3", "mean.S4")] - w$mean.S1)), "\n")
cat("max diff in HDI width ratio:", max(abs(((w$hdi_high.S2 - w$hdi_low.S2) / (w$hdi_high.S1 - w$hdi_low.S1)) - 1),
    abs(((w$hdi_high.S3 - w$hdi_low.S3) / (w$hdi_high.S1 - w$hdi_low.S1)) - 1), abs(((w$hdi_high.S4 - w$hdi_low.S4) / (w$hdi_high.S1 - w$hdi_low.S1)) - 1)), "\n")

# run-length / burn-in / thinning evidence on the main fit
S <- fits$S1$samples
key <- c(paste0(c("alpha", "beta", "q", "g"), rep(paste0("[", 1:6, "]"), each = 4)), "sigma", "sigma_a", "sigma_b")
key <- c(paste0("alpha[", 1:6, "]"), paste0("beta[", 1:6, "]"), paste0("q[", 1:6, "]"), paste0("g[", 1:6, "]"), "sigma", "sigma_a", "sigma_b")
early <- sapply(c(50, 100, 250, 500, 1000), function(n) max(gelman.diag(window(S[, key], start = start(S), end = start(S) + (n - 1) * 2), multivariate = FALSE)$psrf[, 1]))
names(early) <- paste0("first_", c(50, 100, 250, 500, 1000), "_kept")
ac1 <- sapply(key, function(k) mean(sapply(1:3, function(ch) acf(as.numeric(S[[ch]][, k]), plot = FALSE, lag.max = 1)$acf[2])))
ess <- effectiveSize(S[, key]); ratio <- ess / 15000
res <- data.frame(quantity = c(names(early), "max_lag1_acf_thinned", "median_lag1_acf_thinned", "min_ess_per_draw", "median_ess_per_draw"),
                  value = c(early, max(ac1), median(ac1), min(ratio), median(ratio)))
write.csv(res, "out/run_checks.csv", row.names = FALSE); print(res)

