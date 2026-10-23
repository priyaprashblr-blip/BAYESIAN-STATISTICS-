# 04_posterior.R - posterior summaries, MCMC diagnostics, hypotheses, residuals,
# posterior predictive checks and predictions for the main model (normal likelihood, S1 priors)
source("00_common.R")

source("utilities.R")                 
openGraph <- function(...) invisible(NULL)  

d    <- load_data(); J <- max(d$cty)
fit  <- readRDS("out/fit_main.rds"); S <- fit$samples; M <- as.matrix(S)
pn   <- function(b) paste0(b, "[", 1:6, "]")
key  <- c(pn("alpha"), pn("beta"), pn("q"), pn("g"), "sigma", "sigma_a", "sigma_b")
lab  <- function(nm) {
  m <- regmatches(nm, regexec("^(alpha|beta|q|g)\\[([1-6])\\]$", nm))[[1]]
  if (length(m) == 0) return(nm)
  paste0(m[2], " (", conts[as.integer(m[3])], ")")
}
png_open <- function(f, w, h, ps = 10) png(file.path("fig", f), width = w, height = h, units = "in", res = 300, pointsize = ps)
ytxt <- "Life expectancy at birth (years)"

#  1. posterior summary + diagnostics 
ess  <- effectiveSize(S[, key])
gd   <- gelman.diag(S[, key], multivariate = FALSE)$psrf
post <- do.call(rbind, lapply(key, function(k) {
  sp <- summarizePost(M[, k]); sdv <- sd(M[, k])
  gz <- sapply(1:length(S), function(ch) geweke.diag(S[[ch]][, k, drop = FALSE])$z)
  rl <- try(raftery.diag(S[[1]][, k, drop = FALSE], q = 0.025, r = 0.005, s = 0.95)$resmatrix, silent = TRUE)
  hd <- heidel.diag(S[[1]][, k, drop = FALSE])
  data.frame(param = k, label = lab(k), mean = sp["Mean"], sd = sdv, mode = sp["Mode"], median = sp["Median"],
             hdi_low = sp["HDIlow"], hdi_high = sp["HDIhigh"], ess = ess[k], mcse = sdv / sqrt(ess[k]),
             rhat = gd[k, 1], rhat_upper = gd[k, 2], geweke_max_abs_z = max(abs(gz)),
             geweke_n_fail = sum(abs(gz) > 1.96),
             raftery_N = if (inherits(rl, "try-error")) NA else rl[1, "N"],
             raftery_I = if (inherits(rl, "try-error")) NA else rl[1, "I"],
             heidel_stationary = hd[1, "stest"], heidel_halfwidth = hd[1, "htest"], row.names = NULL)
}))
write.csv(post, "out/post_summary_main.csv", row.names = FALSE)

# country-level parameters a_i, b_i (472)
ab <- c(paste0("a[", 1:J, "]"), paste0("b[", 1:J, "]"))
gd_ab  <- gelman.diag(S[, ab], multivariate = FALSE)$psrf[, 1]
ess_ab <- effectiveSize(S[, ab])
gz_ab  <- sapply(ab, function(k) sapply(1:length(S), function(ch) geweke.diag(S[[ch]][, k, drop = FALSE])$z))
write.csv(data.frame(n_params = length(ab), rhat_max = max(gd_ab), rhat_gt_105 = sum(gd_ab > 1.05), rhat_gt_101 = sum(gd_ab > 1.01),
                     ess_min = min(ess_ab), ess_median = median(ess_ab), ess_lt_1000 = sum(ess_ab < 1000),
                     geweke_frac_fail = mean(abs(gz_ab) > 1.96), n_tests = length(gz_ab)),
          "out/diag_country.csv", row.names = FALSE)

# run times
tm <- fit$times; names(tm) <- sub("\\.elapsed$", "", names(tm))
write.csv(data.frame(stage = names(tm), seconds = unname(tm)), "out/runtime_main.csv", row.names = FALSE)
write.csv(data.frame(setting = names(fit$settings), value = unname(fit$settings)), "out/mcmc_settings.csv", row.names = FALSE)
write.csv(data.frame(measure = names(fit$dic), value = unname(fit$dic)), "out/dic_main.csv", row.names = FALSE)

# 2. posterior distribution plots (plotPost)
pp_grid <- function(file, params, titles, xlab, nrow, ncol, h, compVal = NULL, cenTend = "mode") {
  png_open(file, 6.7, h)
  par(mfrow = c(nrow, ncol), mar = c(4.2, 1, 2.6, 1), mgp = c(2.4, 0.8, 0))
  for (i in seq_along(params)) {
    x <- if (is.character(params[[i]])) M[, params[[i]]] else params[[i]]
    plotPost(x, xlab = xlab, main = titles[i], cex = 0.8, cex.lab = 1.0, cex.main = 1.0, compVal = compVal, cenTend = cenTend)
  }
  dev.off()
}
pp_grid("fig10_post_beta.png", pn("beta"), conts, "Trend beta (years per decade, centred at 1986)", 3, 2, 6.6)
pp_grid("fig11_post_g.png",    pn("g"),    conts, "COVID-19 effect g (years)", 3, 2, 6.6, compVal = 0)
pp_grid("fig12_post_alpha.png", pn("alpha"), conts, "Baseline alpha (years, at 1986)", 3, 2, 6.6)
pp_grid("fig13_post_q.png",    pn("q"),    conts, "Curvature q (years per decade squared)", 3, 2, 6.6, compVal = 0)
png_open("fig14_post_sigma.png", 6.7, 2.6)
par(mfrow = c(1, 3), mar = c(4.2, 1, 2.6, 1), mgp = c(2.4, 0.8, 0))
for (i in 1:3) plotPost(M[, c("sigma", "sigma_a", "sigma_b")[i]], cex = 0.8, cex.lab = 1.0, cex.main = 1.0,
     xlab = c("sigma (years)", "sigma_a (years)", "sigma_b (years/decade)")[i],
     main = c("Residual SD", "SD of country baselines", "SD of country slopes")[i])
dev.off()

# 3. hypotheses
Bm <- M[, pn("beta")]; colnames(Bm) <- conts
pr <- combn(6, 2)
contr <- do.call(rbind, lapply(1:ncol(pr), function(i) {
  a <- pr[1, i]; b <- pr[2, i]; x <- Bm[, a] - Bm[, b]; h <- HDIofMCMC(x)
  data.frame(contrast = paste(conts[a], "-", conts[b]), mean = mean(x), hdi_low = h[1], hdi_high = h[2],
             p_pos = mean(x > 0), excludes_zero = (h[1] > 0 | h[2] < 0), row.names = NULL)
}))
write.csv(contr, "out/hyp_pairwise.csv", row.names = FALSE)
cov_t <- do.call(rbind, lapply(1:6, function(k) {
  x <- M[, paste0("g[", k, "]")]; h <- HDIofMCMC(x)
  data.frame(continent = conts[k], mean = mean(x), hdi_low = h[1], hdi_high = h[2], p_neg = mean(x < 0),
             excludes_zero = (h[1] > 0 | h[2] < 0), row.names = NULL)
}))
write.csv(cov_t, "out/hyp_covid.csv", row.names = FALSE)
qq <- do.call(rbind, lapply(1:6, function(k) {
  x <- M[, paste0("q[", k, "]")]; h <- HDIofMCMC(x)
  data.frame(continent = conts[k], mean = mean(x), hdi_low = h[1], hdi_high = h[2], p_neg = mean(x < 0), row.names = NULL)
}))
write.csv(qq, "out/hyp_curvature.csv", row.names = FALSE)
eu_lowest <- mean(apply(Bm, 1, function(r) which.min(r) == 3))
as_highest <- mean(apply(Bm, 1, function(r) which.max(r) == 2))
write.csv(data.frame(quantity = c("P(beta_Europe is the smallest trend)", "P(beta_Asia is the largest trend)"),
                     value = c(eu_lowest, as_highest)), "out/hyp_misc.csv", row.names = FALSE)

png_open("fig15_contrasts.png", 6.7, 5.4)
par(mfrow = c(2, 2), mar = c(4.2, 1, 2.6, 1), mgp = c(2.4, 0.8, 0))
for (cc in list(c("Asia", "Europe"), c("Africa", "Europe"), c("Asia", "Africa"), c("South America", "Europe")))
  plotPost(Bm[, cc[1]] - Bm[, cc[2]], compVal = 0, cex = 0.8, cex.lab = 1.0, cex.main = 1.0,
           xlab = "Difference in trend (years per decade)", main = paste(cc[1], "minus", cc[2]))
dev.off()

# 4. MCMC diagnostic figures
worst_ess  <- key[which.min(ess[key])]; worst_rhat <- key[which.max(gd[key, 1])]
diag_pars  <- unique(c("beta[2]", "g[6]", worst_ess, worst_rhat))[1:min(4, length(unique(c("beta[2]", "g[6]", worst_ess, worst_rhat))))]
for (p in diag_pars) {
  png_open(paste0("fig20_diag_", gsub("[^A-Za-z0-9]", "_", p), ".png"), 7, 5)
  diagMCMC(S, parName = p)
  dev.off()
}
write.csv(data.frame(param = diag_pars), "out/diag_params_shown.csv", row.names = FALSE)

png_open("fig21_geweke.png", 6.7, 5.4)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
for (p in diag_pars) geweke.plot(S[[1]][, p, drop = FALSE], frac1 = 0.1, frac2 = 0.5, nbins = 20, pvalue = 0.05, auto.layout = FALSE)
dev.off()

png_open("fig22_trace_all.png", 6.7, 8.6, 8)
par(mfrow = c(7, 4), mar = c(2.2, 2.2, 1.6, 0.4), mgp = c(1.2, 0.4, 0), tcl = -0.25)
for (k in key) traceplot(S[, k], main = lab(k), cex.main = 0.75, col = c("skyblue", "black", "royalblue"), xlab = "", ylab = "")
dev.off()

dfk <- data.frame(label = sapply(key, lab), ess = ess[key], rhat = gd[key, 1]); dfk$label <- factor(dfk$label, levels = rev(dfk$label))
thm <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5))
p1 <- ggplot(dfk, aes(ess, label)) + geom_point(colour = "#2a78d6", size = 2) + geom_vline(xintercept = 10000, linetype = 2, colour = "#eb6834") +
  labs(title = "(a) ESS, 28 key parameters", x = "Effective sample size", y = NULL) + thm
p2 <- ggplot(dfk, aes(rhat, label)) + geom_point(colour = "#2a78d6", size = 2) + geom_vline(xintercept = 1.01, linetype = 2, colour = "#eb6834") +
  labs(title = "(b) Gelman-Rubin R-hat", x = "R-hat (point estimate)", y = NULL) + thm + theme(axis.text.y = element_blank())
p3 <- ggplot(data.frame(e = ess_ab), aes(e)) + geom_histogram(bins = 30, fill = "#2a78d6", colour = "white", linewidth = .2) +
  labs(title = "(c) ESS, 472 country parameters", x = "Effective sample size", y = "Count") + thm
p4 <- ggplot(data.frame(r = gd_ab), aes(r)) + geom_histogram(bins = 30, fill = "#2a78d6", colour = "white", linewidth = .2) +
  labs(title = "(d) R-hat, 472 country parameters", x = "R-hat (point estimate)", y = "Count") + thm
ggsave("fig/fig23_ess_rhat.png", (p1 + p2) / (p3 + p4), width = 6.7, height = 7.4, dpi = 300, bg = "white")

# 5. residual analysis (posterior-mean fit)
ahat <- colMeans(M[, paste0("a[", 1:J, "]")]); bhat <- colMeans(M[, paste0("b[", 1:J, "]")])
qh <- colMeans(M[, pn("q")]); gh <- colMeans(M[, pn("g")]); sg <- mean(M[, "sigma"])
d$fit  <- ahat[d$cty] + bhat[d$cty] * d$t + qh[d$cont] * d$t2 + gh[d$cont] * d$covid
d$res  <- d$life_expectancy_0 - d$fit; d$sres <- d$res / sg
lag1 <- sapply(split(d$res, d$cty), function(r) cor(r[-1], r[-length(r)]))
rm_c <- sort(sqrt(tapply(d$res^2, d$entity, mean)), decreasing = TRUE)
ry   <- tapply(d$res, d$year, mean)
write.csv(data.frame(quantity = c("sd_std_resid", "mean_resid", "frac_gt2", "frac_gt3", "excess_kurtosis", "skewness",
                                  "lag1_mean", "lag1_median", "mean_resid_min", "mean_resid_max",
                                  "year_of_min", "year_of_max", "fitted_min", "fitted_max", "rmse"),
                     value = c(sd(d$sres), mean(d$res), mean(abs(d$sres) > 2), mean(abs(d$sres) > 3), kurtosis(d$sres, type = 2),
                               skewness(d$sres, type = 2), mean(lag1), median(lag1), min(ry), max(ry),
                               as.numeric(names(ry)[which.min(ry)]), as.numeric(names(ry)[which.max(ry)]),
                               min(d$fit), max(d$fit), sqrt(mean(d$res^2)))), "out/resid_summary.csv", row.names = FALSE)
write.csv(data.frame(country = names(rm_c)[1:8], rmse = unname(rm_c)[1:8],
                     continent = d$continent[match(names(rm_c)[1:8], d$entity)]), "out/resid_worst.csv", row.names = FALSE)
write.csv(data.frame(year = as.numeric(names(ry)), mean_resid = unname(ry)), "out/resid_by_year.csv", row.names = FALSE)
d$continent <- factor(d$continent, levels = conts)
cl <- setNames(cols, conts)
r1 <- ggplot(d, aes(fit, sres)) + geom_point(size = 0.35, alpha = 0.15, colour = "#2a78d6") +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), linetype = c(2, 3, 1, 3, 2), colour = c("#eb6834", "#eb6834", "grey30", "#eb6834", "#eb6834")) +
  labs(title = "(a) Standardised residuals vs fitted", x = "Fitted life expectancy (years)", y = "Standardised residual") + thm
r2 <- ggplot(d, aes(sample = sres)) + stat_qq(size = 0.35, alpha = 0.3, colour = "#2a78d6") + stat_qq_line(colour = "#eb6834", linewidth = 0.8) +
  labs(title = "(b) Normal Q-Q plot", x = "Theoretical quantile", y = "Sample quantile") + thm
mr <- aggregate(res ~ year + continent, d, mean)
r3 <- ggplot(mr, aes(year, res, colour = continent)) + annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf, fill = "#eda100", alpha = .25) +
  geom_hline(yintercept = 0, colour = "grey40") + geom_line(linewidth = 0.7) + scale_colour_manual(values = cl) +
  labs(title = "(c) Mean residual by year and continent", x = "Year", y = "Mean residual (years)", colour = NULL) + thm +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5), legend.key.size = unit(0.3, "cm")) + guides(colour = guide_legend(nrow = 2))
r4 <- ggplot(d, aes(sres)) + geom_histogram(aes(y = after_stat(density)), binwidth = 0.25, fill = "#2a78d6", colour = "white", linewidth = 0.1) +
  stat_function(fun = dnorm, colour = "#eb6834", linewidth = 0.9) + coord_cartesian(xlim = c(-8, 8)) +
  labs(title = "(d) Histogram vs standard normal", x = "Standardised residual", y = "Density") + thm
ggsave("fig/fig30_residuals.png", (r1 + r2) / (r3 + r4), width = 6.7, height = 6.4, dpi = 300, bg = "white")

#  6. posterior predictive checks
ymin <- min(d$life_expectancy_0); ymax <- max(d$life_expectancy_0); ysd <- sd(d$life_expectancy_0)
ppc <- data.frame(statistic = c("minimum", "maximum", "standard deviation", "share beyond 3 SD"),
  observed = c(ymin, ymax, ysd, mean(M[, "prop_big_obs"])),
  rep_mean = c(mean(M[, "yrep_min"]), mean(M[, "yrep_max"]), mean(M[, "yrep_sd"]), mean(M[, "prop_big_rep"])),
  p_value  = c(mean(M[, "yrep_min"] <= ymin), mean(M[, "yrep_max"] >= ymax), mean(M[, "yrep_sd"] >= ysd),
               mean(M[, "prop_big_rep"] >= M[, "prop_big_obs"])))
write.csv(ppc, "out/ppc.csv", row.names = FALSE)
png_open("fig31_ppc.png", 6.7, 5.4)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.6, 1), mgp = c(2.4, 0.8, 0), cex.main = 1.0)
for (i in 1:3) {
  v <- M[, c("yrep_min", "yrep_max", "yrep_sd")[i]]; o <- c(ymin, ymax, ysd)[i]
  hist(v, breaks = 40, col = "skyblue", border = "white", xlim = range(c(v, o)), main = paste0("(", letters[i], ") Replicated ", ppc$statistic[i]),
       xlab = c("Minimum (years)", "Maximum (years)", "Standard deviation (years)")[i], ylab = "Frequency")
  abline(v = o, col = "#eb6834", lwd = 3); legend("topright", sprintf("observed %.2f\nP = %.3f", o, ppc$p_value[i]), bty = "n", cex = 0.85, text.col = "#eb6834")
}
plot(M[, "prop_big_obs"], M[, "prop_big_rep"], pch = 16, cex = 0.3, col = rgb(0.16, 0.47, 0.84, 0.25),
     xlab = "Share beyond 3 SD, observed data", ylab = "Share beyond 3 SD, replicated data", main = "(d) Tail share per draw")
abline(0, 1, col = "#eb6834", lwd = 2); legend("topleft", sprintf("P(rep >= obs) = %.3f", ppc$p_value[4]), bty = "n", cex = 0.85)
dev.off()

# 7. predictions for 2030
jd <- make_jags_data(d)
pc <- c("Australia", "India", "Nigeria", "Germany", "United States", "Brazil")
obs23 <- sapply(pc, function(z) d$life_expectancy_0[d$entity == z & d$year == 2023])
pt <- do.call(rbind, lapply(1:6, function(p) {
  x <- M[, paste0("ypred[", p, "]")]; h <- HDIofMCMC(x); ch <- x - obs23[p]; hc <- HDIofMCMC(ch)
  data.frame(target = pc[p], continent = conts[jd$contpred[p]], observed_2023 = obs23[p], mean = mean(x), hdi_low = h[1], hdi_high = h[2],
             change_mean = mean(ch), change_low = hc[1], change_high = hc[2], p_up = mean(ch > 0), row.names = NULL)
}))
pn2 <- do.call(rbind, lapply(1:6, function(k) {
  x <- M[, paste0("yNew[", k, "]")]; h <- HDIofMCMC(x)
  data.frame(target = paste("New country,", conts[k]), continent = conts[k], observed_2023 = NA, mean = mean(x), hdi_low = h[1], hdi_high = h[2],
             change_mean = NA, change_low = NA, change_high = NA, p_up = NA, row.names = NULL)
}))
write.csv(rbind(pt, pn2), "out/pred_table.csv", row.names = FALSE)
pa <- ggplot(pt, aes(mean, factor(target, levels = rev(pc)))) + geom_errorbarh(aes(xmin = hdi_low, xmax = hdi_high), height = 0.25, colour = "#2a78d6", linewidth = 1) +
  geom_point(colour = "#2a78d6", size = 2.6) + geom_point(aes(observed_2023), shape = 18, size = 3.6, colour = "#eb6834") +
  labs(title = "(a) Named countries, 2030", x = ytxt, y = NULL) + thm
pb <- ggplot(pn2, aes(mean, factor(continent, levels = rev(conts)))) + geom_errorbarh(aes(xmin = hdi_low, xmax = hdi_high), height = 0.25, colour = "#1baf7a", linewidth = 1) +
  geom_point(colour = "#1baf7a", size = 2.6) + labs(title = "(b) A new, unseen country, 2030", x = ytxt, y = NULL) + thm
ggsave("fig/fig50_predictions.png", pa / pb, width = 6.7, height = 5.6, dpi = 300, bg = "white")

#  8. prior versus posterior
png_open("fig40_prior_vs_post.png", 6.7, 5.6)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1), mgp = c(2.4, 0.8, 0), cex.main = 1.0)
pv <- function(x, prior, xlim, main, xlab, lg = "Prior N(0, sd 100)") {
  dd <- density(x); xs <- seq(xlim[1], xlim[2], length.out = 600)
  plot(dd, lwd = 3, col = "#2a78d6", xlim = xlim, main = main, xlab = xlab, ylab = "Density", ylim = c(0, max(dd$y) * 1.08))
  lines(xs, prior(xs), lwd = 3, col = "#eb6834", lty = 2)
  legend("topright", c("Posterior", lg), col = c("#2a78d6", "#eb6834"), lwd = 3, lty = c(1, 2), bty = "n", cex = 0.8)
}
pv(M[, "beta[2]"], function(x) dnorm(x, 0, 100), c(3, 5.5), "(a) Trend, Asia", "beta (years per decade)")
pv(M[, "g[6]"], function(x) dnorm(x, 0, 100), c(-5, 0), "(b) COVID effect, South America", "g (years)")
pv(M[, "alpha[3]"], function(x) dnorm(x, 0, 100), c(70, 78), "(c) Baseline, Europe", "alpha (years)")
fs <- function(s) dgamma(s^-2, 0.01, 0.01) * 2 * s^-3
pv(M[, "sigma"], fs, c(2.3, 2.5), "(d) Residual SD", "sigma (years)", "Prior on sigma implied by Gamma(0.01, 0.01)")
dev.off()

print(post[, c("label", "mean", "hdi_low", "hdi_high", "ess", "rhat", "geweke_max_abs_z", "raftery_I")], digits = 4)
print(read.csv("out/diag_country.csv")); print(contr); print(cov_t); print(ppc); print(rbind(pt, pn2))
print(read.csv("out/resid_summary.csv")); print(tm); print(fit$dic); print(diag_pars)
