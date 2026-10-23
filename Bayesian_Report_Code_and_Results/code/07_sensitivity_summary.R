# 07_sensitivity_summary.R - compare posteriors across prior settings (S1-S4)
source("00_common.R")
suppressMessages({ library(ggplot2) })
source("utilities.R"); openGraph <- function(...) invisible(NULL)
d <- load_data()
fits <- list(S1 = readRDS("out/fit_main.rds"), S2 = readRDS("out/fit_S2.rds"), S3 = readRDS("out/fit_S3.rds"),
             S4 = readRDS("out/fit_S4.rds"))
pn <- function(b) paste0(b, "[", 1:6, "]")
key <- c(pn("alpha"), pn("beta"), pn("q"), pn("g"), "sigma", "sigma_a", "sigma_b")
mats <- lapply(fits, function(f) as.matrix(f$samples))
lab <- function(nm) { m <- regmatches(nm, regexec("^(alpha|beta|q|g)\\[([1-6])\\]$", nm))[[1]]
  if (length(m) == 0) nm else paste0(m[2], " (", conts[as.integer(m[3])], ")") }

rows <- do.call(rbind, lapply(names(mats), function(s) do.call(rbind, lapply(key, function(k) {
  x <- mats[[s]][, k]; h <- HDIofMCMC(x)
  data.frame(setting = s, param = k, label = lab(k), mean = mean(x), sd = sd(x), hdi_low = h[1], hdi_high = h[2], row.names = NULL)
}))))
write.csv(rows, "out/sens_long.csv", row.names = FALSE)

# maximum shift of the posterior mean relative to S1, in S1 posterior SD units
s1 <- rows[rows$setting == "S1", ]
shift <- do.call(rbind, lapply(c("S2", "S3", "S4"), function(s) {
  r <- rows[rows$setting == s, ]; z <- (r$mean - s1$mean) / s1$sd
  data.frame(setting = s, max_abs_shift_sd = max(abs(z)), param_max = r$label[which.max(abs(z))],
             max_abs_shift_years = max(abs(r$mean - s1$mean)), width_ratio_median = median((r$hdi_high - r$hdi_low) / (s1$hdi_high - s1$hdi_low)))
}))
write.csv(shift, "out/sens_shift.csv", row.names = FALSE)

# conclusions that matter- do the main hypotheses survive?
hyp <- do.call(rbind, lapply(names(mats), function(s) {
  M <- mats[[s]]; B <- M[, pn("beta")]
  x <- B[, 2] - B[, 3]; h <- HDIofMCMC(x)
  data.frame(setting = s, asia_minus_europe_mean = mean(x), low = h[1], high = h[2], p_pos = mean(x > 0),
             p_g_all_neg = mean(apply(M[, pn("g")], 1, function(r) all(r < 0))), g6_mean = mean(M[, "g[6]"]), g1_mean = mean(M[, "g[1]"]),
             eu_lowest = mean(apply(B, 1, which.min) == 3))
}))
write.csv(hyp, "out/sens_hyp.csv", row.names = FALSE)

# wide table (mean and HDI) for the report
wide <- reshape(rows[, c("setting", "label", "mean", "hdi_low", "hdi_high")], idvar = "label", timevar = "setting", direction = "wide")
write.csv(wide, "out/sens_wide.csv", row.names = FALSE)

# forest plot for beta and g
sel <- rows[grepl("^(beta|g) ", rows$label), ]
sel$label <- factor(sel$label, levels = rev(unique(sel$label)))
sel$setting <- factor(sel$setting, levels = c("S1", "S2", "S3", "S4"))
pal <- c(S1 = "#14232b", S2 = "#2a78d6", S3 = "#1baf7a", S4 = "#eb6834")
p <- ggplot(sel, aes(mean, label, colour = setting)) + geom_pointrange(aes(xmin = hdi_low, xmax = hdi_high), position = position_dodge(0.7), size = 0.3, linewidth = 0.5) +
  scale_colour_manual(values = pal) + facet_wrap(~ ifelse(grepl("^beta", label), "Trend beta (years per decade)", "COVID effect g (years)"), scales = "free_x") +
  labs(x = "Posterior mean and 95% HDI", y = NULL, colour = NULL) + theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"))
ggsave("fig/fig41_sensitivity.png", p, width = 6.7, height = 5.4, dpi = 300, bg = "white")

# efficiency: run time and ESS per second
tt <- function(f) { x <- f$times; names(x) <- sub("\\.elapsed$", "", names(x)); x }
eff <- do.call(rbind, lapply(names(fits), function(s) {
  f <- fits[[s]]; tms <- tt(f); samp <- if ("sampling" %in% names(tms)) tms["sampling"] else NA
  tot <- sum(tms[intersect(names(tms), c("compile_adapt", "burn_in", "sampling"))])
  e <- effectiveSize(f$samples[, key])
  data.frame(setting = s, kept_draws = nrow(mats[[s]]), total_seconds = unname(tot), sampling_seconds = unname(samp),
             ess_min = min(e), ess_median = median(e), min_ess_per_sec = min(e) / tot, median_ess_per_sec = median(e) / tot,
             rhat_max = max(gelman.diag(f$samples[, key], multivariate = FALSE)$psrf[, 1]))
}))
write.csv(eff, "out/efficiency.csv", row.names = FALSE)
print(shift); print(hyp); print(eff)
