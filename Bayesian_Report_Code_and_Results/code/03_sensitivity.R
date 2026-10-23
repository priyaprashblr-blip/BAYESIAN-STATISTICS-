# 03_sensitivity.R - prior sensitivity runs S2-S4 (same model, different priors)
# S1 (primary priors, normal likelihood) is the main fit from 02_fit_main.R.
source("00_common.R")
d   <- load_data()
dat <- make_jags_data(d)

base <- readLines("model_main.jags")
i0 <- grep("#PRIORS-BEGIN", base); i1 <- grep("#PRIORS-END", base)
ppc <- grep("posterior predictive checks", base)
lean <- c(base[1:(ppc - 1)], "}")                      # drop the (slow) replicate-data block
swap_priors <- function(model, block) {
  j0 <- grep("#PRIORS-BEGIN", model); j1 <- grep("#PRIORS-END", model)
  c(model[1:j0], block, model[j1:length(model)])
}

priors <- list(
  S2 = c("  for (k in 1:K) {",                       # much wider
         "    alpha[k] ~ dnorm(0, 1.0E-6)", "    beta[k]  ~ dnorm(0, 1.0E-6)",
         "    q[k]     ~ dnorm(0, 1.0E-6)", "    g[k]     ~ dnorm(0, 1.0E-6)", "  }",
         "  tau   ~ dgamma(0.001, 0.001)", "  tau_a ~ dgamma(0.001, 0.001)", "  tau_b ~ dgamma(0.001, 0.001)"),
  S3 = c("  for (k in 1:K) {",                       # informed (centred, moderate sd)
         "    alpha[k] ~ dnorm(65, 1/225)",          # sd 15
         "    beta[k]  ~ dnorm(3, 1/4)",             # sd 2
         "    q[k]     ~ dnorm(0, 1)",               # sd 1
         "    g[k]     ~ dnorm(0, 1/4)", "  }",      # sd 2
         "  tau   ~ dgamma(0.01, 0.01)", "  tau_a ~ dgamma(0.01, 0.01)", "  tau_b ~ dgamma(0.01, 0.01)"),
  S4 = c("  for (k in 1:K) {",                       # heavy-tailed t priors + exponential precisions
         "    alpha[k] ~ dt(0, 1/10000, 3)",         # scale 100, 3 df
         "    beta[k]  ~ dt(0, 1/25, 3)",            # scale 5
         "    q[k]     ~ dt(0, 1, 3)",               # scale 1
         "    g[k]     ~ dt(0, 1/25, 3)", "  }",     # scale 5
         "  tau   ~ dexp(0.1)", "  tau_a ~ dexp(0.1)", "  tau_b ~ dexp(0.1)"))

core <- c("alpha", "beta", "q", "g", "sigma", "sigma_a", "sigma_b", "ypred", "yNew")

fit_variant <- function(name, model, monitor, with_dic = FALSE) {
  fit <- run_jags(paste(model, collapse = "\n"), dat, monitor, n_chains = 3, n_adapt = 1000,
                  n_burn = 3000, n_iter = 6000, thin = 3)
  res <- list(samples = fit$samples, times = fit$times, settings = fit$settings)
  if (with_dic) {
    load.module("dic")
    td <- system.time(dic <- dic.samples(fit$model, n.iter = 3000, type = "pD"))["elapsed"]
    res$times <- c(res$times, dic = unname(td))
    res$dic <- c(mean_deviance = sum(dic$deviance), pD = sum(dic$penalty),
                 DIC = sum(dic$deviance) + sum(dic$penalty))
  }
  saveRDS(res, file.path("out", paste0("fit_", name, ".rds")))
  writeLines(paste(model, collapse = "\n"), file.path("out", paste0("model_", name, ".jags")))
  cat(name, "done; total seconds:", res$times["total"], "\n"); flush.console()
}

for (nm in names(priors)) fit_variant(nm, swap_priors(lean, priors[[nm]]), core)

cat("ALL DONE\n")
