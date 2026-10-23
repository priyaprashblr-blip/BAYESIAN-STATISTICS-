# 02_fit_main.R - fit the main hierarchical model (normal likelihood, primary priors)
source("00_common.R")
d   <- load_data()
dat <- make_jags_data(d)
model_text <- paste(readLines("model_main.jags"), collapse = "\n")

monitor <- c("alpha", "beta", "q", "g", "sigma", "sigma_a", "sigma_b", "a", "b",
             "ypred", "yNew", "yrep_min", "yrep_max", "yrep_sd",
             "prop_big_rep", "prop_big_obs")

fit <- run_jags(model_text, dat, monitor, n_chains = 3, n_adapt = 1000,
                n_burn = 3000, n_iter = 10000, thin = 2)

# deviance information criterion 
load.module("dic")
t_dic <- system.time(dic <- dic.samples(fit$model, n.iter = 3000, type = "pD"))["elapsed"]

out <- list(samples = fit$samples, times = c(fit$times, dic = unname(t_dic)),
            settings = fit$settings,
            dic = c(mean_deviance = sum(dic$deviance), pD = sum(dic$penalty),
                    DIC = sum(dic$deviance) + sum(dic$penalty)))
saveRDS(out, "out/fit_main.rds")
print(out$times); print(out$settings); print(out$dic)
