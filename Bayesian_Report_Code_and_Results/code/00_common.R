# 00_common.R  -  shared set-up, data preparation and helper functions
# Bayesian analysis of life expectancy by continent (JAGS via rjags/coda)

suppressMessages({ library(rjags); library(coda) })

conts <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")
cols  <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
names(cols) <- conts

# data 
load_data <- function(path = "clean.csv") {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[order(d$entity, d$year), ]
  d$cty   <- as.integer(factor(d$entity))                     # country index 1..J
  d$cont  <- as.integer(factor(d$continent, levels = conts))  # continent index 1..6
  d$t     <- (d$year - 1986.5) / 10                           # centred decades
  d$t2    <- d$t^2
  d$covid <- as.integer(d$year %in% c(2020, 2021))
  d
}

# list that JAGS needs (country -> continent lookup, prediction inputs)
make_jags_data <- function(d, pred_countries = c("Australia", "India", "Nigeria",
                                                 "Germany", "United States", "Brazil"),
                           pred_year = 2030) {
  cc <- unique(d[, c("cty", "cont", "entity")]); cc <- cc[order(cc$cty), ]
  pidx <- match(pred_countries, cc$entity)
  stopifnot(!any(is.na(pidx)))
  list(N = nrow(d), J = max(d$cty), K = length(conts),
       y = d$life_expectancy_0, t = d$t, t2 = d$t2, covid = d$covid,
       cty = d$cty, cont = d$cont, ccont = cc$cont,
       P = length(pidx), cpred = pidx, contpred = cc$cont[pidx],
       tNew = (pred_year - 1986.5) / 10)
}

# over-dispersed initial values (one list per chain) 
make_inits <- function(chain) {
  base <- list(
    list(alpha = rep(50, 6), beta = rep(1, 6), q = rep(0, 6),    g = rep(0, 6),
         tau = 0.10, tau_a = 0.05, tau_b = 1.0),
    list(alpha = rep(70, 6), beta = rep(5, 6), q = rep(-0.5, 6), g = rep(-1, 6),
         tau = 0.50, tau_a = 0.01, tau_b = 0.2),
    list(alpha = rep(85, 6), beta = rep(3, 6), q = rep(0.3, 6),  g = rep(1, 6),
         tau = 0.02, tau_a = 0.10, tau_b = 3.0))[[chain]]
  c(base, list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2026 + chain))
}

# fit wrapper that records run times 
run_jags <- function(model_text, data, monitor, n_chains = 3, n_adapt = 1000,
                     n_burn = 3000, n_iter = 10000, thin = 2, inits = make_inits) {
  mf <- tempfile(fileext = ".jags"); writeLines(model_text, mf)
  tm <- list()
  tm$compile_adapt <- system.time(
    m <- jags.model(mf, data = data, inits = lapply(1:n_chains, inits),
                    n.chains = n_chains, n.adapt = n_adapt, quiet = TRUE))["elapsed"]
  tm$burn_in <- system.time(update(m, n.iter = n_burn, progress.bar = "none"))["elapsed"]
  tm$sampling <- system.time(
    s <- coda.samples(m, variable.names = monitor, n.iter = n_iter, thin = thin,
                      progress.bar = "none"))["elapsed"]
  tm$total <- tm$compile_adapt + tm$burn_in + tm$sampling
  list(samples = s, model = m, times = unlist(tm),
       settings = c(chains = n_chains, adapt = n_adapt, burn = n_burn,
                    iter = n_iter, thin = thin, kept_total = n_chains * n_iter / thin))
}

#  highest density interval from MCMC draws (Kruschke 2015) 
hdi <- function(x, prob = 0.95) {
  x <- sort(as.numeric(x)); n <- length(x)
  w <- max(1, floor(prob * n))
  gaps <- x[(w + 1):n] - x[1:(n - w)]
  i <- which.min(gaps)
  c(lower = x[i], upper = x[i + w])
}

# posterior mode from a kernel density estimate
post_mode <- function(x) { dd <- density(x); dd$x[which.max(dd$y)] }
