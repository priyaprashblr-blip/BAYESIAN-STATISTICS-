# Generate MCMC samples
genMCMC = function(data, yName = "life_expectancy_0", sName = "entity", cName = "continent",
                   xName = "year", predCountries = c("Australia", "India", "Nigeria",
                                               "Germany", "United States", "Brazil"),
                   predYear = 2030, numSavedSteps = 15000, thinSteps = 2, nChains = 3,
                   adaptSteps = 1000, burnInSteps = 3000, includePPC = TRUE, useRunjags = TRUE,
                   saveName = NULL) {
  require(rjags)
  
  # Prepare the observed data
  y = data[, yName] 
  yr = data[, xName] # Year
  
  if (any(!is.finite(y))) stop("All y values must be finite.")
  
  # Set the continent order
  conts = c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")
  
  # Give each country and continent a number
  sf = factor(data[, sName])
  s = as.integer(sf)
  ct = as.integer(factor(data[, cName], levels = conts))
  if (any(is.na(ct))) stop("Unknown continent label.")
  
  J = max(s) 
  K = length(conts) # Number of countries / continents
  
  # Match each country to its continent
  sc = unique(data.frame(s = s, ct = ct, name = as.character(sf)))
  sc = sc[order(sc$s), ]
  
  # Centre time around 1986.5 and measure it in decades
  t = (yr - 1986.5) / 10
  
  # Identify the COVID-19 years
  cv = as.integer(yr %in% c(2020, 2021))
  
  # Find the countries used for predictions
  pi = match(predCountries, sc$name)
  if (any(is.na(pi))) stop("Unknown prediction country.")
  
  # Send the data to JAGS
  jd = list(
    y = y, t = t, t2 = t^2, covid = cv,
    cty = s, cont = ct, ccont = sc$ct,
    N = length(y), J = J, K = K,
    P = length(pi), cpred = pi, contpred = sc$ct[pi],
    tNew = (predYear - 1986.5) / 10
  )
  
  # Define the Bayesian model
  mod = "
  model {
    # Model the observed life expectancy
    for (n in 1:N) {
      y[n] ~ dnorm(mu[n], tau)
      mu[n] <- a[cty[n]] + b[cty[n]] * t[n] +
               q[cont[n]] * t2[n] + g[cont[n]] * covid[n]
    }
    
    # Country baselines and trends are partially pooled
    # towards the average values for their continent
    for (i in 1:J) {
      a[i] ~ dnorm(alpha[ccont[i]], tau_a)
      b[i] ~ dnorm(beta[ccont[i]], tau_b)
    }
    
    # Priors for the continent-level parameters
    for (k in 1:K) {
      alpha[k] ~ dnorm(0, 1.0E-4)
      beta[k] ~ dnorm(0, 1.0E-4)
      q[k] ~ dnorm(0, 1.0E-4)
      g[k] ~ dnorm(0, 1.0E-4)
    }
    
    # Priors for the precision parameters
    tau ~ dgamma(0.01, 0.01)
    tau_a ~ dgamma(0.01, 0.01)
    tau_b ~ dgamma(0.01, 0.01)
    
    # Convert precision into standard deviations
    sigma <- 1 / sqrt(tau)
    sigma_a <- 1 / sqrt(tau_a)
    sigma_b <- 1 / sqrt(tau_b)
    
    # Predict life expectancy in the selected countries
    for (p in 1:P) {
      mupred[p] <- a[cpred[p]] + b[cpred[p]] * tNew +
                   q[contpred[p]] * pow(tNew, 2)
      ypred[p] ~ dnorm(mupred[p], tau)
    }
    
    # Predict life expectancy for a new country on each continent
    for (k in 1:K) {
      aNew[k] ~ dnorm(alpha[k], tau_a)
      bNew[k] ~ dnorm(beta[k], tau_b)
      yNew[k] ~ dnorm(
        aNew[k] + bNew[k] * tNew +
        q[k] * pow(tNew, 2), tau
      )
    }
  }
  "
  
  # Add simulated data to check how well the model fits
  ppc = "
  # Generate replicated observations
  for (n in 1:N) {
    yrep[n] ~ dnorm(mu[n], tau)
    # Identify values more than 3 standard deviations away
    big_rep[n] <- step(abs(yrep[n] - mu[n]) * sqrt(tau) - 3)
    big_obs[n] <- step(abs(y[n] - mu[n]) * sqrt(tau) - 3)
  }
  
  # Summarise the replicated data
  yrep_min <- min(yrep[])
  yrep_max <- max(yrep[])
  yrep_sd <- sd(yrep[])
  prop_big_rep <- mean(big_rep[])
  prop_big_obs <- mean(big_obs[])
  "
  
  # Include the predictive checks when requested
  mod = paste0(mod, if (includePPC) ppc else "", "\n}\n")
  writeLines(mod, con = "TEMPmodel.txt")
  
  # Start chains from different values to check convergence
  ini = lapply(1:nChains, function(ch) {
    base = list(
      list(alpha = rep(50, 6), beta = rep(1, 6), q = rep(0, 6),
           g = rep(0, 6), tau = 0.10, tau_a = 0.05, tau_b = 1.0),
      list(alpha = rep(70, 6), beta = rep(5, 6), q = rep(-0.5, 6),
           g = rep(-1, 6), tau = 0.50, tau_a = 0.01, tau_b = 0.2),
      list(alpha = rep(85, 6), beta = rep(3, 6), q = rep(0.3, 6),
           g = rep(1, 6), tau = 0.02, tau_a = 0.10, tau_b = 3.0)
    )[[(ch - 1) %% 3 + 1]]
    c(base, list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2026 + ch))
  })
  
  # Choose which parameters to save
  pars = c("alpha", "beta", "q", "g", "sigma", "sigma_a",
           "sigma_b", "a", "b", "ypred", "yNew")
  if (includePPC) {
    pars = c(pars, "yrep_min", "yrep_max", "yrep_sd",
             "prop_big_rep", "prop_big_obs")
  }
  
  # Run the chains and record the computing time
  t0 = proc.time()["elapsed"]
  if (useRunjags && requireNamespace("runjags", quietly = TRUE)) {
    # Run chains in parallel when runjags is available
    out = runjags::run.jags(
      method = "parallel", model = "TEMPmodel.txt",
      monitor = pars, data = jd, inits = ini,
      n.chains = nChains, adapt = adaptSteps,
      burnin = burnInSteps,
      sample = ceiling(numSavedSteps / nChains),
      thin = thinSteps, summarise = FALSE, plots = FALSE
    )
    sam = as.mcmc.list(out)
  } else {
    # Otherwise, run the chains sequentially
    jm = jags.model(
      "TEMPmodel.txt", data = jd, inits = ini,
      n.chains = nChains, n.adapt = adaptSteps, quiet = TRUE
    )
    update(jm, n.iter = burnInSteps, progress.bar = "none")
    sam = coda.samples(
      jm, variable.names = pars,
      n.iter = ceiling(numSavedSteps * thinSteps / nChains),
      thin = thinSteps, progress.bar = "none"
    )
  }
  
  # Store the run time and country names
  attr(sam, "seconds") = unname(proc.time()["elapsed"] - t0)
  attr(sam, "countryNames") = sc$name
  
  # Save the MCMC samples if a file name is provided
  if (!is.null(saveName)) {
    save(sam, file = paste(saveName, "Mcmc.Rdata", sep = ""))
  }
  
  return(sam)
}


# Summarise the posterior results
smryMCMC = function(sam, compVal = 0, saveName = NULL) {
  library(coda)
  conts = c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")
  mat = as.matrix(sam)
  
  # Select the main parameters for the summary
  key = c(paste0("alpha[", 1:6, "]"),
          paste0("beta[", 1:6, "]"),
          paste0("q[", 1:6, "]"),
          paste0("g[", 1:6, "]"),
          "sigma", "sigma_a", "sigma_b")
  
  # Check whether the chains have converged
  rh = gelman.diag(sam[, key], multivariate = FALSE)$psrf[, 1]
  res = NULL
  
  for (p in key) {
    # Compare curvature and COVID effects with zero
    cv = if (grepl("^(q|g)\\[", p)) compVal else NULL
    # Calculate the posterior summaries
    sp = summarizePost(mat[, p], compVal = cv)
    res = rbind(res, c(
      sp[c("Mean", "Median", "Mode", "ESS", "HDIlow", "HDIhigh")],
      Rhat = unname(rh[p]),
      PcntGtCompVal = if (is.null(cv)) NA else unname(sp["PcntGtCompVal"])
    ))
  }
  
  rownames(res) = key
  
  # Save the summary table if requested
  if (!is.null(saveName)) {
    write.csv(res, file = paste(saveName, "SummaryInfo.csv", sep = ""))
  }
  
  return(res)
}


# Plot the posterior distributions
plotMCMC = function(sam, data = NULL, compVal = 0,
                    diffCList = list(
                      c("Asia", "Europe"),
                      c("Africa", "Europe"),
                      c("Asia", "Africa"),
                      c("South America", "Europe")
                    ),
                    saveName = NULL, saveType = "png") {
  conts = c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")
  mat = as.matrix(sam)
  
  # Open a file to save each plot
  openPlot = function(suf, w, h) {
    if (!is.null(saveName)) {
      png(file = paste0(saveName, suf, ".", saveType),
          width = w, height = h, units = "in",
          res = 300, pointsize = 10)
    }
  }
  
  # Close the plot file
  closePlot = function() {
    if (!is.null(saveName)) dev.off()
  }
  
  # Plot one parameter for all six continents
  grid6 = function(pre, suf, xlab, cv = NULL) {
    openPlot(suf, 6.7, 6.6)
    par(mfrow = c(3, 2), mar = c(4.2, 1, 2.6, 1), mgp = c(2.4, 0.8, 0))
    for (k in 1:6) {
      plotPost(mat[, paste0(pre, "[", k, "]")],
               xlab = xlab, main = conts[k],
               cex = 0.8, cex.lab = 1.0,
               cex.main = 1.0, compVal = cv)
    }
    closePlot()
  }
  
  # Plot the baseline life expectancy for each continent
  grid6("alpha", "PostAlpha", "Baseline alpha (years, at 1986)")
  
  # Plot the life expectancy trends
  grid6("beta", "PostBeta", "Trend beta (years per decade)")
  
  # Plot the curvature of the trends
  grid6("q", "PostQ", "Curvature q", cv = compVal)
  
  # Plot the COVID-19 effects
  grid6("g", "PostG", "COVID-19 effect g (years)", cv = compVal)
  
  # Plot the estimated standard deviations
  openPlot("PostSigma", 6.7, 2.6)
  par(mfrow = c(1, 3), mar = c(4.2, 1, 2.6, 1), mgp = c(2.4, 0.8, 0))
  for (p in c("sigma", "sigma_a", "sigma_b")) {
    plotPost(mat[, p], xlab = p, main = p,
             cex = 0.8, cex.lab = 1.0, cex.main = 1.0)
  }
  closePlot()
  
  # Compare trends between pairs of continents
  openPlot("PostContrasts", 6.7, 2.8 * length(diffCList))
  par(mfrow = c(length(diffCList), 2),
      mar = c(4.2, 4.2, 2.6, 1), mgp = c(2.4, 0.8, 0))
  for (cc in diffCList) {
    b1 = mat[, paste0("beta[", match(cc[1], conts), "]")]
    b2 = mat[, paste0("beta[", match(cc[2], conts), "]")]
    
    # Plot the difference between the two trends
    plotPost(b1 - b2, compVal = 0,
             xlab = "Difference in trend (years per decade)",
             main = paste(cc[1], "minus", cc[2]),
             cex = 0.8, cex.lab = 1.0, cex.main = 1.0)
             
    # Show the relationship between the two trends
    plot(b2[1:2000], b1[1:2000],
         pch = 1, cex = 0.5, col = "skyblue",
         xlab = cc[2], ylab = cc[1],
         main = "Joint posterior of the two trends")
    abline(0, 1, lty = "dotted")
  }
  closePlot()
}
