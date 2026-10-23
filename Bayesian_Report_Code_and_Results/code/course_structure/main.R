# Main script for the Bayesian life expectancy model
# Reads the data, runs the MCMC chains, checks convergence,
# calculates summary statistics and produces posterior plots.

# Clear existing plots and objects
graphics.off()
rm(list = ls(all = TRUE))

# Load the model functions
source("main.R")

# Read the cleaned life expectancy data
# The file contains country, year, life expectancy and continent
dat = read.csv("/Users/priya/Downloads/clean.csv")

# Set the names used for saving results and plots
fileNameRoot = "LifeExp-MhierTrend-"
graphFileType = "png"

# Run the Bayesian model using three MCMC chains
# Save 15000 steps and keep every second draw
# Parallel processing is used if runjags is available
mcmcCoda = genMCMC(
  data = dat,
  yName = "life_expectancy_0",
  sName = "entity",
  cName = "continent",
  xName = "year",
  numSavedSteps = 15000,
  thinSteps = 2,
  useRunjags = TRUE,
  saveName = fileNameRoot
)

# Display how long the model took to run
cat("MCMC run time (seconds):", attr(mcmcCoda, "seconds"), "\n")

# Check the convergence and mixing of selected parameters
# These include continent trends, curvature, the COVID effect
# and the overall standard deviation
parNames = varnames(mcmcCoda)
for (parName in c("beta[2]", "g[6]", "q[6]", "q[1]", "sigma")) {
  diagMCMC(
    codaObject = mcmcCoda,
    parName = parName,
    saveName = fileNameRoot,
    saveType = graphFileType
  )
}

# Calculate posterior summary statistics
# These include means, medians, credible intervals and ESS
sumInfo = smryMCMC(
  mcmcCoda,
  compVal = 0,
  saveName = fileNameRoot
)

# Display the summary results
show(sumInfo)

# Plot the posterior distributions and compare continent trends
plotMCMC(
  mcmcCoda,
  data = dat,
  compVal = 0,
  saveName = fileNameRoot,
  saveType = graphFileType
)
