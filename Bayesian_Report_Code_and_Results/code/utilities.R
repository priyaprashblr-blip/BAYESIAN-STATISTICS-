
# #Utility functions for Doing Bayesian Data Analysis
# John K. Kruschke (2015)



# Load packages that are needed for the functions
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Please install the coda package first.")
}

# Open a graphics window on Mac, Windows or Linux
openGraph <- function(width = 7, height = 7, mag = 1, ...) {
  if (Sys.info()[["sysname"]] == "Darwin") {
    grDevices::quartz(width = width * mag, height = height * mag, ...)
  } else if (.Platform$OS.type == "windows") {
    grDevices::windows(width = width * mag, height = height * mag, ...)
  } else {
    grDevices::X11(width = width * mag, height = height * mag,
                   type = "cairo", ...)
  }
}

# Save the current graph
saveGraph <- function(file = "saveGraphOutput", type = "pdf", ...) {
  type <- tolower(type)
  filename <- paste0(file, ".", type)

  if (type == "pdf") {
    grDevices::dev.copy2pdf(file = filename, ...)
  } else if (type == "eps") {
    grDevices::dev.copy2eps(file = filename, ...)
  } else if (type %in% c("png", "jpeg", "jpg", "tiff", "bmp")) {
    device <- switch(type,
      png = grDevices::png,
      jpeg = grDevices::jpeg,
      jpg = grDevices::jpeg,
      tiff = grDevices::tiff,
      bmp = grDevices::bmp
    )
    grDevices::dev.copy(device, filename = filename, ...)
    grDevices::dev.off()
  } else {
    stop("Unsupported file type.")
  }
}

# Calculate the highest density interval from MCMC samples
HDIofMCMC <- function(sampleVec, credMass = 0.95) {
  sampleVec <- sort(sampleVec[is.finite(sampleVec)])

  if (length(sampleVec) < 2) {
    stop("At least two finite samples are required.")
  }
  if (credMass <= 0 || credMass >= 1) {
    stop("credMass must be between 0 and 1.")
  }

  n <- length(sampleVec)
  intervalSize <- max(1, ceiling(credMass * n))
  lastStart <- n - intervalSize + 1

  if (lastStart < 1) {
    return(range(sampleVec))
  }

  widths <- sampleVec[seq_len(lastStart) + intervalSize - 1] -
    sampleVec[seq_len(lastStart)]

  best <- which.min(widths)
  c(sampleVec[best], sampleVec[best + intervalSize - 1])
}

# Calculate an HDI from an inverse cumulative distribution function
HDIofICDF <- function(ICDFname, credMass = 0.95, tol = 1e-8, ...) {
  if (credMass <= 0 || credMass >= 1) {
    stop("credMass must be between 0 and 1.")
  }

  tailMass <- 1 - credMass

  intervalWidth <- function(lowTailPr, ICDFname, credMass, ...) {
    ICDFname(credMass + lowTailPr, ...) -
      ICDFname(lowTailPr, ...)
  }

  result <- optimize(
    intervalWidth,
    interval = c(0, tailMass),
    ICDFname = ICDFname,
    credMass = credMass,
    tol = tol,
    ...
  )

  c(
    ICDFname(result$minimum, ...),
    ICDFname(credMass + result$minimum, ...)
  )
}

# Find the HDI on a probability grid
HDIofGrid <- function(probMassVec, credMass = 0.95) {
  if (any(!is.finite(probMassVec)) || any(probMassVec < 0) ||
      sum(probMassVec) <= 0) {
    stop("Probability masses must be finite, non-negative and non-zero.")
  }

  probMassVec <- probMassVec / sum(probMassVec)
  sortedMass <- sort(probMassVec, decreasing = TRUE)
  cutoff <- sortedMass[which(cumsum(sortedMass) >= credMass)[1]]

  list(
    indices = which(probMassVec >= cutoff),
    mass = sum(probMassVec[probMassVec >= cutoff]),
    height = cutoff
  )
}

# Plot autocorrelation for MCMC chains
DbdaAcfPlot <- function(codaObject,
                        parName = coda::varnames(codaObject)[1],
                        plColors = NULL) {
  if (!parName %in% coda::varnames(codaObject)) {
    stop("parName must be a column name of the coda object.")
  }

  nChain <- length(codaObject)
  if (is.null(plColors)) plColors <- seq_len(nChain)

  acfList <- lapply(seq_len(nChain), function(i) {
    stats::acf(
      as.numeric(codaObject[[i]][, parName]),
      plot = FALSE
    )
  })

  maxLag <- min(vapply(acfList, function(x) length(x$acf) - 1,
                       numeric(1)))
  xMat <- sapply(acfList, function(x) as.numeric(x$lag)[1:(maxLag + 1)])
  yMat <- sapply(acfList, function(x) as.numeric(x$acf)[1:(maxLag + 1)])

  graphics::matplot(xMat, yMat, type = "o", pch = 20,
                    col = plColors, ylim = c(-1, 1),
                    xlab = "Lag", ylab = "Autocorrelation")
  graphics::abline(h = 0, lty = "dashed")
  graphics::text(max(xMat), max(yMat),
                 labels = paste("ESS =", round(coda::effectiveSize(
                   codaObject[, parName]), 1)),
                 adj = c(1, 1))
}

# Plot posterior densities for MCMC chains
DbdaDensPlot <- function(codaObject,
                         parName = coda::varnames(codaObject)[1],
                         plColors = NULL) {
  if (!parName %in% coda::varnames(codaObject)) {
    stop("parName must be a column name of the coda object.")
  }

  nChain <- length(codaObject)
  if (is.null(plColors)) plColors <- seq_len(nChain)

  densList <- lapply(seq_len(nChain), function(i) {
    stats::density(as.numeric(codaObject[[i]][, parName]))
  })

  graphics::matplot(
    do.call(cbind, lapply(densList, `[[`, "x")),
    do.call(cbind, lapply(densList, `[[`, "y")),
    type = "l", col = plColors,
    xlab = "Parameter value", ylab = "Density"
  )

  hdi <- lapply(seq_len(nChain), function(i) {
    HDIofMCMC(as.numeric(codaObject[[i]][, parName]))
  })

  for (i in seq_len(nChain)) {
    graphics::points(hdi[[i]], c(0, 0),
                     col = plColors[i], pch = "|", cex = 2)
  }

  graphics::abline(h = 0)
}

# Show diagnostic plots for an MCMC parameter
diagMCMC <- function(codaObject,
                     parName = coda::varnames(codaObject)[1],
                     saveName = NULL, saveType = "jpg") {
  colors <- c("skyblue", "black", "royalblue", "steelblue")

  openGraph(height = 7, width = 9)
  oldPar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldPar), add = TRUE)

  graphics::par(mfrow = c(2, 2),
                mar = c(4, 4, 2, 1))

  coda::traceplot(codaObject[, parName], main = "",
                  col = colors, ylab = "Parameter value")

  tryCatch(
    coda::gelman.plot(codaObject[, parName],
                      main = "", auto.layout = FALSE,
                      col = colors),
    error = function(e) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Gelman plot unavailable")
    }
  )

  DbdaAcfPlot(codaObject, parName, colors)
  DbdaDensPlot(codaObject, parName, colors)
  graphics::mtext(parName, outer = TRUE, line = -1)

  if (!is.null(saveName)) {
    saveGraph(paste0(saveName, "Diag", parName), saveType)
  }
}

# Convert a Stan fit to coda format and show diagnostic plots
diagStanFit <- function(stanFit, parName,
                        saveName = NULL, saveType = "jpg") {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Please install the rstan package first.")
  }

  draws <- rstan::extract(stanFit, pars = parName,
                          permuted = FALSE)

  if (length(dim(draws)) != 3) {
    stop("Could not extract the selected Stan parameter.")
  }

  codaFit <- coda::mcmc.list(lapply(seq_len(dim(draws)[2]), function(i) {
    coda::mcmc(matrix(draws[, i, 1],
                      ncol = 1,
                      dimnames = list(NULL, parName)))
  }))

  diagMCMC(codaFit, parName, saveName, saveType)
}

# Normalize a vector so that its values sum to one
normalize <- function(v) {
  if (any(!is.finite(v)) || sum(v) == 0) {
    stop("Values must be finite and their sum must not be zero.")
  }
  v / sum(v)
}

# Summarize posterior samples
summarizePost <- function(paramSampleVec, compVal = NULL,
                          ROPE = NULL, credMass = 0.95) {
  x <- as.numeric(paramSampleVec)
  x <- x[is.finite(x)]

  if (length(x) < 2) stop("At least two finite samples are required.")

  dens <- stats::density(x)
  modeVal <- dens$x[which.max(dens$y)]
  hdi <- HDIofMCMC(x, credMass)

  if (!is.null(compVal)) {
    pcgt <- 100 * mean(x > compVal)
  } else {
    compVal <- NA
    pcgt <- NA
  }

  if (!is.null(ROPE)) {
    if (length(ROPE) != 2 || ROPE[1] > ROPE[2]) {
      stop("ROPE must contain a lower and upper limit.")
    }
    pBelow <- 100 * mean(x < ROPE[1])
    pAbove <- 100 * mean(x > ROPE[2])
    pInside <- 100 - pBelow - pAbove
  } else {
    ROPE <- c(NA, NA)
    pBelow <- pInside <- pAbove <- NA
  }

  c(
    Mean = mean(x),
    Median = median(x),
    Mode = modeVal,
    ESS = round(coda::effectiveSize(x), 1),
    HDImass = credMass,
    HDIlow = hdi[1],
    HDIhigh = hdi[2],
    CompVal = compVal,
    PcntGtCompVal = pcgt,
    ROPElow = ROPE[1],
    ROPEhigh = ROPE[2],
    PcntLtROPE = pBelow,
    PcntInROPE = pInside,
    PcntGtROPE = pAbove
  )
}

# Plot posterior samples and display the HDI
plotPost <- function(paramSampleVec,
                     cenTend = c("mode", "median", "mean")[1],
                     compVal = NULL, ROPE = NULL,
                     credMass = 0.95, HDItextPlace = 0.7,
                     xlab = "Parameter value", xlim = NULL,
                     yaxt = "n", ylab = "",
                     main = "", cex = 1.2, cex.lab = 1.3,
                     col = "skyblue", border = "white",
                     showCurve = FALSE, breaks = NULL, ...) {

  x <- as.numeric(paramSampleVec)
  x <- x[is.finite(x)]

  if (length(x) < 2) stop("At least two finite samples are required.")

  hdi <- HDIofMCMC(x, credMass)
  if (is.null(xlim)) xlim <- range(c(x, compVal, ROPE), na.rm = TRUE)

  if (showCurve) {
    histInfo <- graphics::hist(x, plot = FALSE, breaks = breaks)
    dens <- stats::density(x, adjust = 2)
    graphics::plot(dens$x, dens$y, type = "l", lwd = 3,
                   col = col, xlim = xlim, xlab = xlab,
                   ylab = ylab, main = main, ...)
  } else {
    histInfo <- graphics::hist(
      x, probability = TRUE, breaks = breaks,
      col = col, border = border, xlim = xlim,
      xlab = xlab, ylab = ylab, main = main,
      cex = cex, cex.lab = cex.lab, yaxt = yaxt, ...
    )
  }

  peak <- max(histInfo$density)
  mn <- mean(x)
  med <- median(x)
  dens <- stats::density(x)
  mo <- dens$x[which.max(dens$y)]

  center <- switch(cenTend,
                   mode = mo,
                   median = med,
                   mean = mn)

  graphics::abline(v = center, col = "darkgreen", lty = 2)
  graphics::text(center, peak * 0.9,
                 labels = paste(cenTend, "=", round(center, 3)))

  if (!is.null(compVal)) {
    graphics::abline(v = compVal, col = "darkgreen", lty = 3)
    graphics::text(compVal, peak * 0.7,
                   labels = paste(round(100 * mean(x > compVal), 1),
                                  "% above comparison value"),
                   cex = 0.75)
  }

  if (!is.null(ROPE)) {
    graphics::abline(v = ROPE, col = "darkred", lty = 3)
    graphics::abline(v = ROPE, col = "darkred", lty = 3)
  }

  graphics::segments(hdi[1], 0, hdi[2], 0, lwd = 5)
  graphics::text(mean(hdi), 0,
                 labels = paste(100 * credMass, "% HDI"),
                 pos = 3)

  invisible(summarizePost(x, compVal, ROPE, credMass))
}

# Calculate beta distribution parameters from the mean and concentration
betaABfromMeanKappa <- function(mean, kappa) {
  if (mean <= 0 || mean >= 1) stop("Mean must be between 0 and 1.")
  if (kappa <= 0) stop("Kappa must be positive.")

  list(a = mean * kappa, b = (1 - mean) * kappa)
}

# Calculate beta distribution parameters from the mode and concentration
betaABfromModeKappa <- function(mode, kappa) {
  if (mode <= 0 || mode >= 1) stop("Mode must be between 0 and 1.")
  if (kappa <= 2) stop("Kappa must be greater than 2.")

  list(
    a = mode * (kappa - 2) + 1,
    b = (1 - mode) * (kappa - 2) + 1
  )
}

# Calculate beta distribution parameters from the mean and standard deviation
betaABfromMeanSD <- function(mean, sd) {
  if (mean <= 0 || mean >= 1) stop("Mean must be between 0 and 1.")
  if (sd <= 0) stop("Standard deviation must be positive.")

  kappa <- mean * (1 - mean) / sd^2 - 1
  if (kappa <= 0) stop("Invalid combination of mean and standard deviation.")

  list(a = mean * kappa, b = (1 - mean) * kappa)
}

# Calculate gamma distribution parameters from the mean and standard deviation
gammaShRaFromMeanSD <- function(mean, sd) {
  if (mean <= 0 || sd <= 0) stop("Mean and standard deviation must be positive.")

  list(shape = mean^2 / sd^2, rate = mean / sd^2)
}

# Calculate gamma distribution parameters from the mode and standard deviation
gammaShRaFromModeSD <- function(mode, sd) {
  if (mode <= 0 || sd <= 0) stop("Mode and standard deviation must be positive.")

  rate <- (mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2)
  shape <- 1 + mode * rate

  list(shape = shape, rate = rate)
}
