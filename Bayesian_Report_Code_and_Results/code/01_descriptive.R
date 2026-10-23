# 01_descriptive.R - descriptive statistics, covariance/correlation, figures
source("00_common.R")

d <- load_data()
dir.create("out", showWarnings = FALSE); dir.create("fig", showWarnings = FALSE)

th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11.5),
        legend.position = "bottom", strip.text = element_text(face = "bold"))
colscale <- scale_color_manual(values = cols); fillscale <- scale_fill_manual(values = cols)
savefig <- function(p, name, w = 6.7, h = 4) ggsave(file.path("fig", name), p, width = w, height = h, dpi = 300, bg = "white")

# summary statistics
summ <- function(x) c(n = length(x), mean = mean(x), sd = sd(x), min = min(x),
                      q1 = unname(quantile(x, .25)), median = median(x),
                      q3 = unname(quantile(x, .75)), max = max(x), iqr = IQR(x),
                      skew = skewness(x, type = 2), exkurt = kurtosis(x, type = 2))
tab_all  <- data.frame(group = "All country-years", t(summ(d$life_expectancy_0)))
tab_cont <- do.call(rbind, lapply(conts, function(k) {
  x <- d[d$continent == k, ]
  data.frame(group = k, countries = length(unique(x$entity)), t(summ(x$life_expectancy_0)))}))
tab_year <- do.call(rbind, lapply(c(1950, 1970, 1990, 2010, 2019, 2021, 2023), function(y) {
  x <- d$life_expectancy_0[d$year == y]; data.frame(year = y, t(summ(x)))}))
write.csv(tab_all,  "out/desc_all.csv",  row.names = FALSE)
write.csv(tab_cont, "out/desc_continent.csv", row.names = FALSE)
write.csv(tab_year, "out/desc_year.csv", row.names = FALSE)

# raw 2019 -> 2021 change by continent (descriptive evidence for the COVID effect)
cv <- do.call(rbind, lapply(conts, function(k) {
  m19 <- mean(d$life_expectancy_0[d$continent == k & d$year == 2019])
  m21 <- mean(d$life_expectancy_0[d$continent == k & d$year == 2021])
  data.frame(continent = k, y2019 = m19, y2021 = m21, change = m21 - m19)}))
cv <- rbind(cv, data.frame(continent = "All countries",
  y2019 = mean(d$life_expectancy_0[d$year == 2019]), y2021 = mean(d$life_expectancy_0[d$year == 2021]),
  change = mean(d$life_expectancy_0[d$year == 2021]) - mean(d$life_expectancy_0[d$year == 2019])))
write.csv(cv, "out/desc_covid_raw.csv", row.names = FALSE)

# predictors: covariance and correlation 
X <- data.frame(t = d$t, t2 = d$t2, covid = d$covid,
                Asia = as.integer(d$continent == "Asia"), Europe = as.integer(d$continent == "Europe"),
                NorthAmerica = as.integer(d$continent == "North America"),
                Oceania = as.integer(d$continent == "Oceania"),
                SouthAmerica = as.integer(d$continent == "South America"))
write.csv(round(cov(X), 4), "out/cov_predictors.csv")
write.csv(round(cor(X), 4), "out/cor_predictors.csv")
write.csv(data.frame(predictor = names(X), cor_with_y = round(sapply(X, function(z) cor(z, d$life_expectancy_0)), 4)),
          "out/cor_with_y.csv", row.names = FALSE)
pred_sum <- data.frame(variable = names(X), mean = sapply(X, mean), sd = sapply(X, sd),
                       min = sapply(X, min), max = sapply(X, max))
write.csv(pred_sum, "out/desc_predictors.csv", row.names = FALSE)

# country-level gain vs starting level 
w <- reshape(d[d$year %in% c(1950, 2023), c("entity", "continent", "year", "life_expectancy_0")],
             idvar = c("entity", "continent"), timevar = "year", direction = "wide")
names(w)[3:4] <- c("y1950", "y2023"); w$gain <- w$y2023 - w$y1950
fit_gain <- lm(gain ~ y1950, data = w)
write.csv(data.frame(r = cor(w$y1950, w$gain), slope = coef(fit_gain)[2], intercept = coef(fit_gain)[1],
                     r2 = summary(fit_gain)$r.squared, n = nrow(w)), "out/gain_vs_start.csv", row.names = FALSE)
write.csv(w, "out/gain_by_country.csv", row.names = FALSE)

#  figures
d$continent <- factor(d$continent, levels = conts)
w$continent <- factor(w$continent, levels = conts)

# Fig 1: distribution of life expectancy
p1 <- ggplot(d, aes(life_expectancy_0)) + geom_histogram(binwidth = 2, fill = "#2a78d6", colour = "white", linewidth = .2) +
  geom_vline(xintercept = mean(d$life_expectancy_0), colour = "#eb6834", linewidth = .9, linetype = 2) +
  labs(title = "(a) All country-years", x = "Life expectancy at birth (years)", y = "Count") + th
p2 <- ggplot(d, aes(life_expectancy_0, colour = continent)) + geom_density(linewidth = 1.1) + colscale +
  labs(title = "(b) By continent", x = "Life expectancy at birth (years)", y = "Density", colour = NULL) + th +
  theme(legend.text = element_text(size = 8.5), legend.key.size = unit(0.35, "cm")) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(linewidth = 1.5)))
savefig(p1 + p2, "fig01_distribution.png", 6.7, 3.9)

# Fig 2: boxplots 1950 vs 2023
d2 <- d[d$year %in% c(1950, 2023), ]
p <- ggplot(d2, aes(continent, life_expectancy_0, fill = continent)) + geom_boxplot(outlier.size = .8, alpha = .85, linewidth = .4) +
  facet_wrap(~year) + fillscale + labs(x = NULL, y = "Life expectancy at birth (years)") + th +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
savefig(p, "fig02_boxplots.png", 6.7, 4.4)

# Fig 3: country trajectories by continent
mn <- aggregate(life_expectancy_0 ~ year + continent, d, mean)
p <- ggplot() + annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf, fill = "#eda100", alpha = .25) +
  geom_line(data = d, aes(year, life_expectancy_0, group = entity), colour = "grey55", alpha = .35, linewidth = .3) +
  geom_line(data = mn, aes(year, life_expectancy_0, colour = continent), linewidth = 1.4) + colscale +
  facet_wrap(~continent, ncol = 3) + labs(x = "Year", y = "Life expectancy at birth (years)") + th + theme(legend.position = "none")
savefig(p, "fig03_trajectories.png", 6.7, 5.4)

# Fig 4: continent means with COVID window
wm <- aggregate(life_expectancy_0 ~ year, d, mean)
lab <- mn[mn$year == 2023, ]
# spread the end-of-line labels so they never overlap (minimum gap 1.5 years)
lab <- lab[order(lab$life_expectancy_0), ]; lab$ylab <- lab$life_expectancy_0
for (i in 2:nrow(lab)) lab$ylab[i] <- max(lab$ylab[i], lab$ylab[i - 1] + 1.5)
lab$ylab <- lab$ylab - (max(lab$ylab) - max(lab$life_expectancy_0)) / 2
p <- ggplot() + annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf, fill = "#eda100", alpha = .25) +
  geom_line(data = mn, aes(year, life_expectancy_0, colour = continent), linewidth = 1.2) +
  geom_line(data = wm, aes(year, life_expectancy_0), colour = "black", linetype = 2, linewidth = .9) +
  geom_segment(data = lab, aes(x = 2023.2, xend = 2024.6, y = life_expectancy_0, yend = ylab, colour = continent), linewidth = .5) +
  geom_text(data = lab, aes(2025, ylab, label = continent, colour = continent), hjust = 0, size = 3.3, fontface = "bold") +
  annotate("text", x = 2025, y = wm$life_expectancy_0[wm$year == 2023] - 3.6, label = "All countries (dashed)", hjust = 0, size = 3.3) +
  colscale + scale_x_continuous(limits = c(1950, 2046), breaks = seq(1950, 2020, 10)) +
  labs(x = "Year", y = "Mean life expectancy at birth (years)") + th + theme(legend.position = "none")
savefig(p, "fig04_continent_means.png", 6.7, 4.2)

# Fig 5: correlation heatmap of predictors
cm <- cor(X); cmm <- reshape2::melt(cm)
lbl <- c(t = "t (year)", t2 = "t squared", covid = "COVID 2020-21", Asia = "Asia", Europe = "Europe",
         NorthAmerica = "N. America", Oceania = "Oceania", SouthAmerica = "S. America")
cmm$Var1 <- factor(lbl[as.character(cmm$Var1)], levels = rev(lbl)); cmm$Var2 <- factor(lbl[as.character(cmm$Var2)], levels = lbl)
p <- ggplot(cmm, aes(Var2, Var1, fill = value)) + geom_tile(colour = "white") +
  geom_text(aes(label = sub("^-0\\.00$", "0.00", sprintf("%.2f", value))), size = 3.4) +
  scale_fill_gradient2(low = "#eb6834", mid = "white", high = "#2a78d6", limits = c(-1, 1), name = "Correlation") +
  labs(x = NULL, y = NULL) + th + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
savefig(p, "fig05_correlation.png", 6.7, 5.0)

# Fig 6: gain vs starting level
p <- ggplot(w, aes(y1950, gain, colour = continent)) + geom_point(size = 1.8, alpha = .8) + colscale +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, colour = "black", linewidth = .8, fill = "grey80") +
  labs(x = "Life expectancy in 1950 (years)", y = "Gain from 1950 to 2023 (years)", colour = NULL) + th +
  guides(colour = guide_legend(nrow = 1))
savefig(p, "fig06_gain_vs_start.png", 6.7, 4.2)

print(tab_all); print(tab_cont); print(tab_year); print(cv); print(round(cor(X), 3))
print(read.csv("out/cor_with_y.csv")); print(read.csv("out/gain_vs_start.csv"))
