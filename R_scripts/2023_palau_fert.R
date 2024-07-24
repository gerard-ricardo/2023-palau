# 2023 palau fert

# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code
library(lubridate)


# 1 Import data -----------------------------------------------------------
# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data1 <- read.excel() # read clipboard from excel
# save(data1, file = file.path("./Rdata", "2023_palau_fert.RData"))
load("./Rdata/2023_palau_fert.RData")



# 2. Labelling and wrangling -----------------------------------------------
str(data1) # check data type is correct
data1$spoke <- as.factor(as.character(data1$spoke))
data1$id <- as.factor(as.character(data1$id))
data1$obs <- factor(formatC(1:nrow(data1), flag = "0", width = 3)) # unique tank ID for later on
data1$prop <- data1$suc / data1$tot
# Convert to POSIXct, assuming all times are PM
data1$time <- sprintf("%s PM", data1$time)  # Append PM to each time string
data1$time <- as.POSIXct(paste("2000-01-01", data1$time), format="%Y-%m-%d %I:%M %p", tz="UTC")
base_time <- min(data1$time, na.rm = T)
data1$time_from_base <- as.numeric(difftime(data1$time, base_time, units = "secs"))


# 3. Data Exploration ----------------------------------------------------
#





# sum of all embryos
(with(data1, sum(suc, na.rm = T) / sum(tot, na.rm = T)))  #22.4%
#mean of all samples
with(data1, mean(prop, na.rm = T))  #unweighted

# Wrangling
data2 <- data1[25:28, ] # only centre
data1 <- data1[1:24, ] # only radial lines

# Calculate degrees from a baseline (s1) and adjust by subtracting an offset (s5)
deg_from_s1 <- c(26, 52, 78, 104, 130, 156)
(deg_from_s5 <- deg_from_s1 - 104)
data1$deg <- sort(rep(deg_from_s5, 4))

# #quick data input (need to check)
# x  = c(0.7, 3, 10, 30)
# deg_from_s1 = c(26, 52, 78, 104, 130, 156)
# (deg_from_s5 = deg_from_s1 - 104)
# spoke2 = c(0,1,2,NA)
# spoke3 = c(9,1,0,0)
# spoke4 = c(75, 26, 30, 13)
# spoke5 = c(74, NA, 84, 0)
# spoke6 = c(92, NA, 20, 65)
# spoke7 = c(90, 54, 9, 2)
# all_spokes = c(spoke2, spoke3, spoke4, spoke5, spoke6, spoke7)
# length(dist)
# data1 = data.frame(suc = all_spokes, tot = 100, dist = rep(x, 6), deg = sort(rep(deg_from_s5, 4)), spoke = sort(rep(c('s2', 's3', 's4', 's5', 's6', 's7'), 4)))


# 4. Modelling ------------------------------------------------------------
# library(lme4)
# md1 <- glmer(cbind(suc, (tot - suc)) ~ scale(dist) * scale(deg) + (1 | obs), family = binomial, data = data1)

# md1 <- glmer(cbind(suc, (tot - suc)) ~ (dist) * poly(deg, 2) + (1 | obs), family = binomial, data = data1)
# no interactive effect
# AIC(md1)
# md1 <- glmer(cbind(suc, (tot - suc)) ~ (dist) + poly(deg, 2) + (1 | obs), family = binomial, data = data1)
# AIC(md1)
# summary(md1)
# plot(fitted(md1), resid(md1)) # fitted vs residuals
# abline(h = 0)

library(glmmTMB)
md1 <- glmmTMB(cbind(suc, (tot - suc)) ~ (dist) + poly(deg, 2) + (1 | obs), family = "betabinomial", data = data1)
summary(md1)
AIC(md1)
plot(fitted(md1), resid(md1)) # fitted vs residuals
abline(h = 0)
## betabiomial better fitted vs resid

# mes

# 1. Generate a new data frame for predictions
new_data <- expand.grid(
  dist = seq(min(data1$dist), max(data1$dist), length.out = 100),
  deg = seq(min(data1$deg), max(data1$deg), length.out = 100)
)
new_data$spoke <- factor("s5") # Example, choose appropriately based on your model structure
new_data$obs <- factor("001") # Assuming a single observation for prediction purposes

# 2. Predict using the model
new_data$predicted <- predict(md1, newdata = new_data, type = "response")

# 3. Plotting
library(ggplot2)

# Assuming new_data is already created and contains the 'predicted' column
p0 <- ggplot(new_data, aes(x = dist, y = deg)) +
  geom_raster(aes(fill = predicted)) + # Use geom_raster for a continuous field of colors
  scale_fill_gradient(low = "#3B9AB2", high = "#F21A00") + # Color gradient for predictions
  labs(x = "Distance from centre patch (m)", y = "Degrees from downstream", fill = "Fertilisation \n success") +
  geom_contour(aes(z = predicted), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), color = "white") + # Draw contour lines based on 'predicted'
  theme_minimal()
deg_from_s5
p0 <- p0 + geom_hline(yintercept = c(-78, -52, -26, 0, 26, 52), color = "#E1AF00", lty = 2)
p0
#ggsave(p0, filename = '2023palau_fert.tiff',  path = "./plots", device = "tiff",  width = 8, height = 5)  #this often works better than pdf
