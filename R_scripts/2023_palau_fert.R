# 2023 palau fert


# notes -------------------------------------------------------------------
# - gamm might be slightly better fitted vs resid 
# outlier at distnace = 30, deg = 25 is likely from two spokes crossing, check genetics




# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
library(glmmTMB)
library(lubridate)
library(gamm4)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek3") # set theme in code


# 1 Import data -----------------------------------------------------------
# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data1 <- read.excel() # read clipboard from excel
# save(data1, file = file.path("./Rdata", "2023_palau_fert.RData"))
load("./Rdata/2023_palau_fert.RData") #data1

# 2. Labelling and wrangling -----------------------------------------------
str(data1) # check data type is correct
data1$spoke <- as.factor(as.character(data1$spoke))
data1$id <- as.factor(as.character(data1$id))
data1$obs <- factor(formatC(1:nrow(data1), flag = "0", width = 3)) # unique tank ID for later on
data1$prop <- data1$suc / data1$tot
data1$time <- sprintf("%s PM", data1$time)  # Append PM to each time string
data1$time <- as.POSIXct(paste("2000-01-01", data1$time), format="%Y-%m-%d %I:%M %p", tz="UTC")
base_time <- min(data1$time, na.rm = T)
data1$time_from_base <- as.numeric(difftime(data1$time, base_time, units = "secs"))
data1$quality_score = 1
data1 <- data1[!is.na(data1$prop), ]
data1$suc = ifelse(data1$suc == 0, 1, data1$suc)  #try adding small error to work betabinom

# adjust the quality_score based on specific conditions in the dataframe (not on atm)
data1$quality_score[data1$id == "5_05" & data1$spoke == 5] <- 0.5  # Set quality_score to 50 for id 5_05 and spoke 5
data1$quality_score[data1$id == "7_05" & data1$spoke == 7] <- 0.5  # Set quality_score to 50 for id 7_05 and spoke 7


# 3. Data Exploration ----------------------------------------------------

# Wrangling
centre_indiv <- data1[which(data1$spoke == 'c'),]
radial_indiv <- data1[which(data1$spoke != 'c'),]

# Calculate degrees from a baseline (s1) and adjust by subtracting an offset (s5)
# radial_indiv <- radial_indiv %>%
#   mutate(deg_from_s1 = case_when(spoke == 2 ~ 26,spoke == 3 ~ 52,spoke == 4 ~ 78,spoke == 5 ~ 104,spoke == 6 ~ 130,spoke == 7 ~ 156,TRUE ~ NA_real_  # Assign NA for any other values
#   ))

# Calculate degrees from north
radial_indiv <- radial_indiv %>% mutate(deg_from_north = case_when(spoke == 2 ~ 77, spoke == 3 ~ 101, spoke == 4 ~ 124, 
                                                                   spoke == 5 ~ 148, spoke == 6 ~ 173, spoke == 7 ~ 201, TRUE ~ NA_real_))  # Assign NA for any other values
#taken from earth 
#(deg_from_s5 <- deg_from_s1 - 104) 
diff(unique(radial_indiv$deg_from_north))

radial_indiv$deg = radial_indiv$deg_from_north - 146.5   #adjust to water dir
range(radial_indiv$deg )
radial_indiv$deg_rad <- radial_indiv$deg * pi / 180
radial_indiv$sin_deg <- sin(radial_indiv$deg_rad)
radial_indiv$cos_deg <- cos(radial_indiv$deg_rad)
rad_lines = unique(radial_indiv$deg)

## overall mean and weighted means
# sum of all embryos
(with(data1, sum(suc, na.rm = T) / sum(tot, na.rm = T)))  #0.235
#mean of all samples
with(data1, mean(prop, na.rm = T))  #unweighted
#weighted mean
(wei_mean_prop_cent <- sum(data1$prop * data1$tot, na.rm = TRUE) / sum(data1$tot, na.rm = TRUE))
(wei_mean_prop_cent <- sum(centre_indiv$prop * centre_indiv$tot, na.rm = TRUE) / sum(centre_indiv$tot, na.rm = TRUE))  #centre
(wei_mean_prop_rad <- sum(radial_indiv$prop * radial_indiv$tot, na.rm = TRUE) / sum(radial_indiv$tot, na.rm = TRUE))  #spokes


# Modelling ------------------------------------------------------------

## glmm
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
#devtools::install_version("glmmTMB", version = "1.1.08", repos = "http://cran.us.r-project.org")

# str(radial_indiv)
# radial_indiv$cos_deg
# levels(radial_indiv$obs)
# md1 <- glmmTMB(cbind(suc, (tot - suc)) ~ scale(dist) + poly(deg, 2) + (1 | obs), family = "betabinomial", data = radial_indiv)
# summary(md1)
# AIC(md1)
# plot(fitted(md1), resid(md1)) # fitted vs residuals
# abline(h = 0)
# ## betabiomial better fitted vs resid
# 
# 
# ## technically angles should be in radians (sin and cos) or dealt with a gam (see below0
# # library(mgcv)
# # library(gamm4)
# # md1 <- gamm4(cbind(suc, tot - suc) ~ s(deg, bs = "cc") + dist + (1 | obs), family = binomial, data = radial_indiv)
# # md1 <- glmmTMB(cbind(suc, tot - suc) ~ 
# #                  s(deg, bs = "cc", k = 4) + 
# #                  s(dist, k = 4) + 
# #                  (1 | obs), 
# #                family = betabinomial(),
# #                data = radial_indiv)
# 
# 
# # radian
# # angles_deg <- seq(min(radial_indiv$deg), max(radial_indiv$deg), length.out = 100)  # -69.5 to 54.5 degrees
# # angles_rad <- angles_deg * (pi / 180)  # Convert to radians
# # sin_deg <- sin(angles_rad)
# # cos_deg <- cos(angles_rad)
# 
# # Generate new_data for predictions
# new_data <- expand.grid(
#   dist = seq(min(radial_indiv$dist), max(radial_indiv$dist), length.out = 100),
#   deg = seq(min(radial_indiv$deg), max(radial_indiv$deg), length.out = 100),
#   # sin_deg = sin_deg,
#   # cos_deg = cos_deg,
#   quality_score = 1
# )
# 
# #new_data$deg <- angles_deg
# range(new_data$deg)  # Should match radial_indiv$deg range
# new_data$spoke <- factor("s5") # Example, choose appropriately based on your model structure
# new_data$obs <- factor("001") # Assuming a single observation for prediction purposes
# new_data$predicted <- predict(md1, newdata = new_data, type = "response")
# 
# p0 <- ggplot(new_data, aes(x = dist, y = deg)) +
#   geom_raster(aes(fill = predicted)) + # Use geom_raster for a continuous field of colors
#   scale_fill_gradient(low = "#3B9AB2", high = "#F21A00") + # Color gradient for predictions
#   labs(x = "Distance from centre patch (m)", y = "Degrees from downstream", fill = "Fertilisation \n success") +
#   geom_contour(aes(z = predicted), breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), color = "white") + # Draw contour lines based on 'predicted'
#   theme_minimal()
# p0 <- p0 + geom_hline(yintercept = rad_lines, color = "#E1AF00", lty = 2)
# p0


# GAMM --------------------------------------------------------------------

#try gamm
md1 <- gamm(cbind(suc, tot - suc) ~ s(deg, k = 3) + dist,   random = list(obs = ~1), family = binomial,  method = "REML", verbosePQL = F, 
            data = radial_indiv)
md1$gam
coef(md1)
summary(md1)
AIC(md1)  #k = 3 best
plot(fitted(md1$gam), residuals(md1$gam)) # fitted vs residuals
abline(h = 0)
# Generate new_data for predictions
new_data <- expand.grid(
  dist = seq(min(radial_indiv$dist), max(radial_indiv$dist), length.out = 100),
  deg = seq(min(radial_indiv$deg), max(radial_indiv$deg), length.out = 100),
  quality_score = 1  # Assuming this is required for consistency, though not used in the model
)

hist(new_data$predicted)
quantile(new_data$predicted)

# Add a dummy observation level for prediction purposes
new_data$obs <- factor("001")  # Assuming a single observation for prediction
new_data$predicted <- predict(md1$gam, newdata = new_data, type = "response")
radial_indiv$residuals <- residuals(md1$gam, type = "deviance")

p0 <- ggplot(new_data, aes(x = dist, y = deg)) +
  geom_raster(aes(fill = predicted)) +
  scale_fill_gradient(low = "#3B9AB2", high = "#F21A00", breaks = seq(0.1, 0.9, by = 0.1)) +
  labs(x = "Distance from centre patch (m)", y = "Downstream of centre patch (degrees)", fill = "Fertilisation \n success") +
  geom_contour(aes(z = predicted), breaks = seq(0.1, 0.9, by = 0.1), color = "white") +
  theme_minimal() +
  scale_y_reverse() + 
  geom_hline(yintercept = rad_lines, color = "#E1AF00", lty = 2) +
  theme(
    legend.key.height = unit(1.5, "cm"), # Increase the height of legend keys
    axis.title = element_text(size = 13), # Increase axis label size
    axis.text = element_text(size = 12) # Increase tick label size
  )# add deviance
# p0 = p0 +  geom_point(data = radial_indiv, aes(x = dist, y = deg, color = residuals), size = 2)+
#   scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
p0

ggsave(p0, filename = '2023palau_fert.tiff',  path = "./plots", device = "tiff",  width = 8, height = 5)  #this often works better than pdf
