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
library(mgcv)
library(nlme)
library(dplyr)
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
hist(data1$suc / data1$tot)

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
nrow(data1)
#mean of all samples
with(data1, mean(prop, na.rm = T))  #unweighted
with(data1, sd(prop, na.rm = T))  #unweighted
#weighted mean
(wei_mean_prop_cent <- sum(data1$prop * data1$tot, na.rm = TRUE) / sum(data1$tot, na.rm = TRUE))
(wei_mean_prop_cent <- sum(centre_indiv$prop * centre_indiv$tot, na.rm = TRUE) / sum(centre_indiv$tot, na.rm = TRUE))  #centre
#I dont think weighting approporate here, as the low one is over weighted
print(paste0('unweighted centre fert: ' ,mean(centre_indiv$prop)))

#try downsampling
#central
set.seed(123) # Ensure reproducibility
target_n <- min(centre_indiv$tot)# Desired sample size per id
resampled_data <- centre_indiv %>% group_by(id) %>% mutate(prob_suc = suc / tot) %>%  filter(tot >= target_n) %>%  # Only resample when tot is large
  mutate(suc1 = rbinom(1, target_n, prob_suc), fail1 = target_n - suc1, tot1 = target_n, prop1 = suc1 / tot1) %>%
  ungroup()%>% data.frame()
mean(resampled_data$prop1)
sd(resampled_data$prop1)
print(paste0('unweighted downsampled centre fert: ' , 100*round(mean(resampled_data$prop1),4)))


#all
set.seed(123) # Ensure reproducibility
target_n_a <- 50# Desired sample size per id
resampled_data_all <- data1 %>% group_by(id) %>% mutate(prob_suc = suc / tot) %>%  filter(tot >= target_n_a) %>%  # Only resample when tot is large
  mutate(suc1 = rbinom(1, target_n, prob_suc), fail1 = target_n - suc1, tot1 = target_n, prop1 = suc1 / tot1) %>%
  ungroup() %>% data.frame()
mean(resampled_data_all$prop1)
print(paste0('unweighted downsampled centre fert: ' , 100 * round(mean(resampled_data_all$prop1), 4)))
median_prop <- median(resampled_data_all$prop1)
iqr_prop <- IQR(resampled_data_all$prop1)
print(paste0('Median fertilisation: ', 100*round(median_prop, 4), '%, IQR: ', 100*round(iqr_prop, 4), '%'))
range(resampled_data_all$prop1)


#radial
set.seed(123) # Ensure reproducibility
target_n_a <- 50 # Desired sample size per id
resampled_data_rad <- radial_indiv %>% group_by(id) %>% mutate(prob_suc = suc / tot) %>%  filter(tot >= target_n_a) %>%  # Only resample when tot is large
  mutate(suc1 = rbinom(1, target_n, prob_suc), fail1 = target_n - suc1, tot1 = target_n, prop1 = suc1 / tot1) %>%
  ungroup() %>% data.frame()
mean(resampled_data_rad$prop1)
range(resampled_data_rad$prop1)
print(paste0('unweighted downsampled centre fert: ' , 100 * round(mean(resampled_data_rad$prop1), 4))  )



(wei_mean_prop_rad <- sum(radial_indiv$prop * radial_indiv$tot, na.rm = TRUE) / sum(radial_indiv$tot, na.rm = TRUE))  #spokes


plot(radial_indiv$prop ~ radial_indiv$dist)
plot(radial_indiv$prop ~ radial_indiv$deg)
plot(radial_indiv$prop ~ radial_indiv$time_from_base)
text(radial_indiv$time_from_base, radial_indiv$prop, labels = radial_indiv$id, pos = 3, cex = 0.8, col = "red") # Add labels

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

# md1_int <- gamm(cbind(suc, tot - suc) ~ s(deg, k = 3) * dist,   random = list(obs = ~1), family = binomial,  method = "REML", verbosePQL = F, 
#             data = radial_indiv)  #error
md1$gam
coef(md1)
summary(md1)
summary(md1$gam)
AIC(md1)

# Extract coefficient for distance
coef_dist <- summary(md1$gam)$p.table["dist", "Estimate"] # -0.1168
pval_dist <- summary(md1$gam)$p.table["dist", "Pr(>|t|)"] # 0.00182
# Convert to odds ratio
odds_ratio <- exp(coef_dist) # ≈ 0.89
print(paste0("For each 1 m increase in distance, the odds of fertilisation decrease by ", 
             round((1 - odds_ratio) * 100, 1), "% (p = ", signif(pval_dist, 3), ")"))


# Extract fixed effects
fixef <- coef(md1$gam) # This gives intercept and dist coefficient
beta0 <- fixef["(Intercept)"]
beta1 <- fixef["dist"]

# Define a range of distances
dist_vals <- seq(0, 10, by = 1)

# Hold s(deg) = 0 for simplicity (i.e., centring effect)
eta <- beta0 + beta1 * dist_vals
prob <- exp(eta) / (1 + exp(eta)) # Inverse logit

# View results
data.frame(Distance = dist_vals, Probability = prob)

                   

############################
#test  slightly lower sperm conc for some samples (5_30, 3_30, 4_30 from mlg filter)
radial_indiv$clonemate_in_sperm <- ifelse(radial_indiv$egg_clone_in_sperm, 1, 0) # 1 if egg clone found among sperm donors
radial_indiv$effective_sperm_sources <- ifelse(radial_indiv$id %in% c("5_30", "3_30", "4_30"), 14, 15) # Decrease by 1 if clone in sperm
md2 <- gamm(cbind(suc, tot - suc) ~ s(deg, k = 3) + dist + effective_sperm_sources, 
              random = list(obs = ~1), family = binomial, method = "REML", data = radial_indiv)
summary(md2$gam)
AIC(md2) #didn't improve fit and impact negigable.


################
#k-fold CV
k <- 5 # Set number of folds
set.seed(123)  # Ensure reproducibility
# Create k-fold groups
radial_indiv$fold <- sample(rep(1:k, length.out = nrow(radial_indiv)))
# Function to run k-fold CV
cv_results <- lapply(1:k, function(i) {
  # Split data into training and testing
  train_data <- radial_indiv %>% filter(fold != i)
  test_data  <- radial_indiv %>% filter(fold == i)
  # Fit the GAMM model on training data
  model <- gamm(cbind(suc, tot - suc) ~ s(deg, k = 3) + dist, 
                random = list(obs = ~1), family = binomial, 
                method = "REML", data = train_data)
  # Predict on test data
  pred_probs <- predict(model$gam, newdata = test_data, type = "response")
  # Compute performance metric (log-likelihood)
  actual <- test_data$suc / test_data$tot  # Convert to observed proportion
  log_lik <- sum(actual * log(pred_probs) + (1 - actual) * log(1 - pred_probs), na.rm = TRUE)
  return(log_lik)
})
# Calculate mean log-likelihood across folds
(mean_log_lik <- mean(unlist(cv_results)))

#no diff between k = 2 and 3, >3 is worse
################

#AIC(md1)  #k = 3 best
plot(fitted(md1$gam), residuals(md1$gam)) # fitted vs residuals
abline(h = 0)
# Generate new_data for predictions
new_data <- expand.grid(
  dist = seq(min(radial_indiv$dist), max(radial_indiv$dist), length.out = 100),
  deg = seq(min(radial_indiv$deg), max(radial_indiv$deg), length.out = 100),
  quality_score = 1  
)

#hist(new_data$predicted)
quantile(new_data$predicted)

# Add a dummy observation level for prediction purposes
new_data$obs <- factor("001")  # Assuming a single observation for prediction
new_data$predicted <- predict(md1$gam, newdata = new_data, type = "response")
radial_indiv$residuals <- residuals(md1$gam, type = "deviance")

p0 <- ggplot(new_data, aes(x = dist, y = deg)) +
  geom_raster(aes(fill = predicted)) +
  scale_fill_gradient(low = "#3B9AB2", high = "#F21A00", breaks = seq(0.1, 0.9, by = 0.1)) +
  labs(x = "Distance from center patch (m)", y = "Downstream of center patch (°)", fill = "Fertilisation \n success") +
  geom_contour(aes(z = predicted), breaks = seq(0.1, 0.9, by = 0.1), color = "white") +
  theme_minimal() +
  scale_y_reverse() + 
  geom_hline(yintercept = rad_lines, color = "#E1AF00", lty = 2) +
  annotate("text", x = max(new_data$dist)+1, y = rad_lines[1]-5, label = "Radial line 2", 
           color = "#E1AF00", size = 5, hjust = 1.2) +  # Adjust position
  annotate("text", x = max(new_data$dist)+1, y = rad_lines[length(rad_lines)]-5, label = "Radial line 7", 
           color = "#E1AF00", size = 5, hjust = 1.2)+
  theme(
    legend.key.height = unit(1.5, "cm"), # Increase the height of legend keys
    axis.title = element_text(size = 13), # Increase axis label size
    axis.text = element_text(size = 12) # Increase tick label size
  )# add deviance
# p0 = p0 +  geom_point(data = radial_indiv, aes(x = dist, y = deg, color = residuals), size = 2)+
#   scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
p0 <- p0 +
  geom_point(data = radial_indiv, aes(x = dist, y = deg, color = prop), size = 2) +
  geom_text(data = radial_indiv, aes(x = dist, y = deg, label = round(prop, 2)), 
            color = "black", size = 3, vjust = -1)  # Adjust size and position
p0


#ggsave(p0, filename = '2023palau_fert.pdf',  path = "./plots", device = "pdf",  width = 10, height = 6)  #this often works better than pdf
