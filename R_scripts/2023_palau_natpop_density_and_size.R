# 2023 a hya pop density
# 3-5 depth on the bottom

# load libraries ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(truncnorm)
library(tidybayes)
library(spatstat)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code

# import data -------------------------------------------------------------
# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header,...)}
# data1=read.excel() #read clipboard from excel
# save(data1, file = file.path("./Rdata", "2023palau.natural.den.RData"))
load("./Rdata/2023palau.natural.den.RData")


# labeling and wrangling -----------------------------------------------
str(data1) # check data type is correct
# data1$raw.x <- as.numeric(as.character(data1$raw.x))
data1$side <- as.factor(as.character(data1$side)) # transect
data2 <- data1 %>% split(., data1$side)
data2$left$neg <- -data2$left$dist_m_x
data1 <- data.frame(y = data1$dist_m_y, x = c(data2$right$dist_m_x, data2$left$neg))

# Data exploration ------------------------------------------------------

p0 <- ggplot() +
  geom_point(data1, mapping = aes(x = x, y = y), position = position_jitter(width = .02, height = .02), alpha = 0.50, size = 3)
# p0 = p0 + facet_wrap(~trans)#+scale_x_log10(name ="XXXX")#+geom_smooth(data1, mapping = aes(x = raw.x, y = suc/tot))
# p0= p0+ scale_y_continuous( limits = c(0, 1))
p0

# spatial clustering ------------------------------------------------------
# mature

# data2 = split(data1, data1$mature)
# data3 = data2$m

rslt2 <- as.ppp(data1, W = owin(c(0, 50), c(-2.5, 2.5))) # apply ppp to list
# coord_list <- split(data3[,c("x", "y")], data3$trans)  #split by transect into list
# rslt2 <- lapply(coord_list, as.ppp, W = owin(c(0,10), c(-2,2))) #apply ppp to list
plot(as.solist(rslt2))
df <- data.frame(summary1 = c(summary(rslt2$`1`)[2], summary(rslt2$`2`)[2], summary(rslt2$`3`)[2])) # consolidates to a df
summary(rslt2)[2]
# Clake-evens (clusering)
donnelly <- unname(clarkevans(rslt2)[2]) # Clarke-Evans clustering index
donnelly # The clustering index measures the degree of clustering or dispersion of the point pattern,
# with values close to 1 indicating strong clustering and values close to 0 indicating a random or dispersed pattern.
clarkevans(rslt2)

plot(Kest(rslt2, correction = c("best"))) # Ripley’s K-function. At a given distance, how pattern relates to randomness (Kpois).
# If blackline is above, indicates clustering, if blackline below indiates ordering, if on red line, indicated random.

#try and find best cluster
fitted_models <- list(
  Poisson = ppm(rslt2, ~1),  # Null model (Complete Spatial Randomness)
  Thomas = kppm(rslt2 ~1, clusters = "Thomas", method = "clik2"),  # Thomas Process (clustered)
  Matern = kppm(rslt2 ~1, clusters = "MatClust", method = "clik2"),  # Matérn Process (clustered)
  HardCore = kppm(rslt2 ~1, clusters = "Cauchy", method = "clik2")  # Hard-Core Process (competition)
)

# Compare models using AIC (Lower AIC is better)
(model_comparison <- sapply(fitted_models, AIC))

fitted_thomas <- kppm(rslt2 ~1, clusters = "Thomas", method = "clik2") # Fit Thomas process model
fitted_matern <- kppm(rslt2 ~1, clusters = "MatClust", method = "clik2") # Fit MatClust process model

scale_thomas <- sqrt(fitted_thomas$par[2])  # Extract estimated scale parameter

#Thomas Clustering ones the best



# nearest neighbour
dist <- nndist(rslt2) # distance vector
quantile(dist) # quantile  #0.707
median(dist)
n_n.med <- lapply(dist, median) # 0.71m apart
n_n.mean_all <- mean(unlist(n_n.med))
n_n.sd_all <- sd(unlist(n_n.med))
print(paste0('distance: ', median(dist) ))
sd(dist)

# density of colonies
transect_width = 5
transect_length = 50
den1 = nrow(data1) / (transect_width * transect_length) # 0.304 col/m^2
print(paste0('density: ',den1))




# density of neigherest neighbout
plot(density(unlist(dist)))
par(mfrow = c(3, 1), mar = c(2, 2, 2, 2) + 0.1) # bottom, left, right, top
plot(envelope(rslt2, fun = Kest, nsim = 780, nrank = 20)) # Envelopes of K-function: This is a hypothesis test
plot(envelope(rslt2, Kest, correction = "Ripley", verbose = F))
plot(envelope(rslt2, Lest, correction = "Ripley", verbose = F)) # this is often preferred over K test (transformed K)
plot(density(rslt2), main = "Transect 1") # some of the clustering likely from the substrate

# plot -------------------------------------------------------------
# density
p1 <- ggplot() +
  geom_density(aes(dist), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = dist), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
# geom_vline(xintercept = 0, color = "red", lty = 2)+ theme_sleek2()
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.7))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1

# colony size --------------------------------------------------------------------

# 1 Import data -----------------------------------------------------------
# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header, na.strings=c("","-","na"),...)}
# data1 <- read.excel() #read clipboard from excel
# save(data1, file = file.path("./Rdata", "2023_palau_transect_all.RData"))
load("./Rdata/2023_palau_transect_all.RData")

# 2 Labelling and wrangling -----------------------------------------------
str(data1) # check data type is correct
data1$side <- as.factor(as.character(data1$side))
data1$id <- as.factor(as.character(data1$id))
data1$photo_ID <- as.factor(as.character(data1$photo_ID))

## Wrangling
# diameter
quantile(data1$mean, c(0.025, 0.5, 0.975), na.rm = T)
sd(data1$mean, na.rm = T)
# radius
quantile(data1$mean / 2, c(0.025, 0.5, 0.975), na.rm = T)
sd(data1$mean / 2, na.rm = T)

# density
p1 <- ggplot(data1) +
  geom_density(aes(mean), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = mean), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
# geom_vline(xintercept = 0, color = "red", lty = 2)+ theme_sleek2()
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.03))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1

# 2 Labelling and wrangling -----------------------------------------------
str(data1) # check data type is correct
# data1$raw.x <- as.numeric(as.character(data1$raw.x))
data1$side <- as.factor(as.character(data1$side)) # transect
data2 <- data1 %>% split(., data1$side)
data2$left$neg <- -data2$left$dist_m_x
df <- data.frame(y = data1$dist_m_y, x = c(data2$right$dist_m_x, data2$left$neg))
grid_size <- 0.01
df$buffer_radius <- rtruncnorm(nrow(df), a = 0.05, mean = 0.12975, sd = 0.06547)
# mean and sd from above
# Create a grid of points
grid_points <- expand.grid(x = seq(min(df$x), max(df$x), by = grid_size), y = seq(min(df$y), max(df$y), by = grid_size))

# Check for points within the buffer radius
points_in_buffer <- sapply(1:nrow(grid_points), function(i) {
  any(sqrt((df$x - grid_points$x[i])^2 + (df$y - grid_points$y[i])^2) < df$buffer_radius)
})

# Plotting
plot(0, 0, type = "n", xlim = c(min(df$x), max(df$x)), ylim = c(min(df$y), max(df$y)), xlab = "X", ylab = "Y")

# Plot grid points
points(grid_points$x, grid_points$y, pch = 15, col = ifelse(points_in_buffer, "red", "black"))

# Plot random points
points(df$x, df$y, pch = 19)

## Draw circles around grid points
# if(buffer_radius > 0) {
#   apply(df, 1, function(p) {
#     draw.circle(p[1], p[2], buffer_radius, border="blue")
#   })
# }

calculate_nearest_edge_distances <- function(df) {
  distances <- numeric(nrow(df))
  for (i in 1:nrow(df)) {
    other_points <- df[-i, ]
    dist_to_others <- sqrt((df$x[i] - other_points$x)^2 + (df$y[i] - other_points$y)^2)
    nearest_distance <- min(dist_to_others)
    edge_to_edge_distance <- nearest_distance - (df$buffer_radius + df$buffer_radius) # Adjust for varying radii
    distances[i] <- max(edge_to_edge_distance, 0)
  }
  return(distances)
}
(nn <- calculate_nearest_edge_distances(df) %>% quantile(., c(0.025, 0.5, 0.975)))

# CC
points_in_buffer <- sapply(1:nrow(grid_points), function(i) {
  any(sqrt((df$x - grid_points$x[i])^2 + (df$y - grid_points$y[i])^2) < df$buffer_radius) # Adjust condition for varying radii
})
total_points_in_buffer <- sum(points_in_buffer)
grid_area <- length(seq(0, 10, by = grid_size))^2
(percent_cover <- (total_points_in_buffer / grid_area) * 100)

# compare 2022 with 2023 crest slope --------------------------------------
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/allee_experiments/Rdata/2022_adult_nat_intercol.RData") # load 2022 adult
df3 <- data.frame(dist = dist, habitat = "slope")
df4 <- data.frame(dist = ad, habitat = "crest")
df5 <- rbind(df3, df4)
intercol_dist_natural = df5
save(intercol_dist_natural, file = file.path("./Rdata", "intercol_dist_natural.RData"))
load("./Rdata/intercol_dist_natural.RData") #intercol_dist_natural

# Compare  distributions
p1 <- ggplot(df5, aes(x = dist)) +
  geom_density(aes(group = habitat, color = habitat, fill = habitat), alpha = 0.3)
p1 <- p1 + scale_fill_manual(values = c("steelblue4", "orchid4", "red"), name = "Habitat", labels = c("Reef slope", "Reef crest"), guide = guide_legend(override.aes = list(color = NA))) +
  scale_color_manual(values = c("steelblue4", "orchid4", "steelblue1", "steelblue4", "grey", "grey", "grey", "grey"))
p1 <- p1 + geom_vline(xintercept = 0.707, linetype = "dashed", color = "steelblue4", size = 1)
p1 <- p1 + geom_vline(xintercept = 0.632, linetype = "dashed", color = "orchid4", size = 1)
p1 <- p1 + scale_y_continuous(name = "Density")
p1 <- p1 + scale_x_continuous(name = "Intercolonial distance (m)")
p1 <- p1 + theme_sleek2()
p1 <- p1 + theme(
  legend.position = c(0.80, 0.80),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.key.size = unit(1.5, "lines")
)
p1 <- p1 + guides(color = "none")
p1

# ggsave(p1, filename = 'adult_densities.pdf',  path = "./plots", device = 'pdf',  width = 8, height = 5)  #this often works better than pdf

