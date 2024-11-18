
## notes,I think the weighting should reflect the fertilisation % data. Not larvae number
# where is c2 and c4?
# not quite sure about the offset used in the glm yet. 


# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
library(ggmap)
library(units)
library(sf)
library(ggmap)
library(Hmisc)
library(purrr)
library(gridExtra)
library(magick)
library(geosphere)
library(plotrix)
library(networkD3)
library(circular)
library(CircStats)
library(dplyr)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1") # set theme in code


# 1 Import cervus out data -----------------------------------------------------------
(data1 <- read.csv(file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/4 Palau genetics mixing/Cervus", "parentage_out1.csv")))



# 2. Data prep -----------------------------------------------
str(data1) # check data type is correct
data1$offsp_id <- as.factor(as.character(data1$Offspring.ID))
data1$moth_id <- as.factor(as.character(data1$Mother.ID))
data1$fath_id <- as.factor(as.character(data1$Candidate.father.ID))
(data1 <- data1[grep("\\*", data1$Trio.confidence), ]) # Filter rows where 'Trio.confidence' contains '*'
data2 = data1 %>% dplyr::select(., c(offsp_id, moth_id, fath_id))
data2 = data2 %>% mutate(id = moth_id)  #mother


# 1 Import CERVUS out data -----------------------------------------------------------
load("./Rdata/data_gl_filtered.RData")  #data_gl_filtered.RData
meta = data_gl_filtered@other$ind.metrics
meta = meta  %>% rename(id2 = id) %>% mutate(id = paste0("X", id2)) 
meta2 <- meta %>% dplyr::select(c(id, lat, lon, genotype ))
meta2 <- meta2 %>% mutate(genotype = gsub("_5$", "_05", genotype))

# join tables 
join_df1 <- left_join(data2, meta2, by = 'id')  #add mother lat lon
join_df1 <- join_df1 %>% dplyr::select(-id) %>% mutate(id = fath_id)  #reset id to fathers
join_df2 <- left_join(join_df1, meta2, by = 'id')  #add father lat lon
nrow(join_df2)

#try normalizing to weight larvae
larvae_count <- table(join_df2$moth_id)
join_df2$normalised_weight <- 1 / larvae_count[join_df2$moth_id]

# join fert df
load("./Rdata/fert_data.RData") #data1
join_df2 = join_df2 %>% mutate(id = genotype.x)
join_df2  = left_join(join_df2, data1, by = 'id')


# calculate metrics -----------------------------------------------------
# Convert the data frame to two separate sf objects with WGS 84 CRS for  GPS systems
points_x <- st_as_sf(join_df2, coords = c("lon.x", "lat.x"), crs = 4326)
points_y <- st_as_sf(join_df2, coords = c("lon.y", "lat.y"), crs = 4326)
# Using UTM zone 53N for Palau (EPSG: 32653)
points_x_proj <- st_transform(points_x, crs = 32653)
points_y_proj <- st_transform(points_y, crs = 32653)
# calculate distances between parents
join_df2$dist <- st_distance(points_x_proj, points_y_proj, by_element = TRUE)
join_df2$dist_m = as.numeric(join_df2$dist )

# calculate angle between mum and fathers (not weighted)
bin_width <- 20 # Define the bin width (e.g., 20 degrees)
join_df2 <- join_df2 %>% mutate(angle = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) + 360) %% 360)  #angle of male realtive to female
ds_angle = 326
join_df2 <- join_df2 %>% mutate(ang_rel_ds = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) - ds_angle) %% 360, 
                                ang_rel_ds = ifelse(ang_rel_ds > 180, 360 - ang_rel_ds, ang_rel_ds))  # Convert to range (-180, 180)
join_df2 <- join_df2 %>% mutate(angle_binned = cut(angle, breaks = seq(0, 360, by = bin_width), include.lowest = TRUE, labels = seq(bin_width/2, 360 - bin_width/2, by = bin_width))) # Create a new column for the binned angles
angle_counts_binned <- as.data.frame(table(join_df2$angle_binned)) # Create a table of binned angles and their counts
colnames(angle_counts_binned) <- c("angle_binned", "count") # Rename columns for clarity
angle_counts_binned$angle_binned <- as.numeric(as.character(angle_counts_binned$angle_binned)) # Convert the binned angles to numeric for plotting
# to true north
polar.plot(lengths = angle_counts_binned$count, polar.pos = angle_counts_binned$angle_binned, radial.lim = c(0, max(angle_counts_binned$count)), 
           start = 90, lwd = 3, line.col = 4, clockwise = T) # Plot using polar coordinates with binned angles
# to water flow (currently to radial line 5 )
polar.plot(lengths = angle_counts_binned$count, polar.pos = angle_counts_binned$angle_binned, radial.lim = c(0, max(angle_counts_binned$count)), start = 236.5, lwd = 3, line.col = 4) # Plot using polar coordinates with binned angles

# perform circular stats
join_df2$angle_rad <- join_df2$angle * pi / 180
angles_circular <- circular(join_df2$angle_rad, units = "radians")
rayleigh.test(angles_circular)  #strong sig eefect


# analyses ----------------------------------------------------------------
#hist(join_df2$dist_m)  #unweighted
#plot(density(join_df2$dist_m))  #unweighted
(quan <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs=c(0, .25, .5, .83, 1)))
unname(quan[3])
# Units: [m]
# 0%       25%       50%       75%      100% 
# 0.000000  1.621345  3.124454 10.830502 30.690024 

quan_66 <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs = c(0.17, 0.83))
lower_66 <- unname(quan_66[1])
upper_66 <- unname(quan_66[2])
quan_95 <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs = c(0.025, 0.975))
lower_95 <- unname(quan_95[1])
upper_95 <- unname(quan_95[2])

## selfing
selfing_id = join_df2[which(join_df2$dist_m == 0),]  #self fertilisation
selfing_no = nrow(selfing_id)
(selfing_prop = sum(selfing_id$normalised_weight) / sum(join_df2$normalised_weight))  #normalised selfing
#most of these from individual (5/7) from 7_10. Interesting other fragments didnt have this though. 
sum(selfing_id$normalised_weight[c(1,3)]) /sum(join_df2$normalised_weight)  #normalised with Adult 7 removed

## non participants
meta3 <- meta2[complete.cases(meta2), ] # make sure import matches NA type
part1 = data.frame(id = c(join_df2$moth_id , join_df2$fath_id))  #find parental reps participating
part2 = left_join(part1, meta, by  = 'id') %>% dplyr::select(genotype ) %>%   distinct(genotype)
(nrow(part2))  #32 participants

#participants = data.frame(id = c(unique(join_df2$moth_id), unique(join_df2$fath_id)))
anti_join(meta3, part2,  by = 'genotype' ) %>%   distinct(genotype) %>% nrow()
#25 individuals not participating - but this needs to cross check against fert data which may show additional participating mothers (not sequenced)

# no. sires and dams
(distinct_sire = join_df2 %>% distinct(genotype.y) %>% nrow())   #25 sires
join_df2 %>% distinct(genotype.x) %>% nrow()   #19 sires

# median dams sired
join_df2 %>% group_by(genotype.y) %>%                # Group by genotype.y
  summarise(count = n_distinct(genotype.x)) %>% summarise(mean_count = mean(count))  
#mean dams fert per sire was 1.72
#
# how many distinct dams did top 10% sire
data3 = join_df2 %>% dplyr::select(genotype.y, genotype.x) %>% distinct() 

order = data3 %>%    # Select relevant columns                      # Remove duplicate rows
  count(genotype.y) %>%           # Add count of each genotype.y
  arrange(desc(n))  %>%  dplyr::select(-n)     
# Arrange by count in descending order

left_join(order, data3)
#top 3 sires -> 9 of 19 distinct dams = 47.4%. 



# stats -------------------------------------------------------------------

#poisson - to assess sequeced  larvae
unique_dams <- unique(data1$id, na.rm = T)   #all dams sampled for fert
unique_sires <- unique(na.omit(meta2$genotype))
all_pairs <- expand.grid(genotype.x = unique_dams, genotype.y = unique_sires)
#observed_crosses <- join_df2 %>% dplyr::select(genotype.x, genotype.y) %>% mutate(count = 1)
observed_crosses_aggregated <- join_df2 %>% group_by(genotype.x, genotype.y) %>% summarise(count = n(), .groups = 'drop') #agreegate multiple pairwise crosses
complete_data <-left_join(all_pairs, observed_crosses_aggregated, by = c("genotype.x", "genotype.y")) %>% mutate(count = replace_na(count, 0))
meta2_selected <- meta2 %>% mutate(genotype.x = genotype) %>% dplyr::select(genotype.x, lat, lon) %>% distinct() %>% na.omit()
meta2_selected_y <- meta2 %>% mutate(genotype.y = genotype) %>% dplyr::select(genotype.y, lat, lon) %>% distinct() %>% na.omit()
join_df3 = left_join(complete_data, meta2_selected, by = c("genotype.x"), relationship = "many-to-one") %>% rename(lat.x = lat, lon.x = lon)
join_df3 = left_join(join_df3, meta2_selected_y, by = c("genotype.y"), relationship = "many-to-one") %>% rename(lat.y = lat, lon.y = lon)
#remove c2 and c4 for the moment
join_df3 <- join_df3[complete.cases(join_df3), ] # make sure import matches NA type
points_x_proj <- st_as_sf(join_df3, coords = c("lon.x", "lat.x"), crs = 4326) %>% st_transform(crs = 32653)
points_y_proj <- st_as_sf(join_df3, coords = c("lon.y", "lat.y"), crs = 4326) %>% st_transform(crs = 32653)
join_df3$dist_m <- as.numeric(st_distance(points_x_proj, points_y_proj, by_element = TRUE))
#join_df3 <- join_df3 %>% mutate(angle = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) + 360) %% 360)
join_df3 <- join_df3 %>% mutate(ang_rel_ds = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) - ds_angle) %% 360, 
                                ang_rel_ds = ifelse(ang_rel_ds > 180, 360 - ang_rel_ds, ang_rel_ds))  # Convert to range (-180, 180), linear scale with 
data4 = data1 %>% dplyr::select(id, suc, prop) %>% na.omit() %>% rename(genotype.x = id)
#join_df4 = left_join(join_df3, data4, by = 'genotype.x') %>% mutate(p_s  = count/suc) %>% na.omit() 
#remove females with no fert
# Merge all dam-sire pairs with observed crosses without filtering
join_df4 <- join_df3 %>% left_join(data4, by = 'genotype.x') 
# Create structural_zero indicator
join_df4 <- join_df4 %>% mutate(structural_zero = ifelse(suc == 0, 1, 0))
# Calculate sampling fractions and adjust for sequencing proportion
join_df4 <- join_df4 %>%
  mutate(
    p_s = ifelse(suc > 0, count / suc, NA),            # Proportion of crosses
    p_sa = ifelse(!is.na(p_s), p_s * prop, NA),        # Adjusted by sequencing proportion
    p_s2 = ifelse(!is.na(p_sa) & p_sa > 0, p_sa, 1e-10) # Assign small p_sa where p_sa = 0
  )
# Assign weights: 1 for structural zeros, p_sa for sampling zeros
join_df4 <- join_df4 %>%
  mutate(
    weight = ifelse(structural_zero == 1, 1, p_sa)
  )

# Final data cleaning (optional: check for any remaining NAs)
join_df4 <- join_df4 %>%
  filter(!is.na(weight))

## data exploration
hist(join_df4$dist_m)
plot(join_df4$count~ (join_df4$dist_m))  
plot(join_df4$count~ (join_df4$ang_rel_ds))  
join_df4 = join_df4[which(join_df4$genotype.x != join_df4$genotype.y),]  #remove selfing

# #standard poisson
# poisson_model <- glm(count ~ dist_m + angle , family = poisson(link = "log"), data = join_df4)  #without offset
# summary(poisson_model)
# # poisson with offset
# poisson_model_offset <- glm(count ~ dist_m + angle + offset(log(p_s2)), family = poisson(link = "log"), data = join_df4)
# summary(poisson_model_offset)

# zero inflated
library(pscl)
zip_model <- zeroinfl(count ~ scale(dist_m) * scale(ang_rel_ds) | scale(dist_m) * scale(ang_rel_ds) , data = join_df4, dist = "poisson")
summary(zip_model)
zinb_model_negbin <- zeroinfl(count ~ scale(dist_m) * scale(ang_rel_ds) | scale(dist_m) * scale(ang_rel_ds), data = join_df4, dist = "negbin")
summary(zinb_model_negbin)
zip_model_offset_add <- zeroinfl(count ~ scale(dist_m) + scale(ang_rel_ds) + offset(log(p_s2))| scale(dist_m) + scale(ang_rel_ds) , data = join_df4, dist = "poisson")
summary(zip_model_offset_add)
zip_model_offset_int <- zeroinfl(count ~ scale(dist_m) * scale(ang_rel_ds) + offset(log(p_s2))| scale(dist_m) * scale(ang_rel_ds) , data = join_df4, dist = "poisson")
summary(zip_model_offset_int)
AIC(zip_model, zinb_model_negbin, zip_model_offset )

# Overdispersion
E2 <- resid(zip_model_offset_int, type = "pearson")
N  <- nrow(join_df4)
p  <- length(coef(zip_model_offset_int))  
sum(E2^2) / (N - p)


# plots -------------------------------------------------------------------

# weighted density
p1 <- ggplot(join_df2, aes(x = dist_m)) +
  #geom_density(aes(fill= 'steelblue4'), alpha = 0.3) +
  geom_density(aes(weight = normalised_weight, fill = 'steelblue4'), alpha = 0.3) + # Include weight adjustment
  #tidybayes::stat_pointinterval(aes(y = 0, x = dist_m), .width = c(.66, .95))+
  geom_errorbarh(aes(y = 0, xmin = lower_66, xmax = upper_66), height = 0.0, color = "black", linetype = "solid", size = 1.5) +
  geom_errorbarh(aes(y = 0, xmin = lower_95, xmax = upper_95), height = 0.0, color = "black", linetype = "solid", size = .5)+
  geom_point(aes(x = unname(quan[3]), y = 0), color = "black", size = 3, shape = 21, fill = "black")

p1 <- p1 + geom_vline(xintercept = unname(quan[3]), color = "red", lty = 2)
p1 <- p1 + theme_sleek2()
p1 <- p1 + scale_fill_manual(values = c("steelblue4", "white", "steelblue1", "white", "grey", "grey")) +
  scale_color_manual(values = c("steelblue4", "grey", "steelblue1", "steelblue4", "grey", "grey", "grey", "grey")) + theme(legend.position = "none") # nice
p1 <- p1 + scale_y_continuous(name = "Probability density")
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.08))
p1 <- p1 + scale_x_continuous(name = "Distance (m)")
p1 <- p1 + theme(strip.text.x = element_text(colour = "grey30", size = 8, vjust = -7))
p1 <- p1 + theme(panel.spacing.y = unit(-1.5, "lines"))
p1

#ggsave(p1, filename = '2023palau_intercol.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf



# create map --------------------------------------------------------------

# blank map, distance coloured
# create lines connecting dams and sires
lines_sf <- st_sfc(lapply(1:nrow(join_df2), function(i) {
  st_linestring(rbind(c(join_df2$lon.x[i], join_df2$lat.x[i]),
                      c(join_df2$lon.y[i], join_df2$lat.y[i])))}), crs = 4326)
lines_sf <- st_sf(geometry = lines_sf, dist_m = join_df2$dist_m)

ggplot() +
  geom_sf(data = points_x, color = 'blue', size = 2) +
  geom_sf(data = points_y, color = 'red', size = 2) +
  geom_sf(data = lines_sf, aes(color = dist_m), size = 1, alpha = 0.5) +
  scale_color_gradient(low="green", high="red") +
  labs(title="Pairwise Crosses Between Sires and Dams",
       color="Distance (m)") +
  theme_minimal()



# map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "color", maptype = "satellite")
# ggmap(map) +
#   #geom_segment(join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +
#   geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
#   labs(x = "Longitude", y = "Latitude") 
#ggsave(filename = '2023palau_design.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf


#zoomed out
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")

p2 = ggmap(map) +
  geom_segment(data = join_df2, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +  # pairwise parental cross
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +  # all colony positions
  geom_point(data = join_df2, aes(x = lon.x, y = lat.x, color = "Mother"), size = 3) +  # mother, mapped to "Mother"
  geom_point(data = join_df2, aes(x = lon.y, y = lat.y, color = "Father"), size = 3) +  # father, mapped to "Father"
  scale_color_manual(name = "Parent", values = c("Mother" = "blue", "Father" = "green")) +  # manual colour scale for the legend
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()#+
  #annotate("text", x = 134.4953, y = 7.3136, label = paste('test'), hjust = 1.1, vjust = -0.1, size = 5, color = "blue") # 
p2
#ggsave(filename = '2023palau_gen_zoomout.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf

#zoomed out jitter attempt (not sure if better)
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
p2 = ggmap(map) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +  # all colony positions
  geom_segment(data = join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +  # pairwise parental cross
  geom_jitter(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3, shape = 16, width = 0.000005, height = 0.000005) +  # mother, jittered slightly
  geom_jitter(data = join_df2, aes(x = lon.y, y = lat.y), alpha = 0.5, color = "green", size = 3, shape = 17, width = 0.000005, height = 0.000005) +  # father, different shape, jittered slightly
  labs(x = "Longitude", y = "Latitude") +
  theme_sleek1()


#zoomed in
map <- get_googlemap(center = c(lon = join_df2$lon.x[4], lat = join_df2$lat.x[4]), zoom = 22, color = "bw", maptype = "satellite")
p3 = ggmap(map) +
  geom_segment(join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +  #pairwise parental cross
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3)+   # all colony positions
  geom_point(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3)+  #mother
  geom_point(data = join_df2, aes(x = lon.y, y = lat.y), color = "blue", size = 3) +  #father
  labs(x = "Longitude", y = "Latitude") +
  theme_sleek1()
#ggsave(filename = '2023palau_gen_zoomin.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf

#colony positions
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_1 = ggmap(map) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  labs(x = "Longitude", y = "Latitude") +
  theme_sleek1())

#monitoring positions
# map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
# (p4_2 = ggmap(map) +
#     geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
#     geom_point(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3)+
#     labs(x = "Longitude", y = "Latitude") +
#     theme_sleek1())

#mother positions
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_2 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3)+
    labs(x = "Longitude", y = "Latitude") +
    theme_sleek1())

#father positions
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_3 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2, aes(x = lon.y, y = lat.y), color = "red", size = 3)+
    labs(x = "Longitude", y = "Latitude") +
    theme_sleek1())

#selfing positions
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_4 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2[which(join_df2$dist_m == 0),], aes(x = lon.y, y = lat.y), color = "green", size = 3)+
    labs(x = "Longitude", y = "Latitude") +
    theme_sleek1())


#grid.arrange(arrangeGrob(p2, p3, ncol = 1))
grid.arrange(arrangeGrob(p4_1, p4_2, ncol = 2), arrangeGrob(p4_3, p4_4, ncol = 2))


# parental pairwise crosses (based on mother genotypes) -------------------------------------------------------

# Define a function that creates a plot for a given genotype
plot_genotype <- function(genotype) {
  # Subset the data for the current genotype
  subset_df <- join_df2[join_df2$genotype.x == genotype, ]
  # Fetch the map based on the subset data (mother location)
  map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "color", maptype = "satellite")
  p = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +  # all colony positions
    geom_segment(data = subset_df, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +  # pairwise parental cross
    geom_point(data = subset_df, aes(x = lon.x, y = lat.x, color = "Mother"), size = 3) +  # mother, mapped to "Mother"
    geom_point(data = subset_df, aes(x = lon.y, y = lat.y, color = "Father"), size = 3) +  # father, mapped to "Father"
    scale_color_manual(name = "Parent", values = c("Mother" = "blue", "Father" = "green")) +  # manual colour scale for the legend
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal()+
  annotate("text", x = 134.49543, y = 7.3136, label = paste(genotype), hjust = 1.1, vjust = -0.1, size = 8, color = "blue") # annotation for mother genotype
  
  return(p)
}

plots <- map(unique(join_df2$genotype.x), plot_genotype)
# to print each plot:
#walk(plots, print)

#save plots
walk2(plots, unique(join_df2$genotype.x), ~ggsave(path = "./plots/pairwise_crosses", filename = paste0("genotype_", .y, ".jpg"), plot = .x, width = 8, 
                                                  height = 6))

## Add to animation

img_dir <- "./plots/pairwise_crosses"  # Adjust this path
image_files <- list.files(path = img_dir, pattern = "*.jpg", full.names = TRUE)
images <- image_read(image_files)
gif <- image_animate(images, fps = 0.5)  # 

# Save the animated GIF
#image_write(gif, path = "./plots/compiled_maps.gif")



# network sankey ----------------------------------------------------------

#find counts of pairs
links <- join_df2 %>%
  dplyr::select(genotype.x, genotype.y) %>%
  count(genotype.x, genotype.y) %>%
  mutate(source = paste0(genotype.x, "_dam"),  # Add 'x' to genotype.x
         target = paste0(genotype.y, "_sire")) %>%  # Add 'y' to genotype.y
  dplyr::select(source, target, n)  # Select relevant columns

#str(links)

# The unique node names
nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique())
str(nodes)

# match to numbers, not names
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
str(links)

#colour by sire radial
links$group = substr(links$target  , 1, 1)
links$group <- as.factor(links$group)

# Define a colour scale for the links based on the group factor
link_colour_scale <- 'd3.scaleOrdinal()
  .domain(["group1", "group2", "group3"])    # Replace with actual levels of links$group
  .range(["#FF5733", "#33FF57", "#3357FF"])' # Replace with your desired colours


p <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDtarget",
  Target = "IDsource",
  Value = "n",
  NodeID = "name",
  units = "TWh",
  fontSize = 12,
  nodeWidth = 10,
  LinkGroup = "group", 
  iterations = 1)        # ensure node order is as in data
p



# visnetwork --------------------------------------------------------------

library(visNetwork)
library(htmlwidgets)


nodes <- data.frame(id = unique(c(join_df2$genotype.x, join_df2$genotype.y)))
edges <- data.frame(from = join_df2$genotype.y, to = join_df2$genotype.x, value = join_df2$dist_m)

visNetwork = visNetwork(nodes, edges) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)


#saveWidget(visNetwork, "./plots/visNetwork.html", selfcontained = TRUE)



# igraph ------------------------------------------------------------------

library(igraph)
library(ggraph)

# Create edges dataframe
edges <- data.frame(
  from = join_df2$genotype.x,
  to = join_df2$genotype.y,
  weight = join_df2$dist_m
)

# Create graph
g <- graph_from_data_frame(d = edges, directed = TRUE)

# Plot
ggraph(g, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(aes(edge_width = weight), alpha=0.5, color="gray") +
  geom_node_point(size=5, aes(color=name)) +
  geom_node_text(aes(label = name), vjust=1.5, size=3) +
  theme_void() +
  guides(color="none")



# forecenet - not really showing much -------------------------------------


# Create Nodes data frame
Nodes <- links %>%
  dplyr::select(name = source, group) %>%
  bind_rows(links %>% dplyr::select(name = target, group)) %>%
  distinct(name, .keep_all = TRUE) %>%
  arrange(name)  # Optional: Arrange for better readability

# Ensure 'group' is treated as a factor for categorical coloring
Nodes$group <- as.factor(Nodes$group)

# Assign a unique ID to each node (0-based index)
Nodes <- Nodes %>%
  mutate(ID = 0:(nrow(Nodes) - 1))


# Create the force-directed network
forceNetwork(
  Links = links,
  Nodes = Nodes,
  Source = "IDsource",
  Target = "IDtarget",
  NodeID = "name",
  Group = "group",
  opacity = 0.8,
  zoom = TRUE,
  linkDistance = 100,
  charge = -300,
  fontSize = 12,
  linkWidth = 1.5,
  # Optionally, define a color scale
  colourScale = JS("d3.scaleOrdinal()
                    .domain([1,2,3,4,5,6,7,'c'])
                    .range(['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f'])")
)




