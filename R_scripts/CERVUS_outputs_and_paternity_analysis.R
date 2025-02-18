## Cervus outputs

##Proiorites
#add genetic relatedness between each parent pair inot the glm.
#filtering sorted but need to rerun (NOTE: '_05' SHOWNG IN SUMMARY_OUT.CSV LEADING TO LOST LARVAE. NEED TO FIX


## Notes
# Conflict between dDocent and standard. Need to resolve
#Check:Double check sum_count (all larvae) not greater than suc
#INCLUDING ZEROS
#Best approach so far is to include zeros, test for zeroinflation, then sipmply. So far, simple model works best
#I think imputation is working BUTR SEEMS DANGEROUS SO MUCH MISSING DATA, but may need to filter all pairwise crosses for non fragment combinations

#NOT INCLUDING ZEROS
    #no effect of some just because no data points  e.g angle - need to check


# where is c2 and c4?
# fro the network analysis, maybe try keeping all possible pairwise combination to provide structure, but colour the observ


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
library(reshape2)
library(spatstat)
library(glmmTMB)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek3") # set theme in code



# Import cervus out data and prep-----------------------------------------------------------
#(data1 <- read.csv(file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus", "parentage_out1.csv")))
#dDocent (make sure '.' between split words
(data1 <- read.csv(file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus", 
                                    "summary_out_docent.csv"), check.names = T))

data1$Mother.ID <- sub("_5", "_05", data1$Mother.ID)
data1$Candidate.father.ID <- sub("_5", "_05", data1$Candidate.father.ID)
str(data1) # check data type is correc
#View(data1)
data1$offsp_id <- as.factor(as.character(data1$Offspring.ID))
data1$moth_id <- as.factor(as.character(data1$Mother.ID))
data1$fath_id <- as.factor(as.character(data1$Candidate.father.ID))
(data1 <- data1[grep("\\*", data1$Trio.confidence), ]) # Filter rows where 'Trio.confidence' contains '*'
data2 = data1 %>% dplyr::select(., c(offsp_id, moth_id, fath_id))
data2 = data2 %>% mutate(id = moth_id)  #mother
data2 %>% group_by(id) %>% summarise(count = n()) #check number of offspring per mother
nrow(data2)  #108 larvae

# Import meta data -----------------------------------------------------------
load("./Rdata/data_gl.RData")  #data_gl.RData
meta = data_gl@other$ind.metrics
meta = meta  %>% rename(id2 = id) %>% mutate(id = paste0("X", id2)) 
meta2 <- meta %>% dplyr::select(c(id, lat, lon, genotype, total_mean_dia ))
#meta2 <- meta2 %>% mutate(genotype = gsub("_5$", "_05", genotype))
"X7_05_2" %in% meta2$id

# Compute pairwise distances
meta2_2 = meta2 %>% na.omit()
str(meta2_2)
unique_meta <- meta2_2 %>% distinct(genotype, .keep_all = TRUE) %>% mutate(genotype = as.character(genotype))# Keep one row per genotype
pairs_df <- expand.grid(genotype.x = as.character(unique_meta$genotype),
          genotype.y = as.character(unique_meta$genotype), stringsAsFactors=FALSE) %>% 
          filter(genotype.x < genotype.y) # Unique pairwise combinations
str(pairs_df)
pairs_df <- pairs_df %>% left_join(unique_meta %>% dplyr::select(genotype,lat,lon, total_mean_dia), by=c("genotype.x"="genotype")) %>% 
          rename(lat.x=lat,lon.x=lon) # Add lat/lon for genotype.x
pairs_df <- pairs_df %>% left_join(unique_meta %>% dplyr::select(genotype,lat,lon,total_mean_dia), by=c("genotype.y"="genotype")) %>% 
          rename(lat.y=lat,lon.y=lon) # Add lat/lon for genotype.y
points_x <- st_as_sf(pairs_df,coords=c("lon.x","lat.x"),crs=4326) # Create sf object for first point
points_y <- st_as_sf(pairs_df,coords=c("lon.y","lat.y"),crs=4326) # Create sf object for second point
points_x_proj <- st_transform(points_x,crs=32653) # Project first points to UTM Zone 53N
points_y_proj <- st_transform(points_y,crs=32653) # Project second points to UTM Zone 53N
pairs_df$dist <- st_distance(points_x_proj,points_y_proj,by_element=TRUE) # Calculate pairwise distances
pairs_df$dist_m <- as.numeric(pairs_df$dist) # Convert distance to numeric (meters)
ggplot(pairs_df, aes(x=dist_m)) + geom_density(fill="blue",alpha=0.3,colour="blue") + labs(x="Distance (m)",y="Probability Density",title="Pairwise Distance Distribution") + theme_minimal() # Plot density
quantile(pairs_df$dist_m)
meta2_2 <- meta2 %>% na.omit() %>% distinct(genotype, .keep_all=TRUE) 
meta2_2_sf <- st_as_sf(meta2_2, coords = c("lon", "lat"), crs = 4326) # Convert to sf
meta2_2_proj <- st_transform(meta2_2_sf, crs = 32653) # Project to UTM
coords <- st_coordinates(meta2_2_proj)
win <- owin(xrange = range(coords[,1]), yrange = range(coords[,2])) 
points_ppp <- ppp(coords[,1], coords[,2], window=win) 
density_map <- density(points_ppp, sigma = bw.diggle(points_ppp)) # Auto bandwidth
plot(density_map, main="Spatial Density of Population") 
points(points_ppp, pch=20, col="red") # Overlay points
mean_density <- mean(density_map$v) # Compute mean intensity (points per square meter)
#0.02 points per m^2 but not orientated
# join tables 
join_df1 <- left_join(data2, meta2, by = 'id')  #add mother lat lon
join_df1 <- join_df1 %>% dplyr::select(-id) %>% mutate(id = fath_id)  #reset id to fathers
join_df2 <- left_join(join_df1, meta2, by = 'id')  #add father lat lon
nrow(join_df2)

#add genetic relatedness (run 2a first)
adult_colonies1 = adult_colonies_sort
#adult_colonies1$Individual1 <- paste0("X", adult_colonies1$Individual1) # Add 'X' to Individual1
#adult_colonies1$Individual2 <- paste0("X", adult_colonies1$Individual2) # Add 'X' to Individual2
adult_colonies1 = adult_colonies1 %>% rename(genotype.x = Individual1, genotype.y = Individual2) %>% 
  mutate(genotype.x = as.factor(genotype.x), genotype.y = as.factor(genotype.y)) 
str(adult_colonies1)
adult_colonies1 <- adult_colonies1 %>%
  mutate(genotype_x_trimmed = sub("_[^_]*$", "", genotype.x)) # Remove final '_' and trailing characters from genotype.x
adult_colonies1 <- adult_colonies1 %>%
  mutate(genotype_y_trimmed = sub("_[^_]*$", "", genotype.y)) # Remove final '_' and trailing characters from genotype.y

genetic_dist_ave <- adult_colonies1 %>%
  group_by(genotype_x_trimmed, genotype_y_trimmed ) %>% summarise(mean_Distance = mean(Distance, na.rm = TRUE)) %>% # Compute mean Distance
  ungroup() %>% rename(genotype.x = genotype_x_trimmed , genotype.y = genotype_y_trimmed, gen_dist = mean_Distance ) %>% 
  data.frame()# Remove grouping
#add onto join_df3


#try normalizing to weight larvae. This means that multiple larvae from the same mother are given less weight.
# It corrects for more larvae being sequenced by chance.
larvae_count <- table(join_df2$moth_id)
join_df2$normalised_weight <- 1 / larvae_count[join_df2$moth_id]

# join fert df
load("./Rdata/2023_palau_fert.RData") #data1
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
(polar_plot = polar.plot(lengths = angle_counts_binned$count, polar.pos = angle_counts_binned$angle_binned, radial.lim = c(0, max(angle_counts_binned$count)), 
           start = 90, lwd = 3, line.col = 4, clockwise = T)) # Plot using polar coordinates with binned angles
#prepare for panelling
create_polar_plot <- function() {
  par(mar = c(2, 2, 2, 2))  # Set margins
  polar.plot(
    lengths = angle_counts_binned$count,
    polar.pos = angle_counts_binned$angle_binned,
    radial.lim = c(0, max(angle_counts_binned$count)),
    start = 90, 
    lwd = 3, 
    line.col = 4, 
    clockwise = TRUE
  )
}
# to water flow (currently to radial line 5 )
polar.plot(lengths = angle_counts_binned$count, polar.pos = angle_counts_binned$angle_binned, radial.lim = c(0, max(angle_counts_binned$count)), start = 236.5, lwd = 3, line.col = 4) # Plot using polar coordinates with binned angles

# perform circular stats
join_df2$angle_rad <- join_df2$angle * pi / 180
angles_circular <- circular(join_df2$angle_rad, units = "radians")
rayleigh.test(angles_circular)  #strong sig eefect


# Analyses ----------------------------------------------------------------

## Pairwise distances
#unweighted
(quan <- quantile(join_df2$dist_m, probs=c(0, .25, .5, .83, 1)))
unname(quan[3])
#hist(join_df2$dist_m)  #unweighted
#plot(density(join_df2$dist_m))  #unweighted
quan_66 <- wtd.quantile(join_df2$dist_m, probs = c(0.17, 0.83))
quan_95 <- wtd.quantile(join_df2$dist_m, probs = c(0.025, 0.975))
lower_66 <- quan_66[1]; upper_66 <- quan_66[2]
lower_95 <- quan_95[1]; upper_95 <- quan_95[2]
unweight_dist = join_df2$dist_m
#save(unweight_dist, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/2 Allee model and overview/Rdata", "palau_unweight_dist.RData"))

#weighted
(quan_w <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs=c(0, .25, .5, .83, 1)))
unname(quan[3])
# Units: [m]
# 0%       25%       50%       75%      100% 
# 0.000000  1.621345  3.124454 10.830502 30.690024 
quan_66_w <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs = c(0.17, 0.83))
quan_95_w <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs = c(0.025, 0.975))
lower_66_w <- quan_66_w[1]; upper_66_w <- quan_66_w[2]
lower_95_w <- quan_95_w[1]; upper_95_w <- quan_95_w[2]

# Overlay weighted and unweighted density plots with new styling
# Define the colour palette
cols <- c("Unweighted" = "#A2CFE3",  # pastel blue
          "Weighted" = "#E3A2A8")    # pastel red

# Overlay weighted and unweighted density plots with matching error bar colours
pairwise_dist_plot <- ggplot(join_df2, aes(x = dist_m)) +
  geom_density(aes(fill = "Unweighted", color = "Unweighted"), alpha = 0.3) + # Unweighted density
  geom_density(aes(weight = normalised_weight, fill = "Weighted", color = "Weighted"), alpha = 0.3) + # Weighted density
  geom_errorbarh(aes(y = 0, xmin = lower_66, xmax = upper_66, color = "Unweighted"), 
                 height = 0.0, linetype = "solid", size = 1.5, alpha = 0.5) + # Unweighted error bar
  geom_errorbarh(aes(y = 0, xmin = lower_95, xmax = upper_95, color = "Unweighted"), 
                 height = 0.0, linetype = "solid", size = .5, alpha = 0.5) + # Unweighted error bar
  geom_errorbarh(aes(y = 0, xmin = lower_66_w, xmax = upper_66_w, color = "Weighted"), 
                 height = 0.0, linetype = "solid", size = 1.5, alpha = 0.5) + # Weighted error bar
  geom_errorbarh(aes(y = 0, xmin = lower_95_w, xmax = upper_95_w, color = "Weighted"), 
                 height = 0.0, linetype = "solid", size = .5, alpha = 0.5) + # Weighted error bar
  annotate("point", x = unname(quan[3]), y = 0, color = cols["Unweighted"], size = 3, shape = 21, fill = cols["Unweighted"]) + # Unweighted median
  annotate("point", x = unname(quan_w[3]), y = 0, color = cols["Weighted"], size = 3, shape = 21, fill = cols["Weighted"]) + # Weighted median
  geom_vline(xintercept = unname(quan[3]), color = '#A2CFE3', lty = 2) + # Unweighted median line
  geom_vline(xintercept = unname(quan_w[3]), color = "#E3A2A8", lty = 2) + # Weighted median line
  scale_fill_manual(values = cols, name = "Density Type") +
  scale_color_manual(values = cols, name = "Density Type") +
  coord_cartesian(ylim = c(0.0, 0.08)) +
  scale_x_continuous(name = "Distance (m)") +
  scale_y_continuous(name = "Probability density") +
  theme_sleek2() +
  theme(
    legend.position = c(0.8, 0.9),
    legend.text = element_text(size = rel(1), colour = "grey20"),
    strip.text.x = element_text(colour = "grey30", size = 8, vjust = -7),
    panel.spacing.y = unit(-1.5, "lines")
  )

pairwise_dist_plot

#ggsave(p1, filename = '2023palau_intercol.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf

##weighting for possible pairwise combinations (overrepresetns under values)
dens_est <- density(pairs_df$dist_m,from=min(pairs_df$dist_m),to=max(pairs_df$dist_m),n=512) # Estimate density of all available pairwise distances
dens_func <- approxfun(dens_est$x,dens_est$y,rule=2) # Create an interpolation function for density values
join_df2$weight_distance <- 1/dens_func(join_df2$dist_m) # Assign a weight to each observed distance (inverse of estimated density)
weighted_quantiles <- wtd.quantile(join_df2$dist_m, weights = join_df2$weight_distance, probs=c(0,0.25,0.5,0.83,1)) # Compute weighted quantiles adjusting for sampling bias
weighted_quantiles # Display the weighted quantiles
join_df2_2 = join_df2 %>%  dplyr::select(c(dist_m, weight_distance))
ggplot(join_df2, aes(x=dist_m)) + geom_density(aes(fill="blue", alpha=0.3, colour="blue", weight = weight_distance))+
  labs(x="Distance (m)",y="Probability Density",title="Pairwise Distance Distribution") + theme_minimal() # Plot density

join_df2$dist_m



## selfing (see below)
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

# mean dams sired
join_df2 %>% group_by(genotype.y) %>%                # Group by genotype.y
  summarise(count = n_distinct(genotype.x)) %>% summarise(mean_count = mean(count))  
#mean dams fert per sire was 1.72

# mean males per dam (small breeding units)
join_df2 %>% group_by(genotype.x) %>%                # Group by genotype.y
  summarise(count = n_distinct(genotype.y)) %>% summarise(mean_count = mean(count))  
#mean dams fert per sire was 2.26
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

## 3 types of zero. 1) structural (impossible), 2) missing data (only subsampled), 3) true zeros if sampled whole container
# model successful pos cross vs all unsuccessful crosses

## prepare data fro analysis of every pairwise combination as counts. 
unique_dams <- unique(data1$id, na.rm = T)   #all dams sampled for fert
unique_sires <- unique(na.omit(meta2$genotype))
#create all pairwise combinations
all_pairs <- expand.grid(genotype.x = unique_dams, genotype.y = unique_sires)
#observed_crosses <- join_df2 %>% dplyr::select(genotype.x, genotype.y) %>% mutate(count = 1)
observed_crosses_aggregated <- join_df2 %>% group_by(genotype.x, genotype.y) %>% summarise(count = n(), .groups = 'drop') #agreegate multiple pairwise crosses
complete_data <-left_join(all_pairs, observed_crosses_aggregated, by = c("genotype.x", "genotype.y")) %>% mutate(count = replace_na(count, 0))
meta2_selected <- meta2 %>% mutate(genotype.x = genotype) %>% dplyr::select(genotype.x, lat, lon, total_mean_dia) %>% distinct() %>% na.omit()
meta2_selected_y <- meta2 %>% mutate(genotype.y = genotype) %>% dplyr::select(genotype.y, lat, lon, total_mean_dia) %>% distinct() %>% na.omit()
join_df3 = left_join(complete_data, meta2_selected, by = c("genotype.x"), relationship = "many-to-one") %>% rename(lat.x = lat, lon.x = lon)
join_df3 = left_join(join_df3, meta2_selected_y, by = c("genotype.y"), relationship = "many-to-one") %>% rename(lat.y = lat, lon.y = lon)
#remove c2 and c4 for the moment
join_df3 <- join_df3[complete.cases(join_df3), ] # make sure import matches NA type

#add on genetic distances
join_df3 = left_join(join_df3, genetic_dist_ave, by = c("genotype.x", "genotype.y"))
head(join_df3)
join_df3 %>% filter(count > 0 & gen_dist <200)



points_x_proj <- st_as_sf(join_df3, coords = c("lon.x", "lat.x"), crs = 4326) %>% st_transform(crs = 32653)
points_y_proj <- st_as_sf(join_df3, coords = c("lon.y", "lat.y"), crs = 4326) %>% st_transform(crs = 32653)
join_df3$dist_m <- as.numeric(st_distance(points_x_proj, points_y_proj, by_element = TRUE))
join_df3 <- join_df3 %>% mutate(ang_rel_ds = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) - ds_angle) %% 360, 
                                ang_rel_ds = ifelse(ang_rel_ds > 180, 360 - ang_rel_ds, ang_rel_ds),
                                angle = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) + 360) %% 360)  # Convert to range (-180, 180), linear scale with 
join_df3 <- join_df3 %>% mutate(cos_ang = abs(cos(ang_rel_ds * pi / 180)))  # Converts angles into a continuous variable, 0 = perp, 1 = upstream/downstream

data4 = data1 %>% dplyr::select(id, suc, prop) %>% na.omit() %>% rename(genotype.x = id)
#join_df4 = left_join(join_df3, data4, by = 'genotype.x') %>% mutate(p_s  = count/suc) %>% na.omit() 
#remove females with no fert
# Merge all dam-sire pairs with observed crosses without filtering
join_df4 <- join_df3 %>% left_join(data4, by = 'genotype.x') %>% na.omit() #here suc is only use as guidance
join_df4  %>% group_by(genotype.x) %>% summarise(count = sum(count), suc = first(suc)) %>% mutate(diff = suc - count) %>% data.frame()


## data exploration
hist(join_df4$dist_m)
plot(join_df4$count~ (join_df4$dist_m))  
plot(join_df4$count~ (join_df4$ang_rel_ds))  
plot(join_df4$count ~ (join_df4$cos_ang ))
join_df4$count
join_df4 = join_df4[which(as.character(join_df4$genotype.x) != as.character(join_df4$genotype.y)),]  # remove selfing
str(join_df4)
# #standard poisson
# poisson_model <- glm(count ~ dist_m + angle , family = poisson(link = "log"), data = join_df4)  #without offset
# summary(poisson_model)
# # poisson with offset
# poisson_model_offset <- glm(count ~ dist_m + angle + offset(log(p_s2)), family = poisson(link = "log"), data = join_df4)
# summary(poisson_model_offset)

## truncated count (modelling only the magnitude of success once success is possible/observed)
df_pos_only <- subset(join_df4, count > 0)
# try binning to flow
df_pos_only <- df_pos_only %>% mutate(align_cat = case_when(
    (angle >= 281 | ang_rel_ds <= 11) ~ "Aligned",  # Within 45° of 326° (flow direction)
    (angle >= 101 & ang_rel_ds <= 191) ~ "Aligned", # Within 45° of 146° (downstream)
    TRUE ~ "Not Aligned"  # Everything else
  ))
df_pos_only$align_cat <- factor(df_pos_only$align_cat, levels = c("Not Aligned", "Aligned"))  # Reference: "Not Aligned"
table(df_pos_only$align_cat)

plot(df_pos_only$count  ~ df_pos_only$dist_m)
plot(df_pos_only$count  ~ df_pos_only$ang_rel_ds)
plot(df_pos_only$count  ~ df_pos_only$cos_ang)
summary(glm(df_pos_only$count  ~ df_pos_only$cos_ang))
range(df_pos_only$ang_rel_ds)
library(countreg)
trunc_pois <- zerotrunc(count ~ scale(dist_m) + align_cat, data = df_pos_only,  dist = "poisson") 
summary(trunc_pois)
# truncated negative binomial
trunc_nb <- zerotrunc(count ~ scale(dist_m) + cos_ang , data = df_pos_only, dist = "negbin") # truncated negative binomial
summary(trunc_nb) 
AIC(trunc_pois, trunc_nb)

# trunc_pois_wei <- zerotrunc(count ~ scale(dist_m) + scale(ang_rel_ds), data = df_pos_only, 
#                         dist = "poisson", weights = suc)
# summary(trunc_pois_wei)
# trunc_nb_wei  <- zerotrunc(count ~ scale(dist_m) + scale(ang_rel_ds), data = df_pos_only, 
#                       dist = "negbin", weights = suc)
# summary(trunc_nb_wei)
trunc_pois_offset <- zerotrunc(count ~ scale(dist_m) + cos_ang + offset(log(suc)), 
                               data = df_pos_only, dist = "poisson")
summary(trunc_pois_offset)
trunc_nb_offset <- zerotrunc(count ~ scale(dist_m) + cos_ang + offset(log(suc)), 
                             data = df_pos_only, dist = "negbin")
summary(trunc_nb_offset)
AIC(trunc_pois_wei, trunc_nb_wei,trunc_pois_offset, trunc_nb_offset )


## try on colonies sire from central patch
df_pos_only_centre <- df_pos_only[grep('c',df_pos_only$genotype.y),]
trunc_pois <- zerotrunc(count ~ scale(dist_m) + cos_ang, data = df_pos_only_centre,  dist = "poisson") 
summary(trunc_pois)
# truncated negative binomial
trunc_nb <- zerotrunc(count ~ scale(dist_m) + cos_ang , data = df_pos_only_centre, dist = "negbin") # truncated negative binomial
summary(trunc_nb) 
AIC(trunc_pois, trunc_nb)
plot(df_pos_only_centre$count  ~ df_pos_only_centre$dist_m)
plot(df_pos_only_centre$count  ~ df_pos_only_centre$cos_ang)

# zero inflated
# library(pscl)
# # Zero-Inflated Poisson (ZIP) model
# zip_model <- zeroinfl(count ~ scale(dist_m) * scale(ang_rel_ds) | scale(dist_m) * scale(ang_rel_ds) , data = join_df4, dist = "poisson")
# summary(zip_model)
# # Zero-Inflated Negative Binomial (ZINB) model
# zinb_model_negbin <- zeroinfl(count ~ scale(dist_m) * scale(ang_rel_ds) | scale(dist_m) * scale(ang_rel_ds), data = join_df4, dist = "negbin")
# summary(zinb_model_negbin)
# # ZIP model with an offset added to the count component
# zip_model_offset_add <- zeroinfl(count ~ scale(dist_m) + scale(ang_rel_ds) + offset(log(p_s2))| scale(dist_m) + scale(ang_rel_ds) , data = join_df4, dist = "poisson")
# summary(zip_model_offset_add)
# # ZIP model with an offset and interaction in predictors
# zip_model_offset_int <- zeroinfl(count ~ scale(dist_m) * scale(ang_rel_ds) + offset(log(p_s2))| scale(dist_m) * scale(ang_rel_ds) , data = join_df4, dist = "poisson")
# summary(zip_model_offset_int)
# AIC(zip_model, zinb_model_negbin, zip_model_offset)

# Overdispersion
# E2 <- resid(zip_model_offset_int, type = "pearson")
# N  <- nrow(join_df4)
# p  <- length(coef(zip_model_offset_int))  
# sum(E2^2) / (N - p)


## include zeros and imputations
#create three types of zeros. Note ive allowed some flexibility in strcutural since probabiliy is low
join_df5 = join_df4 %>% group_by(genotype.x) %>% mutate(
  # Sum of observed positive counts for this dam
  sum_count = sum(count, na.rm = TRUE),
  # Logic to classify each row:
  final_count = case_when(
    # 1) If suc == 0 OR suc < 5 => structural zero
    suc == 0 | suc < 5 ~ 0,
    # 2) If sum(all positive count) == suc and current count=0 => true zero
    (sum_count == suc & count == 0) ~ 0,
    # 3) Otherwise if count=0 => set NA (missing)
    (sum_count < suc & count == 0) ~ 0,
    # If it's already > 0, keep as is
    TRUE ~ as.numeric(count)
  ),
  zero_type = case_when(
    suc == 0 | suc < 5 ~ "structural_zero", 
    (sum_count == suc & count == 0) ~ "true_zero", 
    (sum_count < suc & count == 0) ~ "missing_zero", 
    count > 0 ~ "positive_count", 
    TRUE ~ NA_character_
  )
) %>% ungroup() %>% data.frame()


table(join_df5$zero_type)
str(join_df5)

#check no mismatches
join_df5 %>% filter(is.na(zero_type)) %>% group_by(genotype.x) %>% summarise(count_na = n()) # Count rows per group
# Keep only rows that are true_zero or positive_count

plot(final_count ~ dist_m, data = join_df5)
plot(final_count ~ cos_ang, data = join_df5)


# zero-inflated models ----------------------------------------------------
join_df5$weight <- join_df5$sum_count / join_df5$suc  # Define the sampling fraction as a weight
join_df5$weight <- ifelse(join_df5$sum_count == join_df5$suc, 
                          1,  # Assign weight of 1 when sum_count == suc
                          join_df5$sum_count / join_df5$suc)  # Otherwise, use the existing formula

plot(join_df5$final_count ~ (join_df5$gen_dist))
plot(join_df5$final_count ~ (join_df5$total_mean_dia.y))
plot(join_df5$final_count ~ (join_df5$cos_ang))

#ZINB
global_model <- glmmTMB(final_count ~ dist_m * cos_ang + total_mean_dia.y + log(gen_dist) + offset(log(suc + 1e-6)) ,  
                        ziformula = ~ suc,  # Include zero-inflation
                        family = nbinom2,
                        data = join_df5,
                        #weights = weight,
                        na.action = na.fail)  # Ensure NA handling works with MuMIn

summary(global_model)

# Generate new data with dist_m varying and all other predictors held at their mean
new_data <- data.frame(
  dist_m = seq(min(join_df5$dist_m), max(join_df5$dist_m), length.out = 100),
  cos_ang = mean(join_df5$cos_ang, na.rm = TRUE),
  total_mean_dia.y = mean(join_df5$total_mean_dia.y, na.rm = TRUE),
  gen_dist = mean(join_df5$gen_dist, na.rm = TRUE),
  suc = mean(join_df5$suc, na.rm = TRUE),  # Offset term
  weight = mean(join_df5$weight, na.rm = TRUE)  # Ensure weight is included
)
preds <- predict(global_model, newdata = new_data, type = "response", se.fit = TRUE)
# Store predictions and confidence intervals
new_data$preds <- preds$fit
new_data$lower <- preds$fit - 1.96 * preds$se.fit  # 95% CI lower bound
new_data$upper <- preds$fit + 1.96 * preds$se.fit  # 95% CI upper bound

ggplot() +
  geom_point(join_df5, mapping = aes(y = final_count, x = dist_m))+
  geom_line(new_data, mapping = aes(x = dist_m, y = preds )) +
  geom_ribbon(new_data, mapping = aes(x = dist_m, ymin = lower, ymax = upper),fill = 'blue', alpha = 0.1) +  # Confidence intervals
  labs(x = "Distance (m)", y = "Final Count") +
  theme_minimal()

hist(join_df5$final_count)




AIC(global_model)
dredged_models <- dredge(global_model, rank = "AIC")
(best_models <- subset(dredged_models, delta < 2))  # Select models within ΔAIC < 2
#best model is dist_m * cos_ang 

global_model_red <- glmmTMB(final_count ~ dist_m * cos_ang,  
                            ziformula = ~ suc,  # Include zero-inflation
                            family = nbinom2,
                            data = join_df5,
                            #weights = weight,
                            na.action = na.fail)  # Ensure NA handling works with MuMIn

summary(global_model_red)
AIC(global_model_red)

global_model_noinf <- glmmTMB(final_count ~ dist_m * cos_ang + total_mean_dia.y + gen_dist ,  
                        family = nbinom2,
                        data = join_df5,
                        weights = weight,
                        na.action = na.fail)  # Ensure NA handling works with MuMIn
summary(global_model)
dredged_models <- dredge(global_model, rank = "AIC")
(best_models <- subset(dredged_models, delta < 2))  # Select models within ΔAIC < 2
#mod_zinb_weighted_nodia best



sum(resid(mod_zinb, type = "pearson")^2) / (nrow(data1) - length(coef(mod_zinb))) # (overdispsrion glm Zuur)
library(RVAideMemoire) # GLMM overdispersion test
library(performance)
# check_overdispersion(md1)  #only for count data
plot(fitted(mod_zinb), resid(mod_zinb)) # fitted vs residuals
abline(h = 0)
performance::r2(mod_zinb)
icc(mod_zinb) # Intraclass Correlation Coefficient
check_zeroinflation(mod_zinb_weighted_nodia)
check_singularity(mod_zinb)
check_model(mod_nb)
library(DHARMa)
sim_res <- simulateResiduals(mod_zinb_weighted_nodia)
plot(sim_res)
check_outliers(mod_zinb_weighted_nodia)

# #remove missing and structural
# join_df4 <- join_df4 %>%
#   filter(!(zero_type %in% c("structural_zero", "missing_zero"))) # Keep only rows that are true_zero or positive_count
# 

# imputation --------------------------------------------------------------
## attempting imputation
# Ensure `count` is numeric
join_df5$count <- as.numeric(join_df5$count)
# Only replace `count` with NA if it's a missing zero
join_df5$count[join_df5$zero_type == "missing_zero"] <- NA
# Check that only missing_zero values are NA
table(join_df5$zero_type, is.na(join_df5$count))
library(mice)
meth <- make.method(join_df5)

meth["count"] <- "rf"  # select method for imputatins
imputed_data <- mice(join_df5, method = meth, m = 5, seed = 123)
summary(imputed_data)
# Extract first imputed dataset
imputed_df <- complete(imputed_data, 1) # First imputed dataset
# Check if `structural_zero` or `true_zero` values were changed (they should NOT be!)
table(imputed_df$zero_type, is.na(imputed_df$count))
# Compare before and after imputation
table(join_df5$zero_type)
table(imputed_df$count)
table(join_df5$zero_type, imputed_df$count)
library(ggplot2)
ggplot(imputed_df, aes(x = count, fill = zero_type)) + geom_bar(position = "dodge") +
  theme_minimal() + labs(title = "Distribution of Imputed vs. Observed Counts",
       x = "Count (Number of Eggs)",
       y = "Frequency",
       fill = "Zero Type")

# Fit negative binomial model on each imputed dataset
nb_models <- with(imputed_data, glm.nb(count ~ dist_m + ang_rel_ds, data = join_df5))
summary(nb_models)
pooled_nb_model <- pool(nb_models)
summary(pooled_nb_model)

#exclude structural
model_data <- subset(imputed_df, zero_type != "structural_zero")
nb_models <- with(model_data, glm.nb(count ~ dist_m + ang_rel_ds, data = model_data))
summary(nb_models)



# # Create structural_zero indicator
# join_df4 <- join_df4 %>% mutate(structural_zero = ifelse(suc == 0, 1, 0)) #think this isnt the best way to assign strucural zeros
# # Calculate sampling fractions and adjust for sequencing proportion
# join_df4 <- join_df4 %>%
#   mutate(
#     p_s = ifelse(suc > 0, count / suc, NA),            # Proportion of crosses
#     p_sa = ifelse(!is.na(p_s), p_s * prop, NA),        # Adjusted by sequencing proportion
#     p_s2 = ifelse(!is.na(p_sa) & p_sa > 0, p_sa, 1e-10) # Assign small p_sa where p_sa = 0
#   )
# # Assign weights: 1 for structural zeros, p_sa for sampling zeros
# join_df4 <- join_df4 %>%
#   mutate(
#     weight = ifelse(structural_zero == 1, 1, p_sa)
#   )
# # Final data cleaning (optional: check for any remaining NAs)
# join_df4 <- join_df4 %>%
#   filter(!is.na(weight))

# ZIB mixture model
###prepare binomial outcome. This crunches multiple observations to 1. 
join_df6 = join_df5 %>% mutate(cross_occurrence = if_else(is.na(final_count), 0L, if_else(final_count > 0, 1L, 0L))) %>% 
  dplyr::select(genotype.x, genotype.y, cross_occurrence, suc, sum_count, zero_type, count, lat.x, lon.x, lat.y, lon.y, dist_m, ang_rel_ds, cos_ang, ang_no_flow, angle)  # Create binary outcome and select key columns
join_df6 <- join_df6 %>% mutate(trials = 1)  # Add a column indicating one trial per observation
library(glmmTMB) # Load glmmTMB package
mod_zib <- glmmTMB(cross_occurrence ~ dist_m + cos_ang, # Conditional model predictors
                   ziformula = ~ suc, # Zero‐inflation component predictor (e.g. fertilisation success. 
                   family = binomial, # Specify binomial family for binary outcome
                   data = join_df6) # Use prepared dataset with binary outcome
summary(mod_zib) # Display model summary


library(brms)  # Load the brms package
# Define the model formula with zero-inflation
zib_formula <- bf(
  cross_occurrence | trials(1) ~ dist_m + cos_ang,  # Conditional model (logit-link Bernoulli)
  zi ~ suc  # Zero-inflation modeled as a function of suc
)
# Set priors for both the conditional and zero-inflation components
priors <- c(
  set_prior("normal(-4, 1)", class = "Intercept"),               # Main model intercept
  set_prior("normal(0, 1)", class = "b"),                        # Priors for dist_m and cos_ang effects
  set_prior("normal(4, 1)", class = "Intercept", dpar = "zi"),   # Zero-inflation intercept (high probability of structural zeros)
  set_prior("normal(-2, 0.5)", class = "b", dpar = "zi")         # Effect of suc on zero-inflation
)
# Fit the Zero-Inflated Binomial (ZIB) model
m <- brm(
  formula = zib_formula,  # Use the defined formula
  data = join_df6,  # Use dataset with binary response (0/1)
  family = zero_inflated_binomial(),  # Zero-inflated binomial family
  prior = priors,  # Apply informed priors
  chains = 4, cores = 4, iter = 2000, control = list(adapt_delta = 0.95)  # Improve sampling stability
)
summary(m)






# plots -------------------------------------------------------------------




# create map --------------------------------------------------------------

# blank map, distance coloured
# create lines connecting dams and sires
lines_sf <- st_sfc(lapply(1:nrow(join_df2), function(i) {
  st_linestring(rbind(c(join_df2$lon.x[i], join_df2$lat.x[i]),
                      c(join_df2$lon.y[i], join_df2$lat.y[i])))}), crs = 4326)
lines_sf <- st_sf(geometry = lines_sf, dist_m = join_df2$dist_m)

distance_plot = ggplot() +
  geom_sf(data = points_x, color = 'blue', size = 2) +
  geom_sf(data = points_y, color = 'red', size = 2) +
  geom_sf(data = lines_sf, aes(color = dist_m), size = 1, alpha = 0.5) +
  scale_color_gradient(low="green", high="red") +
  labs(title="Pairwise Crosses Between Sires and Dams",
       color="Distance (m)") +
  theme_minimal()

## add to google map
# Get the Google map
lon_range <- range(c(join_df2$lon.x, join_df2$lon.y))
lat_range <- range(c(join_df2$lat.x, join_df2$lat.y))

# Add small buffer (adjust the 0.0001 value to fit your needs)
buffer <- 0.0001
map <- get_googlemap(
  center = c(
    lon = mean(lon_range),
    lat = mean(lat_range)
  ),
  zoom = 20,
  color = "color",
  maptype = "satellite",
  # Define bounds explicitly
  bounds = c(
    left = min(lon_range) - buffer,
    bottom = min(lat_range) - buffer,
    right = max(lon_range) + buffer,
    top = max(lat_range) + buffer
  )
)
                 
# Create the base ggmap
base_map <- ggmap(map)
# Overlay the distance plot on the map
# Add jitter to the points data
# First convert sf objects to data frames for jittering if needed
points_y_jittered <- st_jitter(points_y, amount = 0.000002) # Adjust amount as needed
points_x_jittered <- st_jitter(points_x, amount = 0.000002) # Adjust amount as needed

# Overlay the distance plot on the map with jittered points and adjusted alpha
distance_map_plot <- base_map +
  geom_sf(data = lines_sf, 
          aes(color = dist_m), 
          #size = 30,  
          linewidth = 1,# Thicker lines
          inherit.aes = FALSE) + 
  geom_sf(data = points_y_jittered, 
          color = 'red', 
          size = 2.5,           
          alpha = 0.7,          
          inherit.aes = FALSE) + 
  geom_sf(data = points_x, 
          color = 'blue', 
          size = 2.5,           
          alpha = 0.7,          
          inherit.aes = FALSE) + 
  scale_color_gradient(
    low = "#E6E6FA",    # Light purple
    high = "#4B0082"    # Deep/hot purple
  ) +
  labs(title = "Pairwise Crosses Between Sires and Dams", 
       color = "Distance (m)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()

distance_map_plot


map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "color", maptype = "satellite")
ggmap(map) +
  #geom_segment(join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  labs(x = "Longitude", y = "Latitude")
#ggsave(filename = '2023palau_design.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf


#zoomed out
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")

p2 = ggmap(map) +
  geom_segment(data = join_df2, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "grey", size = 1) +  # pairwise parental cross
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +  # all colony positions
  geom_point(data = join_df2, aes(x = lon.x, y = lat.x, color = "Mother"), size = 3) +  # mother, mapped to "Mother"
  geom_point(data = join_df2, aes(x = lon.y, y = lat.y, color = "Father"), size = 3) +  # father, mapped to "Father"
  scale_color_manual(name = "Parent", values = c("Mother" = "blue", "Father" = "red")) +  # manual colour scale for the legend
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
  labs(title = "All positions", x = "Longitude", y = "Latitude") +
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
    labs(title = "Mother", x = "Longitude", y = "Latitude") +
    theme_sleek1())

#father positions
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
join_df2_sire <- join_df2 %>% distinct(lon.y, lat.y, .keep_all = TRUE)
(p4_3 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2_sire, aes(x = lon.y, y = lat.y), color = "red", size = 3, alpha = 0.5)+
    labs(title = "Father", x = "Longitude", y = "Latitude") +
    theme_sleek1())

#selfing positions
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_4 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2[which(join_df2$dist_m == 0),], aes(x = lon.y, y = lat.y), color = "green", size = 3)+
    labs(title = "Selfing", x = "Longitude", y = "Latitude") +
    #annotate("text", x = join_df2$lon.x[8], y = join_df2$lat.x[8] + 0.0002, label = "D", size = 6, colour = "white") +
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
    geom_segment(data = subset_df, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "grey", size = 1) +  # pairwise parental cross
    geom_point(data = subset_df, aes(x = lon.y, y = lat.y, color = "Father"), size = 3) +  # father, mapped to "Father"
    geom_point(data = subset_df, aes(x = lon.x, y = lat.x, color = "Mother"), size = 3) +  # mother, mapped to "Mother"
    scale_color_manual(name = "Parent", values = c("Mother" = "blue", "Father" = "red")) +  # manual colour scale for the legend
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal()+
    theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title = element_text(size = 14)) + 
  annotate("text", x = 134.49543, y = 7.3136, label = paste(genotype), hjust = 1.1, vjust = -0.1, size = 8, color = "blue") # annotation for mother genotype
  
  return(p)
}

plots <- map(unique(join_df2$genotype.x), plot_genotype)
plots[[1]]
# to print each plot:
#walk(plots, print)

#save plots
walk2(plots, unique(join_df2$genotype.x), ~ggsave(path = "./plots/pairwise_crosses", 
                      filename = paste0("genotype_", .y, ".jpg"), plot = .x, width = 6, height = 4))

## Add to animation

img_dir <- "./plots/pairwise_crosses"  # Adjust this path
image_files <- list.files(path = img_dir, pattern = "*.jpg", full.names = TRUE)
images <- image_read(image_files)
gif <- image_animate(images, fps = 0.5)  # 

# Save the animated GIF
image_write(gif, path = "./plots/compiled_maps.gif")
image_write(gif, path = "./plots/compiled_maps_docent.gif")



#==========================================================================================================
#mum = x, sire = y
join_df3 = join_df2 %>% dplyr::select(genotype.x, genotype.y, dist_m, prop) %>%  drop_na()

# Calculate counts of distinct crosses for each genotype.y
cross_counts <- join_df3 %>%
  group_by(genotype.y) %>%
  summarise(distinct_crosses = n_distinct(genotype.x)) %>%
  arrange(desc(distinct_crosses))  %>%  data.frame()

join_df3 <- join_df3 %>%
  left_join(cross_counts, by = "genotype.y") %>%
  arrange(desc(distinct_crosses), genotype.y)


plot(density(cross_counts$distinct_crosses))
quan = quantile(cross_counts$distinct_crosses)
p3 <- ggplot(cross_counts, aes(x = distinct_crosses)) +
  geom_density(aes(fill = 'steelblue4'), alpha = 0.3) + 
  coord_cartesian(ylim = c(0.0, .75)) +
  scale_fill_manual(values = c("steelblue4", "white", "steelblue1", "white", "grey", "grey")) +
  scale_color_manual(values = c("steelblue4", "grey", "steelblue1", "steelblue4", "grey", "grey", "grey", "grey"))+
  tidybayes::stat_pointinterval(aes(y = 0, x = distinct_crosses), .width = c(.66, .95)) +
  #geom_errorbarh(aes(y = 0, xmin = lower_66, xmax = upper_66), height = 0.0, color = "black", linetype = "solid", size = 1.5) +
  #geom_errorbarh(aes(y = 0, xmin = lower_95, xmax = upper_95), height = 0.0, color = "black", linetype = "solid", size = .5)+
  scale_y_continuous(name = "Probability density")+
  scale_x_continuous(name = "Distinct pairwise crosses (no.)") +
  theme_sleek3()+
  geom_point(aes(x = unname(quan[3]), y = 0), color = "black", size = 3, shape = 21, fill = "black")
p3
#ggsave(p3, filename = 'dist_crosses.tiff',  path = "./plots", device = "tiff",  width = 6, height = 5)  #make sure to have .tiff on filename



# breeding units ----------------------------------------------------------

# Load required libraries
library(tidyverse)
library(tidygraph)
library(ggraph)
library(gridExtra)

# Reorganise the order of mothers
mother_list <- unique(join_df2$genotype.x)
mother_list <- sort(mother_list, method = "radix")  # Sort alphanumerically
mothers_with_c <- grep("^c", mother_list, value = TRUE)  # Mothers starting with 'c'
other_mothers <- setdiff(mother_list, mothers_with_c)    # Remaining mothers
sorted_mother_list <- c(mothers_with_c, other_mothers)   # Combine in the desired order

# Create spatial positions for each mother's breeding pairs
spatial_pairs <- map_df(sorted_mother_list, function(mother) {
  # Filter data for current mother
  mother_data <- join_df2 %>%
    filter(genotype.x == mother) %>%
    select(genotype.x, genotype.y, lat.x, lon.x, lat.y, lon.y, dist_m, angle_rad) %>%
    distinct()
  
  # Add mother position at origin
  mother_pos <- tibble(
    genotype.x = mother,
    genotype.y = mother,
    lat.x = mother_data$lat.x[1],  # Mother's latitude
    lon.x = mother_data$lon.x[1],  # Mother's longitude
    lat.y = mother_data$lat.x[1],  # Mother's latitude (matches origin)
    lon.y = mother_data$lon.x[1],  # Mother's longitude (matches origin)
    dist_m = 0,                    # Distance is 0 for self
    angle_rad = 0                  # Angle is 0 for self
  )
  
  bind_rows(mother_data, mother_pos)
})

# Create list of plots, one for each mother
breeding_plots <- map(sorted_mother_list, function(mother) {
  # Filter data for current mother
  current_data <- spatial_pairs %>%
    filter(genotype.x == mother)
  
  # Create edge data
  edges <- current_data %>%
    filter(genotype.y != mother) %>%
    select(from = genotype.x, to = genotype.y, lon.x, lat.x, lon.y, lat.y)
  
  # Create node data
  nodes <- current_data %>%
    select(name = genotype.y, lon = lon.y, lat = lat.y) %>%
    distinct()
  
  # Add the mother as a node
  mother_node <- tibble(
    name = mother,
    lon = current_data$lon.x[1],
    lat = current_data$lat.x[1]
  )
  nodes <- bind_rows(nodes, mother_node)
  
  # Create graph
  graph <- tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = FALSE
  )
  
  # Create plot
  sizee = 0.000015  # Make smaller for larger plots
  p <- ggraph(graph, layout = 'manual', x = nodes$lon, y = nodes$lat) +
    geom_edge_link(alpha = 0.5, color = "grey50") +
    geom_node_point(aes(color = name == mother), size = 2) +
    geom_node_text(aes(label = name), 
                   vjust = -0.3, 
                   hjust = 0.5, 
                   size = 3) +
    scale_color_manual(values = c("red", "blue")) +
    coord_fixed(xlim = c(min(nodes$lon) - sizee, max(nodes$lon) + sizee), 
                ylim = c(min(nodes$lat) - sizee, max(nodes$lat) + sizee)) +
    theme_void() +
    ggtitle(paste("Mother:", mother)) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 9),
      plot.margin = margin(0, 0, 0, 0, "mm")
    )
  
  return(p)
})

breeding_plots[[16]]

# Display all plots in a grid
# For 19 plots, we need a 4x5 layout
total_plots <- length(breeding_plots)  # 19
n_cols <- 4
n_rows <- 5

# Create the grid layout matrix, ensuring left-to-right row-wise order
layout <- t(matrix(seq_len(n_rows * n_cols), nrow = n_cols, ncol = n_rows))  # Transpose the matrix

# Add the plots to a list and fill remaining slots with NULL
plot_list <- c(breeding_plots,           # Your 19 plots
               replicate(n_rows * n_cols - total_plots, NULL))  # Fill remaining slots with NULL

# Create the grid arrangement
arranged_plots <- grid.arrange(
  grobs = plot_list,
  layout_matrix = layout,                # Use the transposed layout
  heights = unit(rep(0.19, n_rows), "npc"),  # Slightly reduced height per row
  widths = unit(rep(0.22, n_cols), "npc"),   # Slightly reduced width per column
  padding = unit(5, "mm")
  # Minimal padding between plots
)

 
# network sankey ----------------------------------------------------------

# Find counts of pairs
links <- join_df3 %>%
  dplyr::select(genotype.x, genotype.y) %>%
  count(genotype.x, genotype.y) %>%
  mutate(source = paste0(genotype.y, "_sire"),  # Add '_sire' to genotype.y (left side)
         target = paste0(genotype.x, "_dam")) %>%  # Add '_dam' to genotype.x (right side)
  dplyr::select(source, target, n) %>%
  arrange(desc(n), target)

# The unique node names ordered by distinct_crosses for sources (sire nodes)
sources_ordered <- cross_counts %>%
  mutate(source = paste0(genotype.y, "_sire")) %>%
  arrange(desc(distinct_crosses)) %>%
  pull(source)

# Order targets (dam nodes) by their first occurrence in links
targets_ordered <- links %>%
  arrange(target) %>%
  distinct(target) %>%
  pull(target)

# Combine ordered nodes
nodes <- data.frame(
  name = c(sources_ordered, targets_ordered) %>%
    unique()
)

# Reorder links$source and links$target to reflect the new node order
links <- links %>%
  mutate(IDsource = match(source, nodes$name) - 1,
         IDtarget = match(target, nodes$name) - 1) %>%
  arrange(IDsource, IDtarget)  # Ensure consistent ordering

# Colour by sire radial
links$group <- substr(links$source, 1, 1)
links$group <- as.factor(links$group)

# Define a colour scale for the links based on the group factor
link_colour_scale <- 'd3.scaleOrdinal()
  .domain(["group1", "group2", "group3"])    # Replace with actual levels of links$group
  .range(["#FF5733", "#33FF57", "#3357FF"])' # Replace with your desired colours

# Create the Sankey network
sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",  # Sire (source, left side)
  Target = "IDtarget",  # Dam (target, right side)
  Value = "n",
  NodeID = "name",
  units = "TWh",
  fontSize = 12,
  nodeWidth = 10,
  LinkGroup = "group",
  iterations = 0 # Disable dynamic node reordering to preserve manual order
)
sankey_plot

#ggsave(p4, filename = 'sankey_crosses.png',  path = "./plots", device = "png",  width = 6, height = 5)  #make sure to have .tiff on filename



# chord diagram -----------------------------------------------------------

# Load libraries
library(tidyverse)
library(circlize)

# Summarise the normalised weights using genotype.x (mother) and genotype.y (father)
cross_matrix <- join_df2 %>%
  group_by(genotype.x, genotype.y) %>%
  summarise(weight = sum(normalised_weight), .groups = "drop") %>%
  pivot_wider(names_from = genotype.y, values_from = weight, values_fill = 0) %>%
  column_to_rownames(var = "genotype.x") %>%
  as.matrix()

# Create a named colour vector for the grid.col argument
all_sectors <- c(rownames(cross_matrix), colnames(cross_matrix))
sector_colours <- setNames(rainbow(length(all_sectors)), all_sectors)

# Generate the chord diagram without the default labels
chordDiagram(
  cross_matrix, 
  transparency = 0.5,           # semi-transparent chords
  annotationTrack = "grid",     # draw only the grid (no default names)
  preAllocateTracks = list(track.height = 0.1), 
  grid.col = sector_colours,
  annotationTrackHeight = c(0.1, 0.02),
  link.sort = TRUE,  
  link.decreasing = FALSE
)

# Add custom labels
circos.trackPlotRegion(
  track.index = 1, 
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, 
      CELL_META$ylim[1] + 1, 
      CELL_META$sector.index, 
      facing = "clockwise", 
      niceFacing = TRUE, 
      cex = 0.6
    )
  },
  bg.border = NA
)



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

# Load necessary libraries
library(igraph)
library(dplyr)

# Summarize the data to count duplicates and get average distances
edge_counts <- join_df3 %>%
  group_by(genotype.x, genotype.y) %>%
  summarise(
    count = n(),                  # Number of times the pair occurs
    dist_m = mean(dist_m)         # Average distance
  ) %>%
  ungroup() %>%
  mutate(
    genotype.x = as.character(genotype.x),
    genotype.y = as.character(genotype.y)
  ) %>% data.frame()

# Create a list of unique genotypes for nodes
nodes <- unique(c(edge_counts$genotype.x, edge_counts$genotype.y))
nodes_df <- data.frame(name = nodes)

# Create the graph
g <- graph_from_data_frame(
  d = edge_counts,
  vertices = nodes_df,
  directed = FALSE
)

# Set edge attributes
E(g)$weight <- edge_counts$count
E(g)$dist_m <- edge_counts$dist_m

# Set edge width proportional to count
max_width <- 5
E(g)$width <- (E(g)$weight / max(E(g)$weight)) * max_width

# Initialize the distance matrix
genotypes <- nodes_df$name
n <- length(genotypes)
dist_mat <- matrix(Inf, nrow = n, ncol = n, dimnames = list(genotypes, genotypes))

# Fill the distance matrix
for (i in 1:nrow(edge_counts)) {
  from <- edge_counts$genotype.x[i]
  to <- edge_counts$genotype.y[i]
  dist <- edge_counts$dist_m[i]
  
  # Update matrix values
  dist_mat[from, to] <- dist
  dist_mat[to, from] <- dist  # Symmetric matrix
}

# Replace Inf with a large number for MDS
max_dist <- max(dist_mat[is.finite(dist_mat)])
dist_mat[is.infinite(dist_mat)] <- max_dist * 2

# Perform MDS to compute layout coordinates
mds_result <- cmdscale(dist_mat, k = 2, eig = TRUE)
layout_coords <- mds_result$points

# Plot the graph without labels
plot(
  g,
  layout = layout_coords,
  vertex.size = 5,
  vertex.label = NA,  # Disable default labels
  vertex.color = "skyblue",
  edge.width = E(g)$width,
  edge.color = 'grey',
  main = "Network"
)

# Add labels under the nodes
text(
  x = layout_coords[, 1],  # X-coordinates of nodes
  y = layout_coords[, 2] - 0.0,  # Y-coordinates of nodes (adjusted slightly below)
  labels = V(g)$name,  # Node labels
  cex = 0.8,           # Text size
  col = "black"        # Text color
)




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




