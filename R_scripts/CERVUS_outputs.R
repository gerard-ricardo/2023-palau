

# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
library(ggmap)
library(units)
library(sf)
library(ggmap)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code


# 1 Import cervus out data -----------------------------------------------------------
data1 <- read.csv(file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/4 Palau genetics mixing/Cervus", "parentage_out1.csv"))



# 2. Data prep -----------------------------------------------
str(data1) # check data type is correct
data1$offsp_id <- as.factor(as.character(data1$Offspring.ID))
data1$moth_id <- as.factor(as.character(data1$Mother.ID))
data1$fath_id <- as.factor(as.character(data1$Candidate.father.ID))
data2 = data1 %>% select(., c(offsp_id, moth_id, fath_id))
data2 = data2 %>% mutate(id = moth_id)


# 1 Import cervus out data -----------------------------------------------------------
meta = data_gl_filtered@other$ind.metrics
meta = meta  %>% rename(id2 = id) %>% mutate(id = paste0("X", id2)) 
meta2 <- meta %>% select(c(id, lat, lon))

# join tables -------------------------------------------------------------
join_df1 <- left_join(data2, meta2, by = 'id')  #add mother lat lon
join_df1 <- join_df1 %>% select(-id) %>% mutate(id = fath_id)  #reset id to fathers
join_df2 <- left_join(join_df1, meta2, by = 'id')  #add father lat lon



# calculate distances -----------------------------------------------------

# Convert the data frame to two separate sf objects with WGS 84 CRS
points_x <- st_as_sf(join_df2, coords = c("lon.x", "lat.x"), crs = 4326)
points_y <- st_as_sf(join_df2, coords = c("lon.y", "lat.y"), crs = 4326)

# Using UTM zone 53N for Palau (EPSG: 32653)
points_x_proj <- st_transform(points_x, crs = 32653)
points_y_proj <- st_transform(points_y, crs = 32653)

join_df2$dist <- st_distance(points_x_proj, points_y_proj, by_element = TRUE)
join_df2$dist_m = as.numeric(join_df2$dist )

# analyses ----------------------------------------------------------------
hist(join_df2$dist_m)
plot(density(join_df2$dist_m))
quan = quantile(join_df2$dist_m)
unname(quan[3])
selfing_no = length(which(join_df2$dist_m == 0))  #self fertilisation
selfing_prop = selfing_no / nrow(join_df2)
# Units: [m]
# 0%       25%       50%       75%      100% 
# 0.000000  1.569776  3.455444 11.996440 30.690024 
# density
p1 <- ggplot(join_df2, aes(x = dist_m)) +
  geom_density(aes(fill= 'steelblue4'), alpha = 0.3) +
  tidybayes::stat_pointinterval(aes(y = 0, x = dist_m), .width = c(.66, .95))
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

map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "color", maptype = "satellite")
ggmap(map) +
  #geom_segment(join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  labs(x = "Longitude", y = "Latitude") 
ggsave(filename = '2023palau_design.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf


#zoomed out
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
ggmap(map) +
  geom_segment(join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
geom_point(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3)+
geom_point(data = join_df2, aes(x = lon.y, y = lat.y), color = "blue", size = 3) +
  labs(x = "Longitude", y = "Latitude") 
ggsave(filename = '2023palau_gen_zoomout.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf


#zoomed in
map <- get_googlemap(center = c(lon = join_df2$lon.x[1], lat = join_df2$lat.x[1]), zoom = 22, color = "bw", maptype = "satellite")
ggmap(map) +
  geom_segment(join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3)+
  geom_point(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3)+
  geom_point(data = join_df2, aes(x = lon.y, y = lat.y), color = "blue", size = 3) +
  labs(x = "Longitude", y = "Latitude") 
ggsave(filename = '2023palau_gen_zoomin.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf

  
