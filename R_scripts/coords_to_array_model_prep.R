#This script applies the corrected GPS coordinates from the metadata (gneetics) onto a fine grid, rotates it, and then applys a perimeter
# to the centre and infils each overlapping cell. Cells (should overlap if two colonies fall in the same cell since in reality no 
#colonies overlapped). Everything is then converted back to a 1 metre grid. The output is fecund cm^2 at present. 

#How to run (currently)
# - resolution at 0.01
# - spemap should be all colonies minus the same spoke if running a spoke. If running a centre, you might have to map spemap up
# individually, but at the moment running all together which includes selffert bias


## Todo

#Note: discretisation seems to increase fecund cm^2 above actual by a fair bit. 0.01 better than 0.1

#NEEDS TO BE ADDED

# - align better with flow direction. Based on releases needs to preference spoke 5. About 115 -120 deg
# - adjust flow speed. See marrote 10-12th for flow speed. Need to match with spawning time. Containers are 0.21m/s
#adjust dispersion constant based on dye
# - add uneven spawning times
# - there seems to be spoke 1 colonies, these may need to be added to centre patch

#Issues
# - sometimes crashes



# load libraries ----------------------------------------------------------
#install_github("envirometrix/plotKML")
#library(devtools)
library(sf)
library(lattice)
library(plotKML)
library(RColorBrewer)
library(dartR)


# inputs ------------------------------------------------------------------
resolution <- 0.01  #grid resolution (only works for 0.01 and 0.1
#colony_diam <- 44  #in cm


# import gmx --------------------------------------------------------------
#gpx_file <- "C:/Users/gerar/OneDrive/1_Work/3_Results/11 Allee effects/3 field experiments/2023_Palau/waypoints/Waypoints_2023-04-04.gpx"
# gpx_file <- "C:/Users/gerar/OneDrive/1_Work/3_Results/11 Allee effects/3 field experiments/2022_12 Heron/13_12_22/Waypoints_13-DEC-22.gpx"
# gpx_data <- readGPX(gpx_file)
# coords <- data.frame(x = gpx_data$waypoints$lon, y = gpx_data$waypoints$lat, ID = gpx_data$waypoints$name)
# save(coords, file = file.path("./Rdata", "2022_heron_coords.RData"))

# add to grid -------------------------------------------------------------
#load("./Rdata/2023_palau_coords.RData") # coords
#load("./Rdata/2022_heron_coords.RData") # coords

#import from sequencing metadata
load("./Rdata/data_gl.RData")  #data_gl_filtered.RData
meta = data_gl@other$ind.metrics
meta = meta  %>% rename(id2 = id) %>% mutate(id = paste0("X", id2)) 
meta2 <- meta %>% dplyr::select(c(lat, lon, genotype, total_mean_dia ))
#meta2 <- meta2 %>% mutate(genotype = gsub("_5$", "_05", genotype))
coords <- meta2 %>% rename(x = lon, y = lat, ID = genotype, colony_diam = total_mean_dia) %>% filter(complete.cases(.)) %>% 
  distinct(ID, .keep_all = TRUE)
head(coords)


# Convert geographic coordinates to local grid coordinates (very simplistic approach)
# Assuming 1 degree of latitude = 111 km, and similarly for longitude at this latitude
# coords_converted <- transform(coords,
#                               x = (x - min(x)) * 111e3,
#                               y = (y - min(y)) * 111e3)
coords_sf <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)

#palau
utm_zone <- 52
crs_utm <- st_crs(sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs", utm_zone))

# #heron
# utm_zone <- 56
# crs_utm <- st_crs(sprintf("+proj=utm +zone=%d +south +datum=WGS84 +units=m +no_defs", utm_zone))

coords_utm <- st_transform(coords_sf, crs = crs_utm)
coords_converted <- as.data.frame(st_coordinates(coords_utm))

# Calculate the dimensions of the grid
grid_dims <- ceiling(as.numeric((sapply(coords_converted, max) - sapply(coords_converted, min))) / resolution) + 1
# 534 513

# Create the 2D array with initial values of 0, with a small buffer
grid <- array(0, dim = c(10^ceiling(log10(grid_dims[1])), 10^ceiling(log10(grid_dims[2]))))

# Convert coordinates to grid indices
coords_idx <- as.data.frame(apply(coords_converted, 2, function(col) ceiling((col - min(col)) / resolution) + 1))
min(coords_idx$X)
min(coords_idx$Y)
coords_idx$ID <- coords$ID
coords_idx$colony_diam <- coords$colony_diam

# Populate the grid
for (i in 1:nrow(coords_idx)) {
  grid[coords_idx[i, "X"], coords_idx[i, "Y"]] <- grid[coords_idx[i, "X"], coords_idx[i, "Y"]] + 1
}

#View(grid)

plot(coords_idx$X, coords_idx$Y, main = " Points Distribution", xlab = "X", ylab = "Y", pch = 19)
text(coords_idx$X, coords_idx$Y, labels = coords_idx$ID, pos = 3, cex = 0.7)

## uses large memory
# coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
# levelplot(t(apply(grid[, ], 2, rev)),
#   col.regions = coul2, xlab = "Transverse",
#   ylab = "Longitudinal", main = "Conc. (cells/m^3)"
# )

# rotate around centre colony 0101, and extract centre ----------------------------------------
# Find the rotation center (coordinates of ID '0101')
rotation_center <- coords_idx[coords_idx$ID == "c17", c("X", "Y")]  #palau
#rotation_center <- coords_idx[coords_idx$ID == "045", c("X", "Y")]  #heron


# Convert 135 degrees to radians for the rotation
theta <- 130 * (pi / 180)  #palau
#theta <- -50 * (pi / 180)  #heron
# Rotate all points around the rotation center
# coords_rotated <- t(apply(coords_idx[, c("X", "Y")], 1, function(coord) {
#   x <- coord[1]
#   y <- coord[2]
#   x_c <- rotation_center$X
#   y_c <- rotation_center$Y
#   # Apply rotation equations
#   x_new <- round(cos(theta) * (x - x_c) - sin(theta) * (y - y_c) + x_c, 0)
#   y_new <- round(sin(theta) * (x - x_c) + cos(theta) * (y - y_c) + y_c, 0)
#   return(c(x_new, y_new))
# }))

rotate_coords <- function(coords_df, center, angle_rad) {
  coords_shifted <- coords_df[, c("X", "Y")] - matrix(unlist(center), nrow = nrow(coords_df), ncol = 2, byrow = TRUE)
  rotation_matrix <- matrix(c(cos(angle_rad), -sin(angle_rad), sin(angle_rad), cos(angle_rad)), nrow = 2, byrow = TRUE)
  coords_rotated <- as.matrix(coords_shifted) %*% rotation_matrix
  coords_rotated <- coords_rotated + matrix(unlist(center), nrow = nrow(coords_df), ncol = 2, byrow = TRUE)
  coords_rotated <- round(coords_rotated, 0)
  coords_rotated_df <- data.frame(X = coords_rotated[, 1], Y = coords_rotated[, 2], coords_df[, !(names(coords_df) %in% c("X", "Y"))])
  coords_rotated_df$X = coords_rotated_df$X + abs(min(coords_rotated_df$X)) #adjust to pos values
  coords_rotated_df$Y = coords_rotated_df$Y + abs(min(coords_rotated_df$Y)) #adjust to pos values
  return(coords_rotated_df)
}


coords_rotated <- rotate_coords(coords_df = coords_idx, center = rotation_center, angle_rad = theta)


# Adjust coords_rotated to fit in the grid dimensions if necessary
# coords_rotated <- apply(coords_rotated, 2, function(col) {
#   round(col - min(col) + 1, 0)
# })
str(coords_rotated)
coords_rotated <- coords_rotated %>% data.frame()
#colnames(coords_rotated) <- c("X", "Y") # replace Y name
coords_rotated$ID <- coords_idx$ID
plot(coords_rotated$X, coords_rotated$Y, main = "Rotated Points Distribution", xlab = "X", ylab = "Y", pch = 19)
text(coords_rotated$X, coords_rotated$Y, labels = coords_rotated$ID, pos = 3, cex = 0.7)
unique(coords_rotated)


# extract spe and egg maps ------------------------------------------------

coords_rotated$ID
spoke_1 <- c("1_05", "1_3", "1_10", "1_30")
spoke_2 <- c("2_05", "2_3", "2_10", "2_30")
spoke_3 <- c("3_05", "3_3", "3_10", "3_30")
spoke_4 <- c("4_05", "4_3", "4_10", "4_30")
spoke_5 <- c("5_05", "5_3", "5_10", "5_30")
spoke_6 <- c("6_05", "6_3", "6_10", "6_30")
spoke_7 <- c("7_05", "7_3", "7_10", "7_30")
spoke_8 <- c("8_05", "8_3", "8_10", "8_30")
all_spokes <- c(spoke_1, spoke_2, spoke_3, spoke_4, spoke_5, spoke_6, spoke_7, spoke_8)
spokes <- coords_rotated %>% dplyr::filter(., ID %in% all_spokes)
centre <- coords_rotated %>% dplyr::filter(!(ID %in% all_spokes))
all_colonies = coords_rotated

## to map a single colony
ID2 = "c19"
target <- coords_rotated %>% dplyr::filter(., ID == ID2)
#'c13' 'c5' 'c10' 'c19'   #these are the random centre


## to map all other or multiple colonies
# manual_centres <- c("0089", "0093", "0098", "0101")  #manual random test
#target <- coords_rotated %>% dplyr::filter(., ID %in% all_colonies$ID)  #all colonies for spemap

# find accurate fecund area
target$fecund_area = pi * (target$colony_diam / 2)^2 * 0.7
sum(target$fecund_area)


myfun1 <- function(full_patch, subset_patch, resolution) {
  coords1 <- full_patch
  coords2 <- subset_patch
  
  
  # Initialize the grid based on the maximum coordinates
  offset_m <- 2 # desired physical offset, e.g. 2 m
  offset <- ceiling(offset_m / resolution) # convert metres to number of fine-grid cells
  grid_dims_x <- 10^ceiling(log10(max(coords1$X))) + offset
  grid_dims_y <- (10^ceiling(log10(max(coords1$Y))) + offset) * 0.8  #0.8 makes grid a bit thinner
  grid <- array(0, dim = c(grid_dims_x, grid_dims_y))
  grid_egg <- array(0, dim = c(grid_dims_x, grid_dims_y))  # define egg grid
  
  dim(grid)
  # Center the coordinates within the grid
  coords1$Y <- round((dim(grid)[2] / 2 - max(coords1$Y) / 2) + coords1$Y, 0)  # Center Y
  coords1$X <- coords1$X + (dim(grid)[1] - max(coords1$X) - 2) - offset # Adjust X
  
  # offset_m <- 2 # 2-metre offset
  # 
  # # 1. Define a real bounding box with offset in metres
  # x_min <- min(coords1$X) - offset_m
  # x_max <- max(coords1$X) + offset_m
  # y_min <- min(coords1$Y) - offset_m
  # y_max <- max(coords1$Y) + offset_m
  # 
  # # 2. Convert that bounding box to integer grid dimensions
  # grid_dims_x <- (x_max - x_min) / resolution
  # grid_dims_y <- (y_max - y_min) / resolution
  # grid <- array(0, dim = c(ceiling(grid_dims_x), ceiling(grid_dims_y)))
  # 
  # # 3. Shift your colony coordinates so that x_min/y_min becomes (1,1)
  # coords1$X <- coords1$X - x_min
  # coords1$Y <- coords1$Y - y_min
  
  
  
  spemap <- coords1[coords1$ID %in% coords2$ID, ]  # Extract not target from full patch
  #eggmap <- coords1[coords1$ID %in% coords2$ID, ]  # Extract target from full patch
  
  
  # Repopulate the grid with the coordinates
  for (i in 1:nrow(spemap)) {
    x_idx <- spemap[i, "X"]
    y_idx <- spemap[i, "Y"]
    grid[x_idx, y_idx] <- grid[x_idx, y_idx] + 1
  }
  
  # for (i in 1:nrow(eggmap)) {
  #   x_idx <- eggmap[i, "X"]  # eggmap X
  #   y_idx <- eggmap[i, "Y"]  # eggmap Y
  #   grid_egg[x_idx, y_idx] <- grid_egg[x_idx, y_idx] + 1  # increment egg map
  # }
  
  n_finecells_colony <- numeric(nrow(spemap)) # all zero initially
  #n_finecells_colony_egg <- numeric(nrow(eggmap))  # all zero initially
  
  
  # Add colony borders using individual colony sizes
  # Prepare grid points only once to improve efficiency
  grid_points <- expand.grid(x = seq(1, nrow(grid), by = 1), y = seq(1, ncol(grid), by = 1))
  
  for (i in 1:nrow(spemap)) {
    fecun_zone = 0.7 # fraction of area
    colony_diam <- spemap$colony_diam[i] # keep original diameter
    if (is.na(colony_diam)) next # skip if colony_diam is NA
    colony_diam_meters <- (colony_diam / 100) * sqrt(fecun_zone) # apply sqrt(0.7) to radius. Sqt adjusts the diameter properly.
    buffer_radius <- (colony_diam_meters / 2) / resolution # final radius
    x_idx <- spemap[i, "X"]
    y_idx <- spemap[i, "Y"]
    # Calculate distances from the current colony center to all grid points
    distances <- sqrt((grid_points$x - x_idx)^2 + (grid_points$y - y_idx)^2)
    # Identify grid points within the buffer radius
    points_in_buffer <- distances < buffer_radius
    n_finecells_colony[i] <- sum(points_in_buffer) #count fine cells
    # Mark these points in the grid
    #grid[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] <- 1
    grid[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] <-
      grid[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] + 1  #allow overlapping points to sum
  }
  
  # for (i in 1:nrow(eggmap)) {
  #   fecun_zone = 0.7  # fraction of area
  #   colony_diam <- eggmap$colony_diam[i]  # diameter from eggmap
  #   if (is.na(colony_diam)) next  # skip if colony_diam is NA
  #   colony_diam_meters <- (colony_diam / 100) * sqrt(fecun_zone)  # scaled diameter
  #   buffer_radius <- (colony_diam_meters / 2) / resolution  # final radius
  #   x_idx <- eggmap[i, "X"]
  #   y_idx <- eggmap[i, "Y"]
  #   distances <- sqrt((grid_points$x - x_idx)^2 + (grid_points$y - y_idx)^2)
  #   points_in_buffer <- distances < buffer_radius
  #   n_finecells_colony_egg[i] <- sum(points_in_buffer)
  #   grid_egg[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] <- 
  #     grid_egg[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] + 1  # allow overlap
  # }
  
  message("Fine-grid cells for each colony:")
  print(data.frame(
    ID = spemap$ID,
    Diam_cm = spemap$colony_diam,
    FineCells = n_finecells_colony,
    Actual = pi * (spemap$colony_diam / 2)^2  * 0.7
  ))
  
  # Convert the fine grid to a coarser resolution
  # Define the block size for aggregation
  res_back <- 1 / resolution  # Aggregate back to 1m resolution
  
  # Calculate new grid dimensions
  coarse_dim_x <- floor(nrow(grid) / res_back)
  coarse_dim_y <- floor(ncol(grid) / res_back)
  
  # Initialize the coarse grid
  grid_coarse_sperm <- matrix(0, nrow = coarse_dim_x, ncol = coarse_dim_y)
  #grid_coarse_egg <- matrix(0, nrow = coarse_dim_x, ncol = coarse_dim_y)  # initialise coarse egg grid
  
  
  # Aggregate the fine grid into the coarse grid
  for (i in 1:coarse_dim_x) {
    for (j in 1:coarse_dim_y) {
      # Define the indices for the fine grid
      x_indices <- ((i - 1) * res_back + 1):(i * res_back)
      y_indices <- ((j - 1) * res_back + 1):(j * res_back)
      
      # Sum the values in the fine grid block
      block_sum <- sum(grid[x_indices, y_indices])  
      block_sum2 = block_sum/(1/resolution^2) * 100 * 100  #this should get the proportion of the total area and convert to cm (100 x 100 = 10000)
      
      # Assign the aggregated value to the coarse grid
      grid_coarse_sperm[i, j] <- block_sum2
    }
  }
  
  # for (i in 1:coarse_dim_x) {
  #   for (j in 1:coarse_dim_y) {
  #     x_indices <- ((i - 1) * res_back + 1):(i * res_back)
  #     y_indices <- ((j - 1) * res_back + 1):(j * res_back)
  #     block_sum <- sum(grid_egg[x_indices, y_indices])  
  #     block_sum2 = block_sum/(1/resolution^2) * 100 * 100  # convert to cm² scale
  #     grid_coarse_egg[i, j] <- block_sum2
  #   }
  # }
  return(grid_coarse_sperm)
  #return(list(grid_coarse_sperm = grid_coarse_sperm, grid_coarse_egg = grid_coarse_egg))  # return both grids
}


# apply function
#out <-  myfun1(full_patch = coords_rotated, subset_patch = target, colony_diam = colony_diam, resolution = resolution) #centre, target
out <- myfun1(full_patch = coords_rotated, subset_patch = target, resolution = resolution)


grid_coarse <-  out
# grid plot
coul <- colorRampPalette(brewer.pal(8, "Purples"))(25)
levelplot(t(apply(grid_coarse, 2, rev)),
  col.regions = coul, xlab = "Transverse",
  ylab = "Longitudinal", main = "Conc. (cells/m^3)"
)
# checking

sum(grid_coarse)
sum(grid_coarse) - sum(target$fecund_area)  #should be zero
target
#View(grid_coarse)
(which_cells <- which(grid_coarse > 1, arr.ind = TRUE)) # get row & column of colonies
grid_coarse[which_cells] # check actual values
dim(grid_coarse)



# saved maps --------------------------------------------------------------

save(grid_coarse, file = file.path("./Rdata", paste0("2023palau_coarse_grid_", ID2, ".RData")))

# save(grid_coarse, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata", 
#                                    "2023palau_coarse_grid_xc19.RData"))

#save(grid_coarse, file = file.path("./Rdata", "2023palau_coarse_grid_rand_centres.RData"))  #grid_coarse

#spemap
# save(grid_coarse, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata", 
#                                    "2023palau_coarse_grid_centre.RData")) #this is what you load into the config file as spemap



# bind eggmaps to list ----------------------------------------------------

load("./Rdata/2023palau_coarse_grid.RData") # grid_coarse
load("./Rdata/2023palau_coarse_grid_centre.RData") # grid_coarse
load("./Rdata/2023palau_coarse_grid_spokes.RData") # grid_coarse
load("./Rdata/2023palau_coarse_grid_rand_centres.RData") # grid_coarse

#random_centre
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc13.RData") # grid_coarse
xc13 = grid_coarse 
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc5.RData") # grid_coarse
xc5 = grid_coarse 
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc10.RData") # grid_coarse
xc10 = grid_coarse 
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc19.RData") # grid_coarse
xc19 = grid_coarse 
rand_center = list('xc13' = xc13, 'xc5' = xc5, 'xc10' = xc10,'xc19' = xc19)
# save(rand_center, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata",
#                                    "2023palau_rand_center_list.RData"))
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_rand_center_list.RData") # grid_coarse


#spke2
load("./Rdata/2023palau_coarse_grid_0106.RData") # grid
x0106 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0107.RData") # grid
x0107 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0108.RData") # grid
x0108 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0109.RData") # grid
x0109 = grid_coarse 
spoke2 = list('x0106' = x0106, 'x0107' = x0107, 'x0108' = x0108,'x0109' = x0109)
#save(spoke2, file = file.path("./Rdata", "2023palau_spoke2_list.RData"))
load("./Rdata/2023palau_spoke2_list.RData") # grid_coarse

#spke3
load("./Rdata/2023palau_coarse_grid_0110.RData") # grid
x0110 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0111.RData") # grid
x0111 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0112.RData") # grid
x0112 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0113.RData") # grid
x0113 = grid_coarse 
spoke3 = list('0110' = x0110, '0111' = x0111, '0112' = x0112, '0113' = x0113)
save(spoke3, file = file.path("./Rdata", "2023palau_spoke3_list.RData"))
load("./Rdata/2023palau_spoke3_list.RData") # spoke3

#spoke4
load("./Rdata/2023palau_coarse_grid_0114.RData") # grid
x0114 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0115.RData") # grid
x0115 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0116.RData") # grid
x0116 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0117.RData") # grid
x0117 = grid_coarse
spoke4 = list('0114' = x0114, '0115' = x0115, '0116' = x0116, '0117' = x0117)
save(spoke4, file = file.path("./Rdata", "2023palau_spoke4_list.RData"))
load("./Rdata/2023palau_spoke4_list.RData") # grid_coarse

#spkoe5
load("./Rdata/2023palau_coarse_grid_0118.RData") # grid
x0118 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0119.RData") # grid
x0119 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0120.RData") # grid
x0120 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0121.RData") # grid
x0121 = grid_coarse
spoke5 = list('0118' = x0118, '0119' = x0119, '0120' = x0120, '0121' = x0121)
save(spoke5, file = file.path("./Rdata", "2023palau_spoke5_list.RData"))
load("./Rdata/2023palau_spoke5_list.RData") # grid_coarse

#spoke6
load("./Rdata/2023palau_coarse_grid_0122.RData") # grid
x0122 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0123.RData") # grid
x0123 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0124.RData") # grid
x0124 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0125.RData") # grid
x0125 = grid_coarse
spoke6 = list('0122' = x0122, '0123' = x0123, '0124' = x0124, '0125' = x0125)
save(spoke6, file = file.path("./Rdata", "2023palau_spoke6_list.RData"))
load("./Rdata/2023palau_spoke6_list.RData") # grid_coarse

#spoke 7
load("./Rdata/2023palau_coarse_grid_0126.RData") # grid
x0126 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0127.RData") # grid
x0127 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0128.RData") # grid
x0128 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0129.RData") # grid
x0129 = grid_coarse
spoke7 = list('0126' = x0126, '0127' = x0127, '0128' = x0128, '0129' = x0129)
save(spoke7, file = file.path("./Rdata", "2023palau_spoke7_list.RData"))
# load("./Rdata/2023palau_spoke7_list.RData") # grid_coarse



