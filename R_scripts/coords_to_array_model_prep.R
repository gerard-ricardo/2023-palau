
## Todo

#BROKEN ATM. IN THE PROCESSES OF ADDING IN COORDS (AND COL SIZES) FROM META (CLEANED LA T LON)  DATA RATHER THAN GPX. THEN INDIVIDUAL COL SIZES 
#NEEDS TO BE ADDED

# - align better with flow direction. Based on releases needs to preference spoke 5. About 115 -120 deg
# - adjust flow speed. See marrote 10-12th for flow speed. Need to match with spawning time. Containers are 0.21m/s
#adjust dispersion constant based on dye
# - adjust colony size for each colony  (added to metadata now).
# - add uneven spawning times
# - there seems to be spoke 1 colonies, these may need to be added to centre patch

#Issues
# - sometimes crashes



# load libraries ----------------------------------------------------------
#install_github("envirometrix/plotKML")
library(devtools)
library(sf)
library(lattice)
library(plotKML)
library(RColorBrewer)
library(dplyr)


# inputs ------------------------------------------------------------------
resolution <- 0.1  #grid resolution
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
grid

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

# Extract spokes
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

# single colony
ID2 = "c19"
target <- coords_rotated %>% dplyr::filter(., ID == ID2)
#rand_center = list('xc13' = xc13, 'xc5' = xc5, 'xc10' = xc10,'xc19' = xc19)


## multiple colonies
# centres <- c("0089", "0093", "0098", "0101")  #manual random test
#target <- coords_rotated %>% dplyr::filter(., ID %in% all_colonies$ID)

## selected subset function
# myfun1 <- function(full_patch, subset_patch, colony_diam, resolution) {
#   coords1 <- full_patch
#   coords2 <- subset_patch
#   # Clear the grid to repopulate it
#   c(max(coords1$X), max(coords1$Y))
#   grid <- array(0, dim = c(10^ceiling(log10( max(coords1$X))), 10^ceiling(log10(max(coords1$Y)))*0.8))
#   dim(grid)
#   coords1$Y <- round((dim(grid)[2]/2 - max(coords1$Y)/2) +  coords1$Y  , 0)  #centre y's
#   coords1$X <-coords1$X + ((dim(grid)[1]) -   max(coords1$X)-2 )
#   coords1 <- coords1[coords1$ID %in% coords2$ID, ]
# 
#   # Repopulate the grid with rotated coordinates
#   for (i in 1:nrow(coords1)) {
#     # Convert to integer indices
#     x_idx <- (coords1[i, 1])
#     y_idx <- (coords1[i, 2])
#     grid[x_idx, y_idx] <- grid[x_idx, y_idx] + 1
#   }
# 
#   dim(grid)
#   # Plotting the rotated grid
#   # coul2 <- colorRampPalette(brewer.pal(8, "Purples"))(25)
#   # levelplot(t(apply(grid, 2, rev)),
#   #   col.regions = coul2, xlab = "Transverse",
#   #   ylab = "Longitudinal", main = "Conc. (cells/m^3)"
#   # )
# 
#   # save(grid, file = file.path("./Rdata", "2023_palau_fine_grid.RData"))
#   # load("./Rdata/2023_palau_fine_grid.RData") # grid
# 
#   # add colony border -------------------------------------------------------
# 
#   colony_diam_meters <- colony_diam / 100  # Convert cm to meters
#   buffer_radius <- (colony_diam_meters / 2) / resolution
#   (row_col_spe <- row(grid[, ])[which(!grid[, ] == 0)]) # finds all row vals
#   (col_col_spe <- col(grid[, ])[which(!grid[, ] == 0)])
#   df <- data.frame(x = row_col_spe, y = col_col_spe)
#   # Generate grid points, considering the resolution
#   grid_points <- expand.grid(x = seq(1, nrow(grid) * resolution, by = resolution),
#                              y = seq(1, ncol(grid) * resolution, by = resolution))
#   # find all cells within the buffer radius using hypot distance
#   points_in_buffer <- sapply(1:nrow(grid_points), function(i) {
#     any(sqrt((df$x - grid_points$x[i])^2 + (df$y - grid_points$y[i])^2) < buffer_radius)
#   })
#   grid1 <- grid
#   # Add 1 in all True grid points
#   for (i in 1:nrow(grid_points)) {
#     # Check if the point is within the buffer radius
#     if (points_in_buffer[i]) {
#       # Convert grid points to matrix indices and mark them
#       grid1[grid_points$x[i], grid_points$y[i]] <- 1
#     }
#   }
# 
#   # save(grid1, file = file.path("./Rdata", "2023pal_fine_border.RData"))
#   # load("./Rdata/2023pal_fine_border.RData") # grid1
# 
#   # plot using df
#   grid1_df <- expand.grid(X = 1:nrow(grid1), Y = 1:ncol(grid1))
#   grid1_df$Value <- as.vector(grid1)
#   grid1_df <- subset(grid1_df, Value != 0) # extract pos values
#   plot(grid1_df$X, grid1_df$Y,
#     col = ifelse(grid1_df$Value == 1, "red", "blue"), pch = 19,
#     main = "Grid Points Distribution", xlab = "Long", ylab = "Trans", xlim = c(1, dim(grid1)[1]), ylim = c(1, dim(grid1)[2])
#   )
# 
#   # matrix plot for subset (slow)
#   dim(grid1)
#   # coul <- colorRampPalette(brewer.pal(8, "Reds"))(25)
#   # high_res_plot <- levelplot(t(apply(grid1[1900:2300, 2800:3400], 2, rev)),
#   #                         col.regions = coul, xlab = "Transverse",
#   #                         ylab = "Longitudinal", main = "Conc. (cells/m^3)"
#   # )
#   # high_res_plot
# 
#   # convert back to 1 m grid ------------------------------------------------
#   # Dimensions of the original fine resolution buffer
#   res_back = 1/ resolution
#   finex_dim <- nrow(grid1) - 1 # Assuming square for simplicity
#   finey_dim <- nrow(grid1) - 1 # Assuming square for simplicity
#   # New grid dimensions for 1m resolution (assuming fine_dim is multiple of 100)
#   coarsex_dim <- finex_dim / res_back
#   coarsey_dim <- finex_dim / res_back
#   # Initialize new matrix for coarse resolution
#   grid1_coarse <- matrix(0, nrow = coarsex_dim, ncol = coarsey_dim)
#   # Function to aggregate 1 values in blocks, dynamically adjusting block size
#   block_size <- res_back # You can adjust this based on the actual dimensions of your matrices
#   aggregate_values <- function(mat, block_size) {
#     dim_x <- nrow(mat)
#     dim_y <- ncol(mat)
#     coarse_dim_x <- ceiling(dim_x / block_size) # new size
#     coarse_dim_y <- ceiling(dim_y / block_size)
#     coarse_mat <- matrix(0, nrow = coarse_dim_x, ncol = coarse_dim_y) # new blank grid
# 
#     for (i in seq(1, dim_x, by = block_size)) {
#       for (j in seq(1, dim_y, by = block_size)) {
#         end_i <- min(i + block_size - 1, dim_x)
#         end_j <- min(j + block_size - 1, dim_y)
#         # Sum the '1' values in the current block
#         block_sum <- sum(mat[i:end_i, j:end_j] == 1)
#         # Assign this sum to the corresponding cell in the coarse resolution matrix
#         coarse_mat[(i - 1) / block_size + 1, (j - 1) / block_size + 1] <- block_sum
#       }
#     }
#     return(coarse_mat)
#   }
#   grid1_coarse <- aggregate_values(grid1, block_size)
# 
#   # plot using df
#   grid2_df <- expand.grid(X = 1:nrow(grid1_coarse), Y = 1:ncol(grid1_coarse))
#   grid2_df$Value <- as.vector(grid1_coarse)
#   grid2_df <- subset(grid2_df, Value != 0) # extract pos values
#   plot(grid2_df$X, grid2_df$Y,
#     col = ifelse(grid2_df$Value == 1, "red", "blue"), pch = 19,
#     main = "Grid Points Distribution", xlab = "X", ylab = "Y",
#   )
# 
#   # grid plot
#   # levelplot(t(apply(grid1_coarse, 2, rev)),
#   #   col.regions = coul, xlab = "Transverse",
#   #   ylab = "Longitudinal", main = "Conc. (cells/m^3)"
#   # )
# 
#   grid_coarse <- grid1_coarse
#   return(grid_coarse)
# }

myfun1 <- function(full_patch, subset_patch, resolution) {
  coords1 <- full_patch
  coords2 <- subset_patch
  
  # Apply margin to avoid boundaries
  coords1$X <- coords1$X 
  coords1$Y <- coords1$Y
  
  # Initialize the grid based on the maximum coordinates
  offset = 20
  grid_dims_x <- 10^ceiling(log10(max(coords1$X))) + offset   #creates some room for offset below
  grid_dims_y <- 10^ceiling(log10(max(coords1$Y))) * 0.8
  grid <- array(0, dim = c(grid_dims_x, grid_dims_y))
  dim(grid)
  
  # Center the coordinates within the grid
  coords1$Y <- round((dim(grid)[2] / 2 - max(coords1$Y) / 2) + coords1$Y, 0)  # Center Y
  coords1$X <- coords1$X + (dim(grid)[1] - max(coords1$X) - 2) - offset # Adjust X
  coords1 <- coords1[coords1$ID %in% coords2$ID, ]  # Subset the coordinates
  
  # Repopulate the grid with the coordinates
  for (i in 1:nrow(coords1)) {
    x_idx <- coords1[i, "X"]
    y_idx <- coords1[i, "Y"]
    grid[x_idx, y_idx] <- grid[x_idx, y_idx] + 1
  }
  
  # Add colony borders using individual colony sizes
  # Prepare grid points only once to improve efficiency
  grid_points <- expand.grid(x = seq(1, nrow(grid), by = 1), y = seq(1, ncol(grid), by = 1))
  
  for (i in 1:nrow(coords1)) {
    fecun_zone = 0.7
    colony_diam <- coords1$colony_diam[i] *  fecun_zone # Use colony_diam from coords1
    if (is.na(colony_diam)) next  # Skip if colony_diam is NA
    
    colony_diam_meters <- colony_diam / 100  # Convert cm to meters
    buffer_radius <- (colony_diam_meters / 2) / resolution  # Calculate buffer radius
    
    x_idx <- coords1[i, "X"]
    y_idx <- coords1[i, "Y"]
    
    # Calculate distances from the current colony center to all grid points
    distances <- sqrt((grid_points$x - x_idx)^2 + (grid_points$y - y_idx)^2)
    
    # Identify grid points within the buffer radius
    points_in_buffer <- distances < buffer_radius
    
    # Mark these points in the grid
    grid[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] <- 1
  }
  
  # Convert the fine grid to a coarser resolution
  # Define the block size for aggregation
  res_back <- 1 / resolution  # Assuming you want to aggregate back to 1m resolution
  
  # Calculate new grid dimensions
  coarse_dim_x <- floor(nrow(grid) / res_back)
  coarse_dim_y <- floor(ncol(grid) / res_back)
  
  # Initialize the coarse grid
  grid_coarse <- matrix(0, nrow = coarse_dim_x, ncol = coarse_dim_y)
  
  # Aggregate the fine grid into the coarse grid
  for (i in 1:coarse_dim_x) {
    for (j in 1:coarse_dim_y) {
      # Define the indices for the fine grid
      x_indices <- ((i - 1) * res_back + 1):(i * res_back)
      y_indices <- ((j - 1) * res_back + 1):(j * res_back)
      
      # Sum the values in the fine grid block
      block_sum <- sum(grid[x_indices, y_indices])  
      block_sum2 = block_sum/(1/resolution^2) * 10000  #this should get the proportion of the total area and convert to cm (100 x 100 = 10000)
      
      # Assign the aggregated value to the coarse grid
      grid_coarse[i, j] <- block_sum2
    }
  }
  
  return(grid_coarse)
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

positive_values <- grid_coarse[grid_coarse > 0]
hist(positive_values)

# saved maps
#save_path <- file.path("./Rdata", paste0("2023palau_coarse_grid_", ID2, ".RData"))
#save(grid_coarse, file.path("./Rdata", paste0("2023palau_coarse_grid_", ID2, ".RData")))
#save(grid_coarse, file = file.path("./Rdata", paste0("2023palau_coarse_grid_", ID2, ".RData")))

save(grid_coarse, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata", 
                                   "2023palau_coarse_grid_xc19.RData"))
#save(grid_coarse, file = file.path("./Rdata", "2023palau_coarse_grid_rand_centres.RData"))  #grid_coarse


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



