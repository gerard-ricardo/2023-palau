#Panel 3 - Map and experimental design?

#pairwise_dist_plot
#polar_plot
#arranged_plots
#sankey_plot
#distance_map_plot

class(pairwise_dist_plot)
class(polar_plot)
class(arranged_plots)
class(sankey_plot)
class(distance_map_plot)

library(gridExtra)
library(grid)
library(plotrix)
library(patchwork)
library(cowplot)
library(htmlwidgets)
library(webshot)
library(png)

# conversion to grobs -----------------------------------------------------
# For sankey_plot (htmlwidget)

#webshot::install_phantomjs()

# Save the Sankey plot as HTML and convert to PNG
# htmlwidgets::saveWidget(sankey_plot, "temp_sankey.html")
# webshot("temp_sankey.html", "temp_sankey.png")
# sankey_plot_grob <- rasterGrob(readPNG("temp_sankey.png"))
# # Clean up temporary files
# file.remove(c("temp_sankey.html", "temp_sankey.png"))
# class(sankey_plot_grob)

# Save the Sankey plot as HTML with adjusted height
htmlwidgets::saveWidget(
  sankey_plot, 
  "temp_sankey.html",
  selfcontained = TRUE
)

# Take screenshot with specific dimensions and no extra space
webshot(
  "temp_sankey.html", 
  "temp_sankey.png",
  selector = ".sankeyNetwork",  # Only capture the actual plot
  zoom = 3,                     # Higher resolution
  delay = 0.2,                 # Small delay to ensure plot renders
  vwidth = 500,                # Viewport width
  vheight = 600               # Viewport height - adjust this to reduce white space
)

# Read and convert to grob, maintaining aspect ratio
png_data <- readPNG("temp_sankey.png", native = TRUE)
sankey_plot_grob <- rasterGrob(
  png_data,
  interpolate = TRUE,
  width = unit(0.8, "npc"),
  height = unit(1, "npc")
)

# Clean up temporary files
file.remove(c("temp_sankey.html", "temp_sankey.png"))

##########################
## polar plot
png(filename = "temp1.png", width = 400, height = 400)
create_polar_plot()
dev.off()
# Read back as a raster and convert to grob
temp_png <- readPNG("temp1.png")
polar_plot_grob <- rasterGrob(temp_png)
# Clean up temporary file
file.remove("temp1.png")
class(polar_plot_grob)


# inset -------------------------------------------------------------------
# First create the simplified inset plot
# inset_pairwise <- pairwise_dist_plot + 
#   #theme_minimal() +
#   theme(
#     plot.background = element_rect(fill = "white", color = "black"),
#     plot.margin = margin(1, 1, 1, 1),
#     text = element_text(size = 8)  # Smaller text for inset
#   )
# 
# # Get the coordinate ranges from your original data
# x_range <- range(c(join_df2$lon.x, join_df2$lon.y))
# y_range <- range(c(join_df2$lat.x, join_df2$lat.y))
# 
# # Calculate inset position (bottom right corner)
# inset_width <- diff(x_range) * 0.65   # Width of inset (30% of plot width)
# inset_height <- diff(y_range) * 0.8  # Height of inset (30% of plot height)
# 
# # Combine main plot with inset
# map_with_inset <- distance_map_plot +
#   annotation_custom(
#     grob = ggplotGrob(inset_pairwise),
#     #edge positon fo plot
#     xmin = 134.49560,     # Move left edge further left
#     xmax = 134.49582,     # Keep right edge near main plot edge
#     ymin = 7.313051,     # Lift bottom edge slightly
#     ymax = 7.313229      # Make height larger
#   )
# map_with_inset


# panel -------------------------------------------------------------------


# panel2 <- grid.arrange(
#   # Left column: 3 plots stacked vertically
#   arrangeGrob(distance_map_plot, pairwise_dist_plot, sankey_plot_grob, ncol = 1),
#   # Right column: 2 plots stacked vertically
#   arrangeGrob(arranged_plots, polar_plot_grob, ncol = 1),
#   ncol = 2,  # Two columns: left and right
#   widths = c(1, 1)  # Adjust width proportions (left column twice as wide as the right)
# )

# Create a layout matrix
# Each number represents a plot position
# Numbers that appear multiple times mean that plot spans multiple cells
lay <- rbind(
  c(1, 1, 1, 2, 2),    # First row: map_with_inset(1) | arranged_plots(2)
  c(1, 1, 1, 2, 2),  
  c(5, 5, 5, 2, 2),# Second row: map_with_inset cont. | arranged_plots cont.
  c(3, 3, 3 ,4, 4),    # Third row: sankey_plot(3) | arranged_plots cont.
  c(3, 3, 3, 4, 4)     # Fourth row: polar_plot(4) | last_plot(5)
)

# Create list of plots in the order they appear in the layout
plot_list <- list(
  distance_map_plot,      # 1
  arranged_plots,      # 2
  sankey_plot_grob,    # 3
  polar_plot_grob,      # 4 (assuming this exists)
  pairwise_dist_plot   #5
)

# Create the final panel
panel2 <- grid.arrange(
  grobs = plot_list,
  layout_matrix = lay,
  widths =  c(1, 1, 1, 1, 1),      # Four columns need four width values
  heights = c(1, 1, 1, 1, 1)      # Four rows need four height values
)

#ggsave(panel2, filename = 'fig2.pdf',  path = "./plots", device = 'pdf',  width = 9, height = 9.5)  #




# will try saving sepeately for Illustrator -------------------------------

# Create plots directory if it doesn't exist
plot_dir <- "plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Save distance map plot
ggsave(
  file.path(plot_dir, "distance_map_plot.pdf"),
  plot = distance_map_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Save arranged plots (model effects)
ggsave(
  file.path(plot_dir, "arranged_plots.pdf"),
  plot = arranged_plots,
  width = 4,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Save pairwise distance plot
ggsave(
  file.path(plot_dir, "pairwise_dist_plot.pdf"),
  plot = pairwise_dist_plot,
  width = 6,      # Standard width - adjust as needed
  height = 4,     # Standard height - adjust as needed
  dpi = 300,
  bg = "white"
)

# Save Sankey plot
htmlwidgets::saveWidget(
  sankey_plot, 
  "temp_sankey.html",
  selfcontained = TRUE
)

webshot(
  "temp_sankey.html", 
  file.path(plot_dir, "sankey_plot.png"),
  selector = ".sankeyNetwork",
  zoom = 4,
  delay = 0.2,
  vwidth = 1000,
  vheight = 1000
)

# Save polar plot
png(
  filename = file.path(plot_dir, "polar_plot.png"),
  width = 1000,
  height = 1000,
  res = 300
)
create_polar_plot()
dev.off()

# Clean up temporary files
file.remove("temp_sankey.html")
