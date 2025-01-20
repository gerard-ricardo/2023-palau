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


# conversion to grobs -----------------------------------------------------
# For sankey_plot (htmlwidget)
library(htmlwidgets)
library(webshot)
webshot::install_phantomjs()

# Save the Sankey plot as HTML and convert to PNG
htmlwidgets::saveWidget(sankey_plot, "temp_sankey.html")
webshot("temp_sankey.html", "temp_sankey.png")
sankey_plot_grob <- rasterGrob(readPNG("temp_sankey.png"))
# Clean up temporary files
file.remove(c("temp_sankey.html", "temp_sankey.png"))
class(sankey_plot_grob)


## polar plot
png(filename = "temp.png", width = 800, height = 800)
create_polar_plot()
dev.off()
# Read back as a raster and convert to grob
library(png)
temp_png <- readPNG("temp.png")
polar_plot_grob <- rasterGrob(temp_png)
# Clean up temporary file
file.remove("temp.png")
class(polar_plot_grob)



# panel -------------------------------------------------------------------


panel2 <- grid.arrange(
  # Left column: 3 plots stacked vertically
  arrangeGrob(distance_map_plot, pairwise_dist_plot, sankey_plot_grob, ncol = 1),
  # Right column: 2 plots stacked vertically
  arrangeGrob(arranged_plots, polar_plot_grob, ncol = 1),
  ncol = 2,  # Two columns: left and right
  widths = c(1, 1)  # Adjust width proportions (left column twice as wide as the right)
)

#ggsave(panel2, filename = 'fig2.pdf',  path = "./plots", device = 'pdf',  width = 9, height = 9.5)  #
