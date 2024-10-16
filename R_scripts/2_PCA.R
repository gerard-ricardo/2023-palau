#PCA



# adult only ---------------------------------------------------------------------

#quick plot
pca = gl.pcoa(data_gl_adult_unique)
gl.pcoa.plot(glPca = pca, data_gl_adult_unique)

# PCA Analysis
pca_data <- tab(data_gl_adult_unique, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete <- data.frame(pca$li, pop = data_gl_adult_unique$pop) # Combine PCA results with population data
#use for adults
#pca_complete <- data.frame(pca$li) # Combine PCA results with population data

# Explained variance
(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()

# Hopkins statistic
set.seed(123) # for reproducibility
#(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))
# Calculated values 0-0.3 indicate regularly-spaced data. Values around 0.5 indicate random data. Values 0.7-1 indicate clustered data.
#PD = 0.14


# K-means clustering
set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 2, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
plot(silhouette_score)
pca_complete$Cluster <- as.factor(kmeans_result$cluster)
#PD: cluster 3 is quite strong, others poor to mod. 

# DBSCAN clustering
# Find the appropriate eps value using kNNdistplot
kNNdistplot(pca_data, k = 5)  #k-nearest neighbour
elbow = 29.3 # Place this at the elbow of the line
abline(h = elbow, col = "red", lty = 2)  
library(dbscan)
# Function to perform DBSCAN clustering and plot results
perform_dbscan <- function(pca_data, pca_complete, eps_value, min_pts = 5) {
  dbscan_result <- dbscan(pca_data, eps = eps_value, minPts = min_pts)
  cluster_col_name <- paste0("Cluster_dbscan_", eps_value)
  pca_complete[[cluster_col_name]] <- as.factor(dbscan_result$cluster)
  plot <- ggplot(pca_complete, aes_string(x = "Axis1", y = "Axis2", color = cluster_col_name)) +
    geom_point(alpha = 0.6) +
    labs(title = paste("PCA Plot with DBSCAN Clusters (eps =", eps_value, ")"),
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal()
  silhouette_score <- silhouette(dbscan_result$cluster, dist(pca_data))
  print(dbscan_result)
  print(summary(silhouette_score))
  return(plot)
}

eps_values <- elbow 
for (eps in eps_values) {
  plot <- perform_dbscan(pca_data, pca_complete, eps)
  print(plot)
}
#no clusters. not sure why noise not picked up


# plotting
pca_complete <- pca_complete %>%
  mutate(
    #Stage = ifelse(str_detect(row.names(pca_complete), "_"), "Adult", "Larva"),
    MumID = str_extract(row.names(pca_complete), "^[^_]+"),
    #RepID = str_sub(row.names(pca_complete), -1, -1),
    NewID = paste0('Adu', "_", MumID)
  )

data1 <- dplyr::arrange(pca_complete, Axis1) # 
pca_complete <- pca_complete %>% mutate(across(c(MumID, NewID), as.factor))
str(pca_complete)
my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato"
)

# Extend the palette if necessary
unique_pops <- length(unique(pca_complete$pop))
if (unique_pops > length(my_palette)) {
  my_palette <- scales::hue_pal()(unique_pops)
}

t2 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = factor(pop), color = factor(pop)), shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = NewID), color = "grey50", size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = scales::alpha(my_palette, 0.1)) + # Adjust the alpha for the fill colours
  scale_color_manual(values = my_palette) +
  #stat_ellipse(aes(x = Axis1, y = Axis2, group = pop, color = pop), level = 0.95, linetype = 2, size = 1) +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Population", fill = "Population", shape = "Stage"
  )
t2


# Convert the ggplot to an interactive plotly plot
# t2_interactive <- ggplotly(t2)
# t2_interactive

# 
# #per cluster
# t2 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
#   geom_point(aes(color = factor(Cluster)), shape = 22, 
#              size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)
#   ) +
#   geom_text_repel(aes(label = NewID), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5) +
#   scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen")) +
#   stat_ellipse(aes(x = Axis1, y = Axis2, group = Cluster, color = Cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
#   theme_sleek2() +
#   labs(
#     x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
#     y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
#     color = "Cluster", fill = "Population", shape = "Stage"
#   ) 
# t2



head(pca_complete)

# t2 <- ggplot(pca_complete, aes(x=Axis1, y=Axis2)) +
#   #geom_point(aes(fill=pop, color=ifelse(grepl("Larva", Stage), "red", NA), shape=ifelse(grepl("Larva", Stage), 22, 21)), size=4, stroke=2) +
#   geom_point(aes(fill=Stage, color=ifelse(grepl("Larva", Stage), "red", 'black')), shape= 21, size=4, stroke=2, alpha = 0.8) +
#   geom_text(aes(label = NewID), hjust=1, vjust=1, size=4) +
#   scale_fill_manual(values=my_palette) +
#   scale_color_manual(values=c("red"= "red", "black"= "black")) +
#   scale_shape_manual(values=c("21", "22")) +
#   theme_sleek2() +
#   labs(x = "PCA1", y = "PCA2")  # Add labels to the axes
# t2


#per cluster
t2 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = factor(Cluster)), shape = 22,
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)
  ) +
  geom_text_repel(aes(label = NewID), size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen")) +
  stat_ellipse(aes(x = Axis1, y = Axis2, group = Cluster, color = Cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Cluster", fill = "Population", shape = "Stage"
  )
t2


# adult and larvae --------------------------------------------------------
#quick plot
# pca = gl.pcoa(data_gl_filtered)
# gl.pcoa.plot(glPca = pca, data_gl_filtered)

# PCA Analysis
pca_data <- tab(data_gl_filtered, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete2 <- data.frame(pca$li, pop = data_gl_filtered$pop, id = data_gl_filtered$ind.names, 
                            genotype = data_gl_filtered$other$ind.metrics$genotype,
                            stage = data_gl_filtered$other$ind.metrics$stage) # Combine PCA results with population data
#use for adults
#pca_complete2 <- data.frame(pca$li) # Combine PCA results with population data

# Explained variance
(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()

# Hopkins statistic
set.seed(123) # for reproducibility
#(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))
# Calculated values 0-0.3 indicate regularly-spaced data. Values around 0.5 indicate random data. Values 0.7-1 indicate clustered data.
#all AH = XXX

# K-means clustering
set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 3, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
plot(silhouette_score)
pca_complete2$Cluster <- as.factor(kmeans_result$cluster)
#PD: cluster 3 is quite strong, others poor to mod. 

#pca_complete2 <- pca_complete2[complete.cases(pca_complete2), ] # make sure import matches NA type
pca_complete2 <- pca_complete2 %>%
  mutate(
    MumID = str_sub(stage, 1, 3),
    NewID = paste0(MumID, genotype),
    color = ifelse(grepl("larvae", stage), "2", "1")
  )

my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato"
)

# Extend the palette if necessary
unique_pops <- length(unique(pca_complete2$genotype))
if (unique_pops > length(my_palette)) {
  my_palette <- scales::hue_pal()(unique_pops)
}

# Plot with ggrepel for label lines
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = genotype, shape = stage, color = color),
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = NewID, color = color),
                  size = 3, max.overlaps = 80, point.padding = 0.7, box.padding = 0.6) +
  #stat_ellipse(aes(x = Axis1, y = Axis2, group = Cluster, color = Cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = c("1" = "grey", "2" = "lightcoral", "3" = "mediumseagreen", "red" = "red", "black" = "black", "lightcoral", 'grey')) +
  scale_shape_manual(values = c("adults" = 22, "larvae" = 21)) + # Set shapes: squares for adults and circles for larvae
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Cluster", fill = "Population", shape = "Stage") # Add labels to the axes and legend
t2


# t2_interactive <- ggplotly(t2)
# t2_interactive
