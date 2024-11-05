#PCA


## notes
# big ddiference betwen kmena sclustering and dbscan
#also difference between regualr and dDocent clustering (bugger)


library(dbscan)


# adult only ---------------------------------------------------------------------

## quick plot
# pca = gl.pcoa(data_gl_adult_unique)
# gl.pcoa.plot(glPca = pca, data_gl_adult_unique)

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
pca_complete$kmeans_cluster <- as.factor(kmeans_result$cluster)
#cluster 2 is  strong, others weak

# DBSCAN clustering
# Find the appropriate eps value using kNNdistplot
kNNdistplot(pca_data, k = 5)  #k-nearest neighbour
elbow = 9.82 # Place this at the elbow of the line
abline(h = elbow, col = "red", lty = 2)  

# # Function to perform DBSCAN clustering and plot results
# perform_dbscan <- function(pca_data, pca_complete, eps_value, min_pts = 3) {
#   dbscan_result <- dbscan(pca_data, eps = eps_value, minPts = min_pts)
#   cluster_col_name <- paste0("Cluster_dbscan_", eps_value)
#   pca_complete[[cluster_col_name]] <- as.factor(dbscan_result$cluster)
#   plot <- ggplot(pca_complete, aes_string(x = "Axis1", y = "Axis2", color = cluster_col_name)) +
#     geom_point(alpha = 0.6) +
#     labs(title = paste("PCA Plot with DBSCAN Clusters (eps =", eps_value, ")"),
#          x = "Principal Component 1",
#          y = "Principal Component 2") +
#     theme_minimal()
#   silhouette_score <- silhouette(dbscan_result$cluster, dist(pca_data))
#   print(dbscan_result)
#   print(summary(silhouette_score))
#   return(plot)
# }
# 
# eps_values <- elbow 
# for (eps in eps_values) {
#   plot <- perform_dbscan(pca_data, pca_complete, eps)
#   print(plot)
# }
#1 cluster, 3 noise

# Single DBSCAN clustering and plotting
dbscan_result <- dbscan(pca_data, eps = elbow, minPts = 3)
pca_complete$Cluster_dbscan <- as.factor(dbscan_result$cluster)

# Plotting
t2_dbscan <- ggplot(pca_complete, aes(x = Axis1, y = Axis2, color = Cluster_dbscan)) +
  geom_point(alpha = 0.6) +
  labs(title = paste("PCA Plot with DBSCAN Clusters (eps =", elbow, ")"),
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

# Silhouette score
silhouette_score <- silhouette(dbscan_result$cluster, dist(pca_data))
(dbscan_result)
(summary(silhouette_score))
(t2_dbscan)


# prepare for plotting ----------------------------------------------------

pca_complete2 <- pca_complete %>%
  mutate(
    MumID = str_extract(row.names(pca_complete), "^[^_]+"),
    NewID = paste0('Adu', "_", MumID)
  )
data1 <- dplyr::arrange(pca_complete2, Axis1) # 
pca_complete2 <- pca_complete2 %>% mutate(across(c(MumID, NewID), as.factor))
str(pca_complete2)
my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  'tomato', 'pink', 'red'
)

#
unique_pops <- length(unique(pca_complete2$pop))
if (unique_pops > length(my_palette)) {
  my_palette <- scales::hue_pal()(unique_pops)
}

#plot by genotype
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = factor(pop), color = factor(pop )), shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = pop), color = "grey50", size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = scales::alpha(my_palette, 0.1)) + # Adjust the alpha for the fill colours
  scale_color_manual(values = my_palette) +
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

## kmean clustering plot
# Get the number of unique clusters in kmeans_cluster
my_palette <- c("red", "blue")

t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = factor(kmeans_cluster), color = factor(kmeans_cluster)), shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = pop), color = "grey50", size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = scales::alpha(my_palette, 0.1)) + # Adjust alpha for fill colours
  scale_color_manual(values = my_palette) + # Apply updated palette
  stat_ellipse(aes(x = Axis1, y = Axis2, group = kmeans_cluster, color = factor(kmeans_cluster)), level = 0.95, linetype = 2, size = 1) +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "K-means Cluster", fill = "K-means Cluster", shape = "Stage"
  )
t2


##dbscan cluster
# Ensure 'Cluster_dbscan_27.4' (replace with actual column name) exists in pca_complete2
# Example: using 'Cluster_dbscan_27.4' as DBSCAN cluster column name

dbscan_cluster_col <- "Cluster_dbscan"  # Replace with your actual DBSCAN cluster column name

my_palette <- c("red", "blue")


# Plot based on DBSCAN clusters
t2_dbscan <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = factor(!!sym(dbscan_cluster_col)), color = factor(!!sym(dbscan_cluster_col))), 
             shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = pop), color = "grey50", size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = scales::alpha(my_palette, 0.1)) + # Adjust alpha for fill colours
  scale_color_manual(values = my_palette) + # Apply updated palette for DBSCAN clusters
  stat_ellipse(aes(x = Axis1, y = Axis2, group = !!sym(dbscan_cluster_col), color = factor(!!sym(dbscan_cluster_col))), 
               level = 0.95, linetype = 2, size = 1) +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "DBSCAN Cluster", fill = "DBSCAN Cluster", shape = "Stage"
  )
t2_dbscan

###################################################################################################################
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
