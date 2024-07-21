# 2023 Acropora hyacinthus ------------------------------------------------

# load libraries ----------------------------------------------------------
#install.packages('dartr')
library(dartR)
library(dartR.popgen)
library(PopGenReport)
library(adegenet)
library(tictoc)
library(HardyWeinberg)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggrepel)
library(hierfstat)
library(ape)
library(poppr)
library(pegas)
library(dbscan)
library(sp)
library(rgdal)
library(clustertend)
library(cluster)
library(plotly)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code

# import data -------------------------------------------------------------


#meta_2023_acro <- read.csv("./data/Report_DAc24-9371_SNP_2_copy.csv", head = T, skip = 6) # make sure samples are in same order as in data_gl
# data_gl <- gl.read.dart(filename = "./data/Report_DAc24-9371_SNP_2_copy.csv", ind.metafile = "./data/2023_palau_meta.csv", topskip = 6)
# 
# # recalculate metrics
# data_gl <- gl.recalc.metrics(data_gl, v = 3) # recalculate loci metrics
# 
# save(data_gl, file = file.path("./Rdata", "2023_Acro_hyac_gl.RData"))
load("./Rdata/2023_Acro_hyac_gl.RData")  #data_gl


data_gl$other$loc.metrics
data_gl$other$loc.metrics$coverage <- data_gl$other$loc.metrics$AvgCountRef + data_gl$other$loc.metrics$AvgCountSnp
mean(data_gl$other$loc.metrics$coverage) # AH = 35.8
min((data_gl$other$loc.metrics$coverage)) #  AH = 5
max((data_gl$other$loc.metrics$coverage)) #  AH = 372
sd(data_gl$other$loc.metrics$coverage) / sqrt(1996) #AH = 0.66

# data filtering ----------------------------------------------------------
data_gl_filtered <- data_gl
#following the dartr suggested order (see tut5)
# gl <-gl.filter.secondaries(gl)
# gl <- gl.filter.rdepth(gl)
# gl <- gl.filter.reproducibility(gl)
# gl <- gl.filter.callrate(gl, method=”loc”)
# gl <- gl.filter.callrate(gl, method=”ind”)
# gl <- gl.filter.monomorphs(gl)

## pre-filtering
#indiv = 218  ind, 50405  loc

##secondaries
gl.report.secondaries(data_gl_filtered)
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method="random", verbose = 3) #remove loci fragment that shared SNPs. Only keep 1
#indiv = 218  ind, 34552   loc

#rdepth
gl.report.rdepth(data_gl_filtered)
#platy = 13.8 Generally 10 is considered min
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered,  lower = 10, v = 3) # filter by loci callrate
#indiv = 218 ind, 21260   loc, med depth = 23

##reproducibility 
gl.report.reproducibility(data_gl_filtered )
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.95, v=3) #filter out loci with limited reproducibility
# at 95%: ind = 218, loci = 20276   

# callrate loci (non missing data)
gl.report.callrate(data_gl_filtered, method = "loc") 
# AH = 0.59 (41% missing data)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
#at 70% AH = 218  ind, 7993  loci

#Minor Allele Frequency (MAF) and Coverage Filter:
list.match <- data_gl_filtered$loc.names[
  which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.01 & 
          data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef > 0.01 & 
          data_gl_filtered$other$loc.metrics$coverage > 4)
]
data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)]
#at 70% AH = 218  ind, 5678   loci

## call rate ind (non missing data). low could indicate poor extract or reference genome or contamination.
#individauls
gl.report.callrate(data_gl_filtered, method = "ind") 
# platy = 87%
pre_filt_ind <- data_gl_filtered@ind.names
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.7, v = 3) # filter by ind callrate
filt_ind <- data_gl_filtered@ind.names
(lost_ind <- setdiff(pre_filt_ind, filt_ind))
# at 70% AH = 204  ind, 5678   loci
#just a few larvae lost

data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3) # recalculate loci metrics


##others - not sure if needed
# data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
# not sure if I need HWE filter because remove inbreeding
# data_gl_filtered <- gl.filter.hwe(data_gl_filtered, alpha_val = 0.05, subset = "each", multi_comp_method = 'bonferroni',v=3) #filter out loci that depart from H-W proportions
# list.match <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.05 & data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef > 0.05 & data_gl_filtered$other$loc.metrics$coverage > 5)] #remove loci based on minor allele frequency and low data coverage
# data_gl_filtered <- data_gl_filtered[,match(list.match, data_gl_filtered$loc.names)]#keep only loci in the list above

# population filtering and objects ----------------------------------------


# look into genotype as population
data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "id")
data_gl_filtered


# Convert GENIND OBJECT all indiv
data_genind <- gl2gi(data_gl_filtered)

# Filter out eggs and larvae to keep only adults
adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)

# Convert genind adults only
data_genind_adult <- gl2gi(data_gl_filtered_adult)





# adult only ---------------------------------------------------------------------

#quick plot
pca = gl.pcoa(data_gl_filtered_adult)
gl.pcoa.plot(glPca = pca, data_gl_filtered_adult)

# PCA Analysis
pca_data <- tab(data_gl_filtered_adult, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete <- data.frame(pca$li, pop = data_gl_filtered_adult$pop) # Combine PCA results with population data
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
elbow = 8.7 # Place this at the elbow of the line
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
#The clustering contains 2 cluster(s) and 2 noise points (id = 12).


# plotting
pca_complete <- pca_complete %>%
  mutate(
    Stage = ifelse(str_detect(row.names(pca_complete), "_"), "Adult", "Larva"),
    MumID = str_sub(row.names(pca_complete), 1, -2),
    RepID = str_sub(row.names(pca_complete), -1, -1),
    NewID = paste0(Stage, "_", MumID, "_", RepID)
  )

data1 <- dplyr::arrange(pca_complete, Axis1) # 
pca_complete <- pca_complete %>% mutate(across(c(Stage, MumID, RepID, NewID), as.factor))
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

#color individuals
t2 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = factor(pop)), shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = pop ), size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
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
t2_interactive <- ggplotly(t2)
t2_interactive

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
    NewID = paste0(MumID, "_", genotype)
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
head(pca_complete2)
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = genotype, shape = stage, color = ifelse(grepl("larvae", stage), "red", "black")),
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = NewID), size = 3, max.overlaps = 100, point.padding = 0.5, box.padding = 0.5) +
  #stat_ellipse(aes(x = Axis1, y = Axis2, group = Cluster, color = Cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen", "red" = "red", "black" = "black")) +
  scale_shape_manual(values = c("adults" = 22, "larvae" = 21)) + # Set shapes: squares for adults and circles for larvae
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Cluster", fill = "Population", shape = "Stage") # Add labels to the axes and legend
t2

t2_interactive <- ggplotly(t2)
t2_interactive
