
##notes
# use batch to workout clones/fragments
# COLONY output (palau 2) for clones aswell to confirm (did all ther filtering remove indiv?)
#    - c11_1, c11_2 mismatch  (fine to join group 11, COLONY false negative). C11 may be somewhat continatmed (lots of hetero)
#    -5_10_2, 5_10_1 mismatch   - likely mislabelled 3_5_1 with 5_10_1
#    - 3_5_2 (no rep)
#    - 2_10_1 (no rep)
#    - 1_10_1 (no rep)
#   -  c5_1, c5_2  mismatch - does seem weird, not clearly mislabelled but much worse than others. Possible contamination?
#    - 1_5_  (no rep)
#    - 1_3_  (no rep)
#   - 2_5_1  (no rep)
#   - 2_3_  (no rep)
#    - 1_30_1 (no rep)

# after corrections, 23 groups (assuming c5_1 is actually a rep

# CloneIndex      Prob      Member1,Member2,Member3,Member4,Member5,Member6,Member7,Member8 
# 1     1.000      c1_1,c1_2  (good)
# 2     1.000      7_30_1,7_30_2, 7_3_1, 7_3_2, 7_5_1, 7_5_2, 7_10_1, 7_10_2 (good)
# 3     1.000      c7_1,c7_2 (good)

# 5     1.000      c16_1,c16_2 (good)
# 6     1.000      c20_1,c20_2 (good)
# 7     1.000      5_3_1,5_3_2,5_5_1, 5_5_2, 5_10_2 (missing rep 5_10_2)
# 8     1.000      4_3_1, 4_3_2, 4_10_1, 4_10_2, 4_5_1, 4_5_2 (good)
# 9     1.000      6_30_1,6_30_2,6_3_1,6_3_2, 6_5_1, 6_5_2,6_10_1,6_10_2 (good)
# 10     1.000      3_3_1, 3_3_2, 3_10_1, 3_10_2, 3_5_2,3_5_1
# 11     1.000      c11_2, c17_1,c17_2, 5_30_1,5_30_2, c11_1 
# 12     1.000      2_30_1,2_30_2, 2_10_1 (good)
# 13     1.000      c3_1,c3_2,c13_1,c13_2, c5_2 (good)
# 14     1.000      c8_1,c8_2  (good)
# 15     1.000      c12_1,c12_2  (good)
# 16     1.000      1_10_1 (missing reps)
# 17     1.000      3_30_1, 3_30_2, 4_30_1, 4_30_, c10_1, c10_2 (good)

# 19     1.000      c9_1,c9_2 (good)
# 20     1.000      c18_1,c18_2  (good)
# 21     1.000      1_5_, c14_1, c14_2, 1_30_1 (missing reps)

# 23     1.000      c6_1,c6_2 (good)
# 24     1.000      c19_1,c19_2 (good)
# 25     1.000      1_3_ (missing reps)
# 26     1.000      2_5_1,2_3_ (missing reps)



# genetic distance grousps
#c11_1, c11_2, c17_1, c17_2, 5_30_1,5_30_2



# load libraries ----------------------------------------------------------
library(poppr)
library(magrittr)
library(adegenet) 
library(pegas)  
library(dplyr)
library(purrr)
library(tidyverse)
library(igraph)
library(ggraph)

# clones and genetic relatedness------------------------------------------------------------------

# Perform MLG analysis
mlg_analysis <- mlg(data_genind_adult)
#82 distinct individuals - doesnt work again

mlg_analysis <- mll(data_genind_adult, threshold = 0.02) # Adjust threshold value as needed


#Compare relatedness on single individual vs all
data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) #combine all
genetic_dist_matrix <- gd.smouse(data_genind_adult, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
genetic_dist_df <- as.data.frame(as.matrix(genetic_dist_matrix))
genetic_dist_df <- tibble::rownames_to_column(genetic_dist_df, "Individual1")
adult_colonies <- pivot_longer(genetic_dist_df, cols = -Individual1, names_to = "Individual2", values_to = "Distance") %>% data.frame()
adult_colonies_sort <- adult_colonies %>% arrange(Individual1, Distance)
plot(density(adult_colonies_sort$Distance, main = "Genetic Distance Distribution", xlab = "Genetic Distance", ylab = "Frequency"))

#
unique(adult_colonies_sort$Individual1)
first_group = 'c5_1'
first_group_data <- adult_colonies_sort %>%
  filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) 

(first_group_data_names <- adult_colonies_sort %>%
  filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) %>% 
  filter(Distance < 2000)  %>% reframe(Individual2) %>% as.vector())
#individual by genetic distance
p1 <- ggplot(first_group_data, aes(x = Individual2, y = Distance)) +
  geom_point() +
  facet_wrap(~ Individual1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Distances for Adult Colonies", x = "Individual2", y = "Genetic Distance")
p1  #some error between reps. Typically around 500 for clones

## batch run for each Indiv1
# Get all unique Individual1 values
individual_list <- unique(adult_colonies_sort$Individual1)
# Iterate through each Individual1 and extract the Individual2 where Distance < 2000
result <- map_df(individual_list, function(first_group) {
  first_group_data <- adult_colonies_sort %>%
    filter(Individual1 == first_group) %>%
    mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) %>%
    filter(Distance < 3000)  # Filter for Distance < 2000
  # Store the results for each Individual1
  tibble(
    Individual1 = first_group,
    close_rel = list(first_group_data$Individual2)  # Store as list
  )
})

result$close_rel[result$Individual1 == first_group]
result$close_rel

# visualise groups 
# Expand the list column into individual rows
edges <- result %>%
  unnest(close_rel) %>%
  filter(!is.na(close_rel))  # Remove rows with NA values in close_rel
graph <- graph_from_data_frame(edges, directed = FALSE)
p0 <- ggraph(graph, layout = "fr") +  # fr = Fruchterman-Reingold layout for better separation
  geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) +  # Draw edges
  geom_node_point(color = "dodgerblue", size = 5) +  # Draw nodes
  geom_node_text(aes(label = name), vjust = 1.5, hjust = 0.5, size = 3) +  # Add labels to nodes
  theme_void()  # Clean up the background for a minimalist look
p0

#  checking hetero for evidence of contamination
Hobs <- function(x) {
  apply(tab(x), 1, function(ind) {
    # calculate the number of loci that are heterozygous 
    heterozygous_loci <- sum(ind == 1, na.rm = TRUE)
    non_missing_loci <- sum(!is.na(ind))
    heterozygous_loci / non_missing_loci
  })
}
hetero <- Hobs(data_genind_adult)
filt_hetero <- hetero[hetero > 0.17]
barplot(hetero, las = 2) # Adjust the colour as needed


# Calculate Nei's distance for additional insights (doesnt work)
# nei_dist <- nei.dist(data_genind,warning = TRUE)
# nei_dist_df <- as.data.frame(as.matrix(nei_dist))
# print(nei_dist_df)
# 
# # Calculate Roger's distance for additional insights (deosnt work)
# rogers_dist <- rogers.dist(data_genind)
# rogers_dist_df <- as.data.frame(as.matrix(rogers_dist))
# print(rogers_dist_df)
