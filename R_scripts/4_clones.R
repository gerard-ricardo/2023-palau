
##notes
# use batch to workout clones/fragments





# load libraries ----------------------------------------------------------
library(poppr)
library(magrittr)
library(adegenet) 
library(pegas)  
library(dplyr)
library(purrr)

# clones and genetic relatedness------------------------------------------------------------------

# Perform MLG analysis
mlg_analysis <- mlg(data_genind_adult, quiet = FALSE)
#82 distinct individuals - doesnt work again

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
first_group = 'c20_2'
first_group_data <- adult_colonies_sort %>%
  filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)]))
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
    filter(Distance < 1000)  # Filter for Distance < 2000
  # Store the results for each Individual1
  tibble(
    Individual1 = first_group,
    close_rel = list(first_group_data$Individual2)  # Store as list
  )
})

result$close_rel[result$Individual1 == first_group]




# Calculate Nei's distance for additional insights (doesnt work)
nei_dist <- nei.dist(data_genind,warning = TRUE)
nei_dist_df <- as.data.frame(as.matrix(nei_dist))
print(nei_dist_df)

# Calculate Roger's distance for additional insights (deosnt work)
rogers_dist <- rogers.dist(data_genind)
rogers_dist_df <- as.data.frame(as.matrix(rogers_dist))
print(rogers_dist_df)