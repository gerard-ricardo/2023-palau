
#clones
mlg_analysis <- mlg(data_genind_adult, quiet = FALSE)  #not very good this function



#Compare relatedness on single individual vs all
data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) #combine all
genetic_dist_matrix <- gd.smouse(data_genind_adult, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
genetic_dist_df <- as.data.frame(as.matrix(genetic_dist_matrix))
genetic_dist_df <- tibble::rownames_to_column(genetic_dist_df, "Individual1")
adult_colonies <- pivot_longer(genetic_dist_df, cols = -Individual1, names_to = "Individual2", values_to = "Distance") %>% data.frame()
adult_colonies_sort <- adult_colonies %>% arrange(Individual1, Distance)
hist(adult_colonies_sort$Distance, main = "Genetic Distance Distribution", xlab = "Genetic Distance", ylab = "Frequency")

#
unique(adult_colonies_sort$Individual1)
first_group = '7_5_1'
first_group_data <- adult_colonies_sort %>%
  filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)]))
#indivudal by genetic distance
p1 <- ggplot(first_group_data, aes(x = Individual2, y = Distance)) +
  geom_point() +
  facet_wrap(~ Individual1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Distances for Adult Colonies", x = "Individual2", y = "Genetic Distance")
p1  #some error between reps
