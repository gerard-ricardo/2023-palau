# clones and genetic relatedness------------------------------------------------------------------

#SEE COLONY BESTCLONE FOR CLONE ANALYSIS (Palau2).  SAME AS SMOUSE BUT MORE DEFENSIABLE (EXCEPT C5_2)




library(igraph)
library(ggraph)

# Perform MLG analysis
mlg_analysis <- mlg(data_genind_adult)
print(paste("MLG Analysis:", mlg_analysis, "- clearly incorrect since there are def reps"))


# Compare relatedness on single individual vs all
data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) # combine all
genetic_dist_matrix <- gd.smouse(data_genind_adult, verbose = TRUE) # Smouse and Peakall (1999) is a method used to quantify the
genetic_dist_df <- as.data.frame(as.matrix(genetic_dist_matrix))
genetic_dist_df <- tibble::rownames_to_column(genetic_dist_df, "Individual1")
adult_colonies <- pivot_longer(genetic_dist_df, cols = -Individual1, names_to = "Individual2", values_to = "Distance") %>% data.frame()
adult_colonies_sort <- adult_colonies %>% arrange(Individual1, Individual2)
str(adult_colonies_sort)
plot(density(adult_colonies_sort$Distance, main = "Genetic Distance Distribution", xlab = "Genetic Distance",
             ylab = "Frequency"))
# heatmap
heatmap_data <- adult_colonies_sort %>%
  mutate(Individual1 = factor(Individual1, levels = unique(Individual1)),
         Individual2 = factor(Individual2, levels = unique(Individual2)))
ggplot(heatmap_data, aes(x = Individual1, y = Individual2, fill = Distance)) +
  geom_tile() + scale_fill_gradient(low = "blue", high = "red") + # Blue = close relatives, Red = distant
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Relatedness Heatmap", x = "Individual 1", y = "Individual 2", fill = "Genetic Distance")

#reference to single colony
threshold = 200  #find this manually based on genetic distance of a few individuals
unique(adult_colonies_sort$Individual1)
first_group = '4_30_'
first_group_data <- adult_colonies_sort %>% filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) 
#relatedness per indviduals
p1 <- ggplot(first_group_data, aes(x = Individual2, y = Distance)) +
  geom_point() +
  facet_wrap(~ Individual1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Distances for Adult Colonies", x = "Individual2", y = "Genetic Distance")
p1  #some error between reps. Typically around XXX for clones
#their groups
(first_group_data_names <- adult_colonies_sort %>%
    filter(Individual1 == first_group) %>%
    mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) %>% 
    filter(Distance < threshold)  %>% reframe(Individual2) %>% as.vector())


## batch run for each Indiv1
# Get all unique Individual1 values
individual_list <- unique(adult_colonies_sort$Individual1)
# Iterate through each Individual1 and extract the Individual2 where Distance < XXX
result <- map_df(individual_list, function(first_group) {
  first_group_data <- adult_colonies_sort %>%
    filter(Individual1 == first_group) %>%
    mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) %>%
    filter(Distance < threshold)  # Filter for Distance < XXX
  # Store the results for each Individual1
  tibble(
    Individual1 = first_group,
    close_rel = list(first_group_data$Individual2)  # Store as list
  )
})

result$close_rel[result$Individual1 == first_group]
result$close_rel


# Convert list into a dataframe where each row represents a group
group_list <- lapply(result$close_rel, as.character) # Convert factors to characters
edges <- do.call(rbind, lapply(group_list, function(g) if (length(g) > 1) combn(g, 2, simplify = FALSE))) %>%
  do.call(rbind, .) %>% as.data.frame(stringsAsFactors = FALSE) %>% setNames(c("Individual1", "Individual2"))
graph <- graph_from_data_frame(edges, directed = FALSE)
clusters <- components(graph)
grouped_individuals <- data.frame(id = names(clusters$membership), Cluster = clusters$membership)
genotype_data <- grouped_individuals %>% group_by(Cluster) %>% summarise(id = paste(unique(id), collapse = ", ")) %>%
  data.frame() %>% mutate(genotype2 = paste0('x', Cluster)) %>% dplyr::select(-Cluster)
nrow(genotype_data)

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
  theme_void()  # 
p0

#how many groups
unique_genotypes <- V(graph)$name  # extract unique genotype names
clusters <- components(graph)  # get connected components
group_list <- split(V(graph)$name, clusters$membership)  # group genotypes by component
length(group_list)
#remove reps
#add a '?' to fill short labels
group_list <- lapply(group_list, function(x) ifelse(grepl("_$", x), paste0(x, "?"), x)) # Add '?' to strings ending with '_'
group_list1 <- lapply(group_list, function(x) substr(x, 1, nchar(x) - 2)) # Remove last two characters from each element
group_list2 <- lapply(group_list1, function(x) {
  unique_elements <- unique(x)  # Remove duplicates within the vector
  paste(unique_elements, collapse = " ")  # Recombine into a single string
})  #remove duplicates
# Convert each list element to a proper character vector
group_list2 <- lapply(group_list2, function(x) unlist(strsplit(x, " "))) # Split string into vector



#  checking high hetero for evidence of contamination
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



# #Compare relatedness on single individual vs all
# data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) #combine all
# genetic_dist_matrix <- gd.smouse(data_genind_adult, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
# genetic_dist_df <- as.data.frame(as.matrix(genetic_dist_matrix))
# genetic_dist_df <- tibble::rownames_to_column(genetic_dist_df, "Individual1")
# adult_colonies <- pivot_longer(genetic_dist_df, cols = -Individual1, names_to = "Individual2", values_to = "Distance") %>% data.frame()
# adult_colonies_sort <- adult_colonies %>% arrange(Individual1, Distance)
# hist(adult_colonies_sort$Distance, main = "Genetic Distance Distribution", xlab = "Genetic Distance", ylab = "Frequency")
# 
# #
# unique(adult_colonies_sort$Individual1)
# first_group = '7_5_1'
# first_group_data <- adult_colonies_sort %>%
#   filter(Individual1 == first_group) %>%
#   mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)]))
# #indivudal by genetic distance
# p1 <- ggplot(first_group_data, aes(x = Individual2, y = Distance)) +
#   geom_point() +
#   facet_wrap(~ Individual1, scales = "free_x") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "Genetic Distances for Adult Colonies", x = "Individual2", y = "Genetic Distance")
# p1  #some error between reps
