# clones and genetic relatedness------------------------------------------------------------------

#SEE COLONY BESTCLONE FOR CLONE ANALYSIS (Palau 15fullrun clones).  SAME AS SMOUSE BUT MORE DEFENSIABLE (EXCEPT C5_2)
#note: further runs in COLONy hasn't reproduce palau2, which  might be because of increased filtering



# load libraries ----------------------------------------------------------

library(igraph)
library(ggraph)
library(dplyr)



# from colony -------------------------------------------------------------

text <- '1     1.000      c1_1,c1_2
2     1.000      7_30_1,7_3_2,7_05_2,7_3_1,7_10_2,7_30_2,7_05_1,7_10_1
3     1.000      c7_1,c7_2
4     1.000      c11_1
5     1.000      c16_1,c16_2
6     1.000      c20_1,c20_2
7     1.000      5_3_1,5_3_2,5_05_2,5_05_1,5_10_2
8     1.000      4_3_2,4_10_2,4_05_1,4_10_1,4_3_1,4_05_2
9     1.000      6_30_1,6_30_2,6_3_1,6_05_2,6_3_2,6_10_2,6_10_1,6_05_1
10     1.000      3_3_2,3_3_1,3_10_2,3_05_2,3_10_1,3_05_1
11     1.000      c11_2,c17_1,5_30_1,c17_2,5_30_2
12     1.000      2_30_1,2_10_1,2_30_2
13     1.000      c3_1,c3_2,c13_1,c13_2
14     1.000      c8_1,c8_2
15     1.000      c12_1,c12_2
16     1.000      1_10_1
17     1.000      3_30_2,4_30_,3_30_1,4_30_1,c10_1,c10_2
18     1.000      c9_1,c9_2
19     1.000      c18_1,c18_2
20     1.000      1_05_,c14_1,c14_2,1_30_1
21     1.000      c5_2
22     1.000      c6_1,c6_2
23     1.000      c19_1,c19_2
24     1.000      1_3_
25     1.000      2_05_1,2_3_'

df1 <- read.table(text = text, header = FALSE, sep = "", fill = TRUE, col.names = c("CloneID", "Prob", "group")) # basic table
df1$group1 <- sapply(strsplit(df1$group, ","), function(x) x) # split comma-separated groups
df1$group <- NULL # drop raw string if unnecessary
genotype_data_long <- do.call(rbind, lapply(1:nrow(df1), function(i) data.frame(genotype2 = df1$CloneID[i], Prob = df1$Prob[i], id = df1$group[[i]])))
genotype_data_long$genotype2 = paste0('X', genotype_data_long$genotype2)
genotype_data_long = genotype_data_long %>% dplyr::select(-Prob )
str(genotype_data_long)


# # Perform MLG analysis
# mlg_analysis <- mlg(data_genind_adult)
# print(paste("MLG Analysis:", mlg_analysis, "- clearly incorrect since there are def reps"))
# # Convert genind to genclone to use @mlg slot
# data_genclone_adult <- as.genclone(data_genind_adult) # convert to genclone
# data_genclone_adult@mlg <- mlg.filter(data_genclone_adult, dist=bitwise.dist, threshold=0.05) # assign clone groups
# mlg(data_genclone_adult) # now returns clone clusters
# mlg_vec <- data_genclone_adult@mlg
# mlg_table <- data.frame(Individual=indNames(data_genclone_adult),MLG=mlg_vec)
# clone_groups <- split(mlg_table$Individual,mlg_table$MLG)
# clone_groups[sapply(clone_groups,length)>0]
# #this works much better
# #pssible mismatches C11_2; c3_2, 
# #possible clones: c11 & c17 & 5_30; c14 and 1_05/1_30; 3_30, 4_30, c10; c3 & c13
# genotype_data_long <- stack(clone_groups[sapply(clone_groups, length) > 0])
# colnames(genotype_data_long) <- c("id", "group")
# group_labels <- setNames(paste0("X", seq_along(unique(genotype_data_long$group))), unique(genotype_data_long$group))
# genotype_data_long$genotype2 <- group_labels[as.character(genotype_data_long$group)]
# 

#
real_geno_df <- stack(clone_groups) # collapse list into data.frame with values and group names
real_geno_df <- real_geno_df %>% mutate(across(c(values ), ~ gsub("05", "0.7", .)))
real_geno_df <- real_geno_df %>% mutate(values = paste0("X", values))
real_geno_df.x <- real_geno_df %>% rename(moth_id = values, real_geno.x = ind)
real_geno_df.y = real_geno_df %>% rename(fath_id = values, real_geno.y = ind)



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
# Temporary fix by adding leading zeroes
heatmap_data <- adult_colonies_sort %>%
  # mutate(Individual1 = gsub("_05_", "_0.5_", Individual1), # Replace '05' with '0.5'
  #        Individual1 = gsub("_3_", "_03_", Individual1), # Replace '3' with '03'
  #        Individual2 = gsub("_05_", "_0.5_", Individual2), # Same for Individual2
  #        Individual2 = gsub("_3_", "_03_", Individual2)) %>%
  mutate(Individual1 = factor(Individual1, levels = unique(Individual1)),
         Individual2 = factor(Individual2, levels = unique(Individual2)))# %>%
  # After sorting, revert the temporary labels back to original format
  # mutate(Individual1 = gsub("_0.5_", "_05_", Individual1), # Revert '0.5' back to '05'
  #        Individual1 = gsub("_03_", "_3_", Individual1), # Revert '03' back to '3'
  #        Individual2 = gsub("_0.5_", "_05_", Individual2), # Same for Individual2
  #        Individual2 = gsub("_03_", "_3_", Individual2)) # Same for Individual2
ggplot(heatmap_data, aes(x = Individual1, y = Individual2, fill = Distance)) +
  geom_tile() + scale_fill_gradient(low = "blue", high = "red") + # Blue = close relatives, Red = distant
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Relatedness Heatmap", x = "Individual 1", y = "Individual 2", fill = "Genetic Distance")

#reference to single colony
threshold = 200  #find this manually based on genetic distance of a few individuals
unique(adult_colonies_sort$Individual1)
first_group = 'c16_1'
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


##old grouping useing 200 threshold
# # Convert list into a dataframe where each row represents a group  (
# group_list <- lapply(result$close_rel, as.character) # Convert factors to characters
# edges <- do.call(rbind, lapply(group_list, function(g) if (length(g) > 1) combn(g, 2, simplify = FALSE))) %>%
#   do.call(rbind, .) %>% as.data.frame(stringsAsFactors = FALSE) %>% setNames(c("Individual1", "Individual2"))
# graph <- graph_from_data_frame(edges, directed = FALSE)
# clusters <- components(graph)
# grouped_individuals <- data.frame(id = names(clusters$membership), Cluster = clusters$membership)
# genotype_data <- grouped_individuals %>% group_by(Cluster) %>% summarise(id = paste(unique(id), collapse = ", ")) %>%
#   data.frame() %>% mutate(genotype2 = paste0('x', Cluster)) %>% dplyr::select(-Cluster)
# nrow(genotype_data)
# genotype_data_long <- genotype_data %>%
#   separate_rows(id, sep = ",")  %>% data.frame()
# genotype_data_long$id <- trimws(genotype_data_long$id)  #trims leading white space from labels

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
