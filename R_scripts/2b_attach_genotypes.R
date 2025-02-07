
## Adds genotype numbers to meta data

# Notes for standard filtering --------------------------------------------
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

# genetic distance groups
#c11_1, c11_2, c17_1, c17_2, 5_30_1,5_30_2


# Notes for dDocent filtering ---------------------------------------------




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


# attach clone id (genotype to data) --------------------------------------

# genotype_data <- data.frame(
#   genotype2 = paste0("x", c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 
#                             12, 13, 14, 15, 16, 17, 19, 20, 21, 23, 
#                             24, 25, 26)),  # Add 'x' prefix
#   id = c("c1_1,c1_2", "7_30_1,7_30_2,7_3_1,7_3_2,7_5_1,7_5_2,7_10_1,7_10_2", 
#          "c7_1,c7_2", "c16_1,c16_2", "c20_1,c20_2", 
#          "5_3_1,5_3_2,5_5_1,5_5_2,5_10_2", 
#          "4_3_1,4_3_2,4_10_1,4_10_2,4_5_1,4_5_2", 
#          "6_30_1,6_30_2,6_3_1,6_3_2,6_5_1,6_5_2,6_10_1,6_10_2", 
#          "3_3_1,3_3_2,3_10_1,3_10_2,3_5_2,3_5_1", 
#          "c11_2,c17_1,c17_2,5_30_1,5_30_2,c11_1", 
#          "2_30_1,2_30_2,2_10_1", 
#          "c3_1,c3_2,c13_1,c13_2,c5_2", 
#          "c8_1,c8_2", "c12_1,c12_2", "1_10_1", 
#          "3_30_1,3_30_2,4_30_1,4_30_,c10_1,c10_2", 
#          "c9_1,c9_2", "c18_1,c18_2", 
#          "1_5_,c14_1,c14_2,1_30_1", 
#          "c6_1,c6_2", "c19_1,c19_2", 
#          "1_3_", "2_5_1,2_3_")  # Member1, Member2, etc.
# )

# Step 2: Convert to long format, splitting each 'id' entry by commas
genotype_data_long <- genotype_data %>%
  separate_rows(id, sep = ",")  %>% data.frame()
genotype_data_long$id <- trimws(genotype_data_long$id)  #trims leading white space from labels


# Add the genotype2 information to the ind.metrics in data_gl_filtered_adult
# Convert 'id' in both dataframes to character type
data_gl_filtered_adult@other$ind.metrics$id <- as.character(data_gl_filtered_adult@other$ind.metrics$id)
genotype_data_long$id <- as.character(genotype_data_long$id)

data_gl_filtered_adult@other$ind.metrics <- data_gl_filtered_adult@other$ind.metrics %>% left_join(genotype_data_long, 
                                                                                                   by = "id")  # Join by the genotype2 column
data_gl_filtered_adult@other$ind.metrics$genotype2
tt = data_gl_filtered_adult@other$ind.metrics

# # replace id name with genotype2 label 
# if (length(data_gl_filtered_adult@other$ind.metrics$genotype2) == length(data_gl_filtered_adult@ind.names)) {
#   # Replace ind.names with genotype2
#   data_gl_filtered_adult@ind.names <- data_gl_filtered_adult@other$ind.metrics$genotype2
# } else {
#   warning("Length of genotype2 does not match length of ind.names")
# }


# population filtering and objects ----------------------------------------

# look into genotype as population
data_gl_filtered_adult@pop <- as.factor(data_gl_filtered_adult@other$ind.metrics$genotype2) # Convert genotype2 to factor and assign
data_gl_filtered_adult
#81 gentopyes

## temporarily remove missing genotypes
# Filter out individuals with missing genotype2
valid_indices <- !is.na(data_gl_filtered_adult@other$ind.metrics$genotype2)
data_gl_filtered_adult <- data_gl_filtered_adult[valid_indices, ] # Keep only valid rows
data_gl_filtered_adult@pop <- as.factor(data_gl_filtered_adult@other$ind.metrics$genotype2[valid_indices]) # Assign genotype2 as factor to pop



# Convert GENIND OBJECT all indiv
data_genind <- gl2gi(data_gl_filtered)

# Convert genind adults only
data_genind_adult <- gl2gi(data_gl_filtered_adult)


#select the best rep (unique adults) based on best callrate
genotype_matrix <- data_genind_adult@tab
(callrate <- rowMeans(!is.na(genotype_matrix)))
ind_names <- indNames(data_genind_adult)
genotypes <- data_genind_adult@other$ind.metrics$genotype2  # Adjust if necessary
geno_df <- data.frame(individual = ind_names, genotype = genotypes, callrate = callrate, stringsAsFactors = FALSE)
best_geno_df <- geno_df %>% group_by(genotype) %>% slice_max(order_by = callrate, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% data.frame()
best_ind_names <- best_geno_df$individual
best_indices <- match(best_ind_names, indNames(data_genind_adult))
data_genind_adult_unique <- data_genind_adult[best_indices, ]
data_gl_adult_unique = gi2gl(data_genind_adult_unique, parallel = FALSE, verbose = NULL)
data_gl_adult_unique@pop


