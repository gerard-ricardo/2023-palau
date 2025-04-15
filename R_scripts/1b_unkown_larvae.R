##DdOCENT FILTERING PLUS ADDED

#added some addition following orginal filtering, Jenkins et al 2024, and Premachandra 2019

#ATTEMPT TO FIND UNKOWN LARVAE INC X1-X4 USIGN PARENTAL (NOT PATERNITY)

#some mismatch between seq and col sizes. 
# - c2 missing in 96 well palte, -ok issue is c1-c4 were artibuatry Ids for the centre patch. Dont actually know positions 
#could possibly run cervus for these without know paternity (for the moment remove, removes 17 indiduls)
#the images show clearly c4 and c16 missing, also no image of c2
# sequencing is missing c4, c15, so likely c15 should be c16 

#single adults results in almost same results

#x1 likely c19 and x2 and x2 likely c18, all crossing with c12. Unkown larvae seem unerliable (selfs)

# Load required libraries ----------------------------------------------------
library(dartR)
library(adegenet)
library(tidyverse)
library(HardyWeinberg)
library(pegas)
library(dartR.popgen)
library(PopGenReport)
library(tictoc)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)
library(hierfstat)
library(ape)
library(poppr)
library(dbscan)
library(sp)
library(rgdal)
library(clustertend)
library(cluster)
library(plotly)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code

# Import data ----------------------------------------------------------------
# Read in your DArT data and metadata
# data_gl <- gl.read.dart(filename = "./data/Report_DAc24-9371_SNP_2_copy.csv",
#                         ind.metafile = "./data/2023_palau_meta_2.csv", topskip = 6)
# 
# # Clean up sample names if necessary
# indNames(data_gl) <- gsub("5_10_1", "3_5_1", indNames(data_gl))
# data_gl <- data_gl[!indNames(data_gl) %in% "c5_1"]
# 
# # Recalculate metrics after initial cleaning
# data_gl <- gl.recalc.metrics(data_gl, v = 3)
# 
# save(data_gl, file = file.path("./Rdata", "2023_Acro_hyac_gl_dDocent.RData"))
load("./Rdata/2023_Acro_hyac_gl_dDocent.RData")  #data_gl



# Fix labels ------------------------------------------------------------

# Create a filtered object from raw data
data_gl_filtered <- data_gl
# Initial dataset: 217 individuals, 50405 loci

#replace '5 with 05 labels
data_gl_filtered@other$ind.metrics$genotype <- gsub("(?<=_)5(?=$)", "05", data_gl_filtered@other$ind.metrics$genotype,
                                                    perl = TRUE) # Replace '_5' with '_05'

data_gl_filtered@other$ind.metrics$id <- gsub("(_5)(?![0-9])", "_05", data_gl_filtered@other$ind.metrics$id, perl = TRUE) 

data_gl_filtered@ind.names <- gsub("(_5)(?![0-9])", "_05", data_gl_filtered@ind.names, perl = TRUE) # Replace '_5' with '_05' when not followed by a digit

tail(data_gl_filtered@other$ind.metrics, 10)

all_id = data_gl_filtered@other$ind.metrics$genotype %>% unique()
centreid = all_id %>% grep('c', .) 
centreid %>% length()
all_id[centreid]   #missing C2, c4, c15, 

## c16 should be c15 (there is no c16)
data_gl_filtered@other$ind.metrics$genotype <- gsub("c16", "c15", data_gl_filtered@other$ind.metrics$genotype)
data_gl_filtered@other$ind.metrics$genotype %>% unique()

# ##try single reps
# data_gl_filtered  
# # Get IDs of all larvae
# df_ind_metrics <- data_gl_filtered@other$ind.metrics
# larvae_ids <- df_ind_metrics %>%
#   filter(stage == 'larvae') %>%
#   pull(id)
# # Get IDs of one adult per distinct genotype
# adult_unique_ids <- df_ind_metrics %>%
#   filter(stage == 'adults') %>%
#   distinct(genotype, .keep_all = TRUE) %>% # Keep only the first row for each unique genotype
#   pull(id)
# keep_ids <- c(adult_unique_ids, larvae_ids)
# data_gl_filtered <- data_gl_filtered[data_gl_filtered@ind.names %in% keep_ids, ]
# df_ind_metrics2 <- data_gl_filtered@other$ind.metrics

##remove all larvae with c1-c4 (they were not 100% known)
# Get IDs of individuals to keep (i.e., not larvae from c1–c4)
data_gl_filtered  #217
keep_ids <- data_gl_filtered@other$ind.metrics %>%
  filter(
    (stage == 'larvae' & genotype %in% c('x1', 'x2', 'x3', 'x4', 'unkn1', 'unkn2', 'unkn3', 'unkn4')) |
      (stage == 'adults')
  ) %>%
  pull(id)
length(keep_ids)
#data_gl_filtered <- data_gl_filtered[na.omit(match(keep_ids, data_gl_filtered@ind.names)), ]

keep_ids <- unique(keep_ids)
data_gl_filtered <- data_gl_filtered[data_gl_filtered@ind.names %in% keep_ids, ]

tt = arrange(data_gl_filtered@other$ind.metrics, stage , genotype )

data_gl_filtered
#100  individuals, 50405  loci

# Filter ------------------------------------------------------------------

# Step 4: Remove Monomorphic Loci --------------------------------------------
data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, verbose = 2)
# After filtering: 200  individuals, 50405  loci


# Step 1: Remove Individuals with High Missing Data --------------------------
ind_before <- indNames(data_gl_filtered)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.5, v = 3)
ind_after <- indNames(data_gl_filtered)
(removed_individuals <- setdiff(ind_before, ind_after))
# After filtering: 185 individuals, 50344  loci
#"160" "66"  "61"  "154" "169" "53"  "43"  "110" "210" "108" "202" "213" "87"  "206" "78"  "111
#3 of the 8  previous filt picked up as selfs were dropped here

# Step 2: Remove Loci with High Missing Data ---------------------------------
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.5, v = 3)
# After filtering: 185 individuals, 32097  loci

# Step 3: Filter by Minor Allele Frequency (MAF) i.e rare alleles----------------------------
data_gl_filtered <- gl.filter.maf(data_gl_filtered, threshold = 0.05, v = 3)
mean(data_gl_filtered@other$loc.metrics$maf)
# After filtering: 185 individuals, 19563 loci


# Step 5: Filter Loci by Read Depth ------------------------------------------
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = 3, v = 3)
# After filtering: 185 individuals, 16504 loci

# Recalculate metrics after initial cleaning
data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3)
# After filtering: 185 individuals, 19,452 loci

# Step 6: Apply Stricter Filtering on Loci -----------------------------------
# Remove loci with call rate < 95%
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.95, v = 3)
# After filtering: 185 individuals, 3,988 loci

# Filter loci with mean read depth < 20
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = 20, v = 3)
# After filtering: 185 individuals, 2644 loci

# Step 7: Remove Loci with Excessively High Depth ----------------------------
mean_depth <- mean(data_gl_filtered$other$loc.metrics$rdepth)
sd_depth <- sd(data_gl_filtered$other$loc.metrics$rdepth)
max_depth <- mean_depth + 3 * sd_depth
min_depth <- mean_depth - 3 * sd_depth

data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = min_depth, upper = max_depth, v = 3)
# After filtering: 185 individuals, 2615  loci

# Step 8: Remove Secondary Loci ----------------------------------------------
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method = "random", v = 3)
# After filtering: 185 individuals, 2386 loci

gl.report.reproducibility(data_gl_filtered)
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t = 0.95, v = 3)
# After filtering: 185 individuals, 2144  loci

# Step 9: Filter Loci Deviating from Hardy-Weinberg Equilibrium (HWE) -------
# (May not be suitable for selfing populations)
data_gl_filtered <- gl.filter.hwe(
  x = data_gl_filtered,
  subset = "each",              # Perform HWE filtering per population
  n.pop.threshold = 1,          # Filter SNPs that fail in at least 1 population
  test.type = "Exact",          # Use Exact test for small sample sizes
  mult.comp.adj = TRUE,         # Apply multiple testing correction
  mult.comp.adj.method = "fdr", # Use False Discovery Rate (FDR) correction
  alpha = 0.01,                 # Stricter p-value threshold for filtering
  n.min = 5,                    # Ignore populations with <5 individuals
  verbose = 2                   # Set verbosity for tracking
)
# After filtering: 185 individuals, 1348  loci


#Final filtering to match Premachandra 2019 paternity recommendation
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.99, v = 3)
data_gl_filtered <- gl.filter.maf(data_gl_filtered, threshold = 0.07, v = 3)
mean(data_gl_filtered@other$loc.metrics$maf)



# Final Recalculation of Metrics ---------------------------------------------

data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3)
#ind = 201     , loci = 1328 

#replace '5 with 05 labels
data_gl_filtered@other$ind.metrics$genotype <- gsub("(?<=_)5(?=$)", "05", data_gl_filtered@other$ind.metrics$genotype,
                                                    perl = TRUE) # Replace '_5' with '_05'
data_gl_filtered@ind.names <- gsub("(_5)(?![0-9])", "_05", data_gl_filtered@ind.names, perl = TRUE) # Replace '_5' with '_05' when not followed by a digit

data_gl_filtered@other$ind.metrics %>% group_by(stage) %>% summarise(n = n()) 



# Convert GENIND OBJECT all indiv
data_genind_pre <- gl2gi(data_gl_filtered)


# Filter for adults
adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)
data_gl_filtered_adult@ind.names
# 81 genotypes  
# Convert genind adults only
data_genind_adult <- gl2gi(data_gl_filtered_adult)




# cervus prep -------------------------------------------------------------

# Cervus Platy mapping extraction (working)------------------------------------------------------

##
####NOTE: Seems to be too many canidate reps n offspring file
data1 <- data_genind_pre@tab  # This assumes SNP data is stored in the 'tab' slot of the genind object
head(data1)

data1 = t(data1) %>% data.frame()
#data1$SNP = data_genind@other$loc.metrics$SNP
colnames(data1) <- gsub("\\.", "_", colnames(data1))
colnames(data1) <- ifelse(grepl("^X", colnames(data1)), colnames(data1), paste0("X", colnames(data1)))

data1$rownames <- rownames(data1)
#
data2 <- data1 %>% mutate(CloneID = sub("-.*", "", rownames), allele = sub(".*\\.(.)$", "\\1", rownames)) %>% 
  dplyr::select(CloneID, allele, everything()) %>% dplyr::select(-rownames) %>% arrange(., CloneID)
#data1 %>%  dplyr::select(., CloneID, pd13_l_14_11, pd13_a_1) %>% filter(., CloneID == check )  #

#this is correct, so there must be an issue with extracting from the genlight

#data2 = data1[8:16, 1:5]
data1_long = data2 %>% tidyr::pivot_longer(-c(CloneID, allele) ,  names_to = "id" ,values_to = "counts") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)


nrow(data1_long)  #91413
row_counts = data1_long %>%
  group_by(CloneID, id) %>%
  summarise(row_count = n()) %>%
  ungroup()
min(row_counts$row_count)

# Filter out groups that do not have exactly 2 rows
data1_long = data1_long %>%
  group_by(CloneID, id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  data.frame()
nrow(data1_long)  #81144
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe

#data1_long = data1_long[129:140,]

# Function to repeat allele with a separator
repeat_with_separator <- function(allele, counts, sep = "_") {
  if (is.na(allele) | is.na(counts)) {
    return("NA")
  } else {
    return(paste(rep(allele, counts), collapse = sep))
  }
}

# Apply the function to the dataframe
data1_long <- data1_long %>%
  mutate(allele_multiplied = mapply(repeat_with_separator, allele, counts, sep = ",")) # you can change sep to "_" if you want underscores
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe

# Group and summarise the data
summary_data <- data1_long %>%
  group_by(CloneID, id) %>%
  summarise(
    id = first(id),
    base = paste(allele_multiplied, collapse = ",")
  ) %>%
  arrange(id) %>%
  data.frame()

# Remove leading or trailing commas in summary_data$base
summary_data$base <- gsub("^,+|,+$", "", summary_data$base)


# Separate the 'base' column into 'base1' and 'base2'
summary_data <- summary_data %>%
  separate(base, into = c("a", "b"), sep = ",", extra = "drop", fill = "right")

data1_long = summary_data %>% tidyr::pivot_longer(-c(CloneID, id) ,  names_to = "base" ,values_to = "code") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)
data1_long$LocusID = paste0(data1_long$CloneID, data1_long$base)
#result <- data1_long %>% group_by(CloneID) %>% summarise(star_count = sum(code == "NA"))
result <- data1_long %>%
  group_by(CloneID) %>%
  summarise(missing_prop = sum(code == "NA") / n())
hist(result$missing_prop)
quantile(result$missing_prop)
nrow(data1_long)  #599382, 
length(unique(data1_long$CloneID))  #1491
good_loci = result %>% filter(missing_prop < 0.01)  # Adjust threshold as needed
nrow(good_loci)  #863
print(paste('You have', nrow(good_loci), 'good loci'))
data1_long <- data1_long %>%
  semi_join(good_loci, by = "CloneID")
length(unique(data1_long$CloneID))  #2621, 955
data1_long = data1_long %>% dplyr::select(-c(CloneID , base))
data1_long <- data1_long %>%
  mutate(code = ifelse(code == "NA", "*", code))


data_wide <- data1_long %>% tidyr::pivot_wider(names_from = LocusID, values_from = code, names_prefix = "X") %>% 
  data.frame()#year goes to columns, their areas go as the values, area is the prefix



write.csv(data_wide, row.names = FALSE,
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus",
                           "a_hyac_map_letters_code1_unknowns.csv"))


# offspring file ----------------------------------------------------------

df2 = data_gl_filtered$other$ind.metrics
df2$id = paste0('X', df2$id)
mismatches <- setdiff(colnames(data1), df2$id)

df2_lar = subset(df2, stage!= 'adults')
df2_adult = subset(df2, stage!= 'larvae')

df1a = df2_lar %>% dplyr::select(., id, genotype )


# Prepare the offspring IDs
offspring_ids <- df2_lar %>% dplyr::select(id)

# Prepare the candidate parent IDs
candidate_parents <- unique(df2_adult$id)
num_candidates <- length(candidate_parents)

# Create the offspring data frame
offspring_df <- offspring_ids

# Add candidate parent columns
for (i in 1:num_candidates) {
  offspring_df[[paste0("candidate_", i)]] <- candidate_parents[i]
}

offspring_df

# Write the CSV file
write.csv(offspring_df, row.names = FALSE,
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus",
                           "offspring_a_hyac_both_unknown.csv"))




