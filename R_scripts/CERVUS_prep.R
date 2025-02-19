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
                           "a_hyac_map_letters_code1_docent.csv"))


#DONT USE BELOW - SEE MOST RECENT HERON WRANGLING AS data_gl_filtered EXTRACTING CAUSED MISMATCHES

# #trying to extract from filtered object rather than raw (below) so already filtered.
# snp_data_list <- data_gl_filtered@gen
# snp_data_matrix <- do.call(cbind, lapply(snp_data_list, as.integer))
# data1 <- as.data.frame(snp_data_matrix)
# #rownames(data1) <- data_gl_filtered@loc.names  # Locus names as rownames
# colnames(data1) <- data_gl_filtered@ind.names  # Individual IDs as colnames
# data1$AlleleID = data_gl_filtered@loc.names
# data1$CloneID = data_gl_filtered@other$loc.metrics$CloneID
# data1$SNP = data_gl_filtered@other$loc.metrics$SNP
# data_gl_filtered@other$ind.metrics$stage
# head(data1)
# 
# 
# #wrangling
# names(data1) <- gsub("\\.", "_", names(data1))  #use underscores instead
# 
# # remove all eggs
# # columns_with_e <- grep("_e_", names(data1), value = FALSE)
# # data1 <- data1[, -columns_with_e]
# # head(data1)
# 
# 
# # clean and simplyfy dataframe
# colnames(data1)
# 
# data2 <- data1 %>%
#   mutate(LocusID = CloneID) %>%
# dplyr::select(-c(AlleleID , CloneID))
#   # dplyr::select(., c(LocusID, SNP), starts_with(c("pd")))
# 
# # Replacing '-' with NA and converting to character to prevent conversion to factor
# data2[data2 == "-"] <- NA
# data2 <- mutate_all(data2, as.character)
# data2$ref <- substr(data2$SNP, start = nchar(data2$SNP) - 2, stop = nchar(data2$SNP) - 2) # extract ref
# data2$var <- substr(data2$SNP, start = nchar(data2$SNP), stop = nchar(data2$SNP)) # extract var
# nrow(data2) # no of alleles
# length(unique(data2$LocusID)) # distint alleles
# head(data2)
# duplicated_rows <- data2[duplicated(data2$LocusID), ]
# # Remove duplicate rows based on LocusID.Not sure why duplicates but will remove
# data2 <- data2 %>% distinct(LocusID, .keep_all = TRUE)
# nrow(data2)
# 
# # long format prep for two columns
# data3a <- data2 %>%
#   pivot_longer(-c(SNP, LocusID, ref, var),
#     names_to = "sample",
#     values_to = "genotype"
#   ) %>%
#   mutate(., rowid = "a") %>%
#   data.frame()
# 
# data3b <- data2 %>%
#   pivot_longer(-c(SNP, LocusID, ref, var),
#                names_to = "sample",
#     values_to = "genotype"
#   ) %>%
#   mutate(., rowid = "b") %>%
#   data.frame()
# nrow(data3b)
# 
# #description: SNP 1 Row Mapping Format: "0" = Reference allele homozygote, "1" = SNP allele homozygote, "2"= heterozygote and "-" = double null/null allele homozygote (absence of fragment with SNP in genomic representation)
# data3a$base <- ifelse(data3a$genotype == 0, data3a$ref,
#                       ifelse(data3a$genotype == 1, data3a$var,
#                              ifelse(data3a$genotype == 2, data3a$ref,
#                                     data3a$genotype
#                              )
#                       )
# )
# data3a
# 
# data3b$base <- ifelse(data3b$genotype == 0, data3b$ref,
#                       ifelse(data3b$genotype == 1, data3b$var,
#                              ifelse(data3b$genotype == 2, data3b$var,
#                                     data3b$genotype
#                              )
#                       )
# )
# 
# data4 <- rbind(data3a, data3b) # join back
# data4$base <- ifelse(is.na(data4$base), 0, data4$base) # need missing values to be 0 for nalysis.
# head(data4)
# data4$stage = data_gl_filtered@other$ind.metrics$stage
# data5 <- data4 %>%
#   select(c(LocusID, sample, base, rowid, stage)) %>%
#   dplyr::arrange(., sample, LocusID)
# head(data5)
# unique(data5$sample)
# 
# data5$LocusID <- paste0(data5$LocusID, data5$rowid)
# data6 <- data5 %>% dplyr::select(., c(LocusID, sample, base, stage))
# data6$base
# data_wide <- data6 %>%
#   tidyr::pivot_wider(names_from = LocusID, values_from = base) %>%
#   data.frame()
# 
# rownames(data_wide) <- NULL
# head(data_wide)
# data_wide$sample
# (no_loc <- (ncol(data_wide) - 1) / 2) # no of distict locii
# # data_wide$sample <- gsub("\\.", "_", data_wide$sample)
# 
# # min typed loci
# typed_loci_per_individual <- apply(data_wide, 1, function(x) sum(x != 0))
# (min_typed_loci <- min(typed_loci_per_individual))
# data_wide <- data_wide[typed_loci_per_individual >= 500, ] # remove all <500 loci
# data_wide$sample
# 
# #add x to 
# data_wide$sample
# data_wide$sample = paste0('x', data_wide$sample)
# 
# # subsample to access run speed
# (no_loc <- (ncol(data_wide) - 1) / 2) # no of distict locii
# data_wide = data_wide[1:4002]  #CERVUS struggles around 8000. NOTE YOU SHOULD FILTER (RATHER THAN FIRST) THIS ON QUALITY LATER ON
# (no_loc <- (ncol(data_wide) - 2) / 2) # no of distict locii
# 
# # split larvae in ~half (there might be some issues running altogether)
# #data_wide1 <- data_wide[grep("_a_|^pd13|^pd14", data_wide$sample), ]
# #data_wide2 <- data_wide[grep("_a_|^pd5|^pd9|^pd15", data_wide$sample), ]
# #str(data_wide2)
# #nrow(data_wide2)
# #(ncol(data_wide2)-1)/2  #no of distict locii
# #data_wide2$sample
# # trimmed test
# # quarter_loc = ceiling(no_loc / 4)
# # data_wide = data_wide[, 1:(quarter_loc+1)]
# 
# 
# #old
# # write.csv(data_wide, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code.csv"))
# # write.csv(data_wide1, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code1.csv"))
# # write.csv(data_wide2, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code2.csv"))
# #write.table(data_wide, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code.txt"))
# 
# #new
# write.csv(data_wide, row.names = FALSE,
#           file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/4 Palau genetics mixing/Cervus",
#                            "a_hyac_map_letters_code1.csv"))
# 
# 
# 
# 
# 
# 
# 
# offspring file ----------------------------------------------------------


df2 = data_gl_filtered$other$ind.metrics
df2$id = paste0('X', df2$id)
mismatches <- setdiff(colnames(data1), df2$id)


df2_lar = subset(df2, stage!= 'adults')
df2_adult = subset(df2, stage!= 'larvae')


df1a = df2_lar %>% select(., id, genotype )
# df1 <- subset(data_wide, stage!= 'adults')  #remove factor treatment level. Use '%in%' to keep.
# df1$id = df1$sample
# df1a = left_join(df1, df2_lar, by  = 'id')


# Assuming df2_adult is your dataframe with columns 'genotype' and 'id'

# Extract the first replicate for each genotype
df2_adult$id
known_dam <- df2_adult %>%
  group_by(genotype) %>%
  slice(1) %>%
  ungroup() %>%
  select(genotype, id) %>%
  rename(known_dam = id) %>% data.frame()

# Merge the first replicate back into the original dataframe
df1b = left_join(df1a, known_dam, by = "genotype")
nrow(df1b)
df1 = select(df1b, c(id, known_dam))
df1 <- df1[complete.cases(df1), ] # make sure import matches NA type



# indices_with_a <- grep("_a_", data_wide$sample)
# labels_with_a <- data_wide$sample[indices_with_a]
# first_matching_label <- sapply(dam, function(dam) {
#   matching_labels <- grep(dam, labels_with_a, value = TRUE)
#   if (length(matching_labels) > 0) {
#     return(matching_labels[1])
#   } else {
#     return(NA)
#   }
# })
# df1$known_dam  = first_matching_label
# find_error_dam = which(df1$known_dam == 'pd5_a_1')
# df1[find_error_dam, 2] = 'pd13_a_1'

# not_known_dam <- labels_with_a[!labels_with_a %in% first_matching_label]
# length(not_known_dam)

#add candiatate parents
len = length(unique(df2_adult$id))
cands <- rep(df2_adult$id, length(df1$known_dam)) %>% sort(.)
cands_df <- matrix(cands, nrow = length(df1$known_dam), ncol = len, byrow = FALSE) %>% data.frame()
colnames(cands_df) <- rep("candidate", len)
candidate_indices <- which(colnames(cands_df) == "candidate")
colnames(cands_df)[candidate_indices] <- paste0("candidate_", seq_along(candidate_indices))

offspring_df <- cbind(df1, cands_df)
# offspring_df <- offspring_df %>%
#   mutate_at(vars(-1), ~ paste0("x", .))
#colnames(cands_df) <- rep("candidate", len)

# for (col in names(cands_df)) {
#   # Add leading zero to single-digit numbers at the end of the strings
#   offspring_df[[col]] <- gsub("_(\\d)$", "_0\\1", offspring_df[[col]])
#
#   # Add leading zero to single-digit numbers following 'pd'
#   offspring_df[[col]] <- gsub("pd(\\d)_", "pd0\\1_", offspring_df[[col]])
# }
offspring_df


#some reason i removed this
#offspring_df <- subset(offspring_df, offspring != "pd15_l_14_10") # remove factor treatment level. Use '%in%' to keep.


# split
#offspring_df1 <- offspring_df[grep("^pd13|^pd14", offspring_df$Offspring), ]
#offspring_df2 <- offspring_df[grep("^pd5|^pd9|^pd15", offspring_df$Offspring), ]

#str(offspring_df)
# write.csv(offspring_df, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "offspring_platy.csv"))
# write.csv(offspring_df1, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "offspring_platy1.csv"))
# write.csv(offspring_df2, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "offspring_platy2.csv"))

write.csv(offspring_df, row.names = FALSE, 
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus", 
                           "offspring_a_hyac.csv"))

