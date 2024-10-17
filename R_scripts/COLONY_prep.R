# COLONY formatting -------------------------------------------------------
data_genind_adult  #82 x 10717 

# (optional) filtering
hist(data_genind_adult@other$loc.metrics$CallRate)
hist(data_genind_adult@other$loc.metrics$coverage)
high_filt = which(data_genind_adult@other$loc.metrics$CallRate > 0.8 & data_genind_adult@other$loc.metrics$coverage > 20)
data_genind_adult <- data_genind_adult[, high_filt]  #82 x 4221


genotype_data_means <- data_genind_adult@other$loc.metrics %>%
  dplyr::select(CallRate:coverage) %>%  # Select the desired range of columns
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  # Use anonymous function to pass na.rm

genotype_data_means_sub1 <- data_genind_adult_subset1@other$loc.metrics %>%
  dplyr::select(CallRate:coverage) %>%  # Select the desired range of columns
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  # Use anonymous function to pass na.rm


# Calculating the genotyping error rate based on replicate consistency
genotyping_error_rate <- 1 - genotype_data_means$CallRate
med_genotyping_error_rate <- median(genotyping_error_rate, na.rm = TRUE)
#0, mean = 0.01



colony_data <- function(genind_obj) {
  
  input_name <- deparse(substitute(genind_obj))
  
  
  
  
  ## Calculating allelic dropout rate - these estimates are off because of subpop structure
  # Where allelic dropout is inferred from reduced FreqHets compared to homozygous frequencies
  # allelic_dropout_rate <- 1 - genotype_data$FreqHets / (genotype_data$FreqHomRef + genotype_data$FreqHomSnp + genotype_data$FreqHets) 
  # med_allelic_dropout_rate <- median(allelic_dropout_rate, na.rm = TRUE)
  
  hist(genind_obj@other$loc.metrics$AvgPIC )
  
  
  
  
  genotype_matrix <- as.matrix(tab(genind_obj))
  ncol(genotype_matrix)
  column_names <- colnames(genotype_matrix)
  locus_identifiers <- sapply(strsplit(column_names, "/"), `[`, 1)
  locus_counts <- table(locus_identifiers)
  valid_loci <- names(locus_counts[locus_counts == 2])
  valid_columns <- locus_identifiers %in% valid_loci
  genotype_matrix <- genotype_matrix[, valid_columns]
  ncol(genotype_matrix)
  num_columns_to_keep = ncol(genotype_matrix)
  
  #str(genotype_matrix)
  
  #str(genind_obj@other$loc.metrics)
  
  ##  make marker file ##
  (num_loci <- ncol(genotype_matrix)  /2)
  marker_names <- paste0("locus", 1:num_loci)
  
  marker_types <- rep(0, num_loci)
  
  # Use default error rates for allelic dropout and other errors (can be modified)
  allelic_dropout_rates <- rep(0.01, num_loci)   #based on REPave  being >99
  genotyping_error_rates <- rep(0.01, num_loci)  #checked
  
  # Create a matrix for easier handling
  marker_info <- rbind(
    marker_names,
    marker_types,
    allelic_dropout_rates,
    genotyping_error_rates
  )
  
  
  write.table(marker_info, file = paste0("./data/marker_info_", input_name, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  
  marker_info  #this creates a col for each loci. Should be half colony_data_lines -1
  ncol(marker_info)
  
  ## make offspring file ##
  
  # Extract individual names
  individual_names <- indNames(genind_obj)
  
  # Prepare for COLONY format: allele pairs per locus
  colony_data_lines <- character(nrow(genotype_matrix))
  for (i in 1:nrow(genotype_matrix)) {
    individual_data <- c(individual_names[i])  # Start with individual ID
    
    # Loop through loci, taking two columns at a time
    for (j in seq(1, ncol(genotype_matrix), by = 2)) {
      # Extract the two columns representing the dosage for alleles at each locus
      allele1_dosage <- genotype_matrix[i, j]
      allele2_dosage <- genotype_matrix[i, j + 1]
      
      # Determine the COLONY-compatible alleles based on dosage information
      if (is.na(allele1_dosage) || is.na(allele2_dosage)) {
        # Missing data
        alleles <- c(0, 0)
      } else if (allele1_dosage == 2 && allele2_dosage == 0) {
        # Homozygous for first allele, i.e., both alleles are type "C"
        alleles <- c(1, 1)
      } else if (allele1_dosage == 1 && allele2_dosage == 1) {
        # Heterozygous, one "C" and one "G"
        alleles <- c(1, 2)
      } else if (allele1_dosage == 0 && allele2_dosage == 2) {
        # Homozygous for second allele, i.e., both alleles are type "G"
        alleles <- c(2, 2)
      } else {
        # Handle any unexpected cases by assigning missing (failsafe)
        alleles <- c(0, 0)
      }
      
      # Add alleles to the individual's data
      individual_data <- c(individual_data, alleles)
    }
    
    # Combine all elements into a space-separated string
    colony_data_lines[i] <- paste(individual_data, collapse = " ")
  }
  
  # Write to a text file for COLONY input
  writeLines(colony_data_lines, con = paste0("./data/col_input_", input_name, ".txt"))
  
  #str(colony_data_lines)
  length(unlist(strsplit(colony_data_lines[1], " ")))   #this creates 200 +1 col. 2 col per loci
  ((length(unlist(strsplit(colony_data_lines[1], " ")))) - 1)/2
}

colony_data(data_genind_adult)  #1511 (after more filtering)
colony_data(data_genind_adult_subset1)
colony_data(data_genind_adult_subset2)

# split_lines <- strsplit(colony_data_lines, " ")
# colony_df <- do.call(rbind, lapply(split_lines, function(x) x))
# colony_df <- as.data.frame(colony_df, stringsAsFactors = FALSE)



## adjust marker info
data1 <- read.table(file = "./data/marker_info_data_genind_adult.txt", header = TRUE, dec = ",", na.strings = c("", ".", "na")) ## replace XXX to document name.
data2 <- read.table(file = "C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Colony/heron_sub1_inb_clones/heron_sub1_inb_clones.ErrorRate", 
                    header = TRUE, sep = ",", na.strings = c("", ".", "na")) ## replace XXX to document name.

error_rates <- data2 %>% select(MarkerID, DropRateEst)

data1_t <- as.data.frame(t(data1))
data1_t <- tibble::rownames_to_column(data1_t, "MarkerID")

updated_data <- data1_t %>% left_join(error_rates, by = "MarkerID") %>% mutate(V2 = ifelse(!is.na(DropRateEst), DropRateEst, V2)) %>% 
  select(-DropRateEst)


data1_updated <- as.data.frame(t(updated_data))
colnames(data1_updated) <- NULL
rownames(data1_updated) <- NULL
write.table(data1_updated, file = './data/marker_info_data_genind_adult_subset1_corr.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")


