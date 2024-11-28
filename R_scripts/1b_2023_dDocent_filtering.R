## following the vcf link (vcf files complex to recreate from given inputs)
#try using https://rdrr.io/cran/dartR/man/gl2vcf.html



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
#                         ind.metafile = "./data/2023_palau_meta.csv", topskip = 6)
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

# Begin Filtering ------------------------------------------------------------

data_gl_filtered <- data_gl
#ind =  217  , loci = 50405 

# Initial Filtering with Call Rate and MAF ---------------------------

# 1a. Filter loci with call rate less than 50% (missing in more than 50% of individuals)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.5, v = 3) # Equivalent to --max-missing 0.5
#ind =  217  , loci = 30726 

# 1b. Filter loci with Minor Allele Frequency (MAF) less than 0.05
data_gl_filtered <- gl.filter.maf(data_gl_filtered, threshold = 0.05, v = 3)
# This is equivalent to --maf 0.05 in vcftools
#ind = 217   , loci = 18611 

# Step 2: Filter Genotypes by Minimum Read Depth -----------------------------

# Since individual genotype depth filtering is not directly available, we filter loci with mean read depth less than 3
# Assuming 'AvgDepth' represents average depth per locus
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = 3, v = 3)
# This approximates filtering genotypes with DP < 3
#ind = 217   , loci = 16540  

# Step 3: Remove Individuals with High Missing Data --------------------------

# Filter individuals with call rate less than 50% (more than 50% missing data)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind",
                                       threshold = 0.5, v = 3)
#ind = 204    , loci = 16540  

# Step 4: Apply Stricter Filters on Loci -------------------------------------

# 4a. Filter loci with call rate less than 95% (missing in more than 5% of individuals)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc",
                                       threshold = 0.95, v = 3)
#ind = 204    , loci = 2852  

# 4b. Ensure MAF ≥ 0.05 (already applied in Step 1b)

# 4c. Filter loci with mean depth less than 20
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = 20, v = 3)
#ind = 204    , loci = 1543 

# Step 5: Remove Loci with High Missing Data within Populations --------------

# # If population information is available in your metadata, assign populations
# # For example, if your metadata has a column 'Population'
# #data_gl_filtered <- gl.reassign.pop(data_gl_filtered, pop.assign = "Population")
# 
# # Calculate missing data per locus within each population
# # This requires custom processing since dartR doesn't have a direct function
# # We'll use the 'gl.filter.missloc' function and apply it per population
# 
# # Get list of populations
# populations <- levels(data_gl_filtered$pop)
# 
# # Initialize list to store loci to keep
# loci_to_keep <- list()
# 
# # Loop through each population and filter loci
# for (pop in populations) {
#   gl_pop <- gl.keep.pop(data_gl_filtered, pop.list = pop, mono.rm = FALSE, v = 0)
#   # Calculate locus call rate within the population
#   loci_callrate <- gl.report.callrate(gl_pop, method = "loc", v = 0)
#   # Identify loci with missing data less than 10% (call rate > 90%)
#   loci_keep <- names(which(loci_callrate > 0.9))
#   loci_to_keep[[pop]] <- loci_keep
# }
# 
# # Find loci that meet the threshold in all populations
# loci_keep_all_pops <- Reduce(intersect, loci_to_keep)
# 
# # Subset the data to keep only those loci
# data_gl_filtered <- data_gl_filtered[, loci_keep_all_pops]

# Step 6: Filter Based on Allele Balance (AB) --------------------------------

# Allele balance filters are specific to VCF INFO fields and may not be directly applicable
# If your data includes allele balance metrics, you can filter accordingly
# Since this may not be available, we'll skip this step

# Step 7-9: Filters Requiring Strand and Mapping Quality Information ---------

# These steps require detailed sequencing information not typically available in DArT data
# We will skip Steps 7-9 due to data limitations

# Step 10: Filter Based on Quality-to-Depth Ratio and Mean Depth -------------

# 10a. Filter loci where the quality score is less than 0.25 times the depth
# If quality scores are available in 'data_gl_filtered$other$loc.metrics$RepAvg'
# We can compute QUAL/DP and filter loci accordingly

# Calculate QUAL/DP ratio if possible  (not possible)
# if ("RepAvg" %in% colnames(data_gl_filtered$other$loc.metrics)) {
#   data_gl_filtered$other$loc.metrics$QUAL_DP_ratio <- data_gl_filtered$other$loc.metrics$RepAvg / data_gl_filtered$other$loc.metrics$coverage
#   # Filter loci where QUAL/DP < 0.25
#   loci_to_keep <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$QUAL_DP_ratio > 0.25)]
#   data_gl_filtered <- data_gl_filtered[, loci_to_keep]
# }

# 10b. Remove loci with excessively high depth (mean depth > mean + 3*SD)
hist(data_gl_filtered$other$loc.metrics$rdepth )
mean_depth = mean(data_gl_filtered$other$loc.metrics$rdepth)
sd_depth = sd(data_gl_filtered$other$loc.metrics$rdepth)
max_depth = mean_depth + 3 * sd_depth
min_depth = mean_depth - 3 * sd_depth

# Filter loci based on total read depth
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = min_depth, upper = max_depth, v = 3)
#ind = 204    , loci = 1543 


# Step 11: Filter Loci Deviating from Hardy-Weinberg Equilibrium (HWE) ------- (not sure if nessasry)

data_gl_filtered <- gl.filter.hwe(
  x = data_gl_filtered,         # The genlight object with SNP data
  subset = "all",               # Ignore population structure, treat as a single population
  n.pop.threshold = 1,          # Only one "population" considered, so set threshold to 1
  #method_sig = "Exact",         # Use the Exact test, which is more accurate for HWE
  #multi_comp = TRUE,            # Apply multiple testing correction
  #multi_comp_method = "fdr",    # Use False Discovery Rate (FDR) correction for multiple comparisons
  #alpha_val = 0.01,             # p-value threshold for HWE violations after FDR correction
  #min_sample_size = 5,          # Minimum sample size to apply the HWE test (irrelevant here with "all" but kept as default)
  verbose = 2                   # Set verbosity to track progress
)
#ind = 204    , loci = 908 


# # Perform HWE test per locus per population using 'pegas' package
# # Convert genlight object to genind
# data_genind <- gl2gi(data_gl_filtered)
# 
# # Initialize vector to store loci to remove
# loci_to_remove <- c()
# 
# # Perform HWE test for each population
# for (pop in populations) {
#   genind_pop <- seppop(data_genind)[[pop]]
#   # Compute HWE exact tests
#   hwe_results <- hw.test(genind_pop, B = 0) # B=0 for exact test
#   # Adjust p-values for multiple testing using Bonferroni correction
#   p_values <- hwe_results@p.value
#   p_adj <- p.adjust(p_values, method = "bonferroni")
#   # Identify loci with p-value less than 0.01
#   loci_violating_hwe <- names(p_adj[p_adj < 0.01])
#   # Add to loci to remove
#   loci_to_remove <- c(loci_to_remove, loci_violating_hwe)
# }
# 
# # Remove loci that deviate from HWE in any population
# loci_to_remove <- unique(loci_to_remove)
# data_gl_filtered <- data_gl_filtered[, !data_gl_filtered$loc.names %in% loci_to_remove]
# 
# # Step 12: Filter Based on Haplotype Structure -------------------------------
# 
# # Haplotype-based filtering requires phased genotype data and haplotype reconstruction
# # Since this is complex and may not be feasible with current data, we'll skip this step

# Final Recalculation of Metrics ---------------------------------------------

data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3)



# Filter for adults
adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)
data_gl_filtered_adult@ind.names
  

## BLAST prep
# Check the lengths of TrimmedSequence in @other$loc.metrics
trimmed_lengths <- nchar(as.character(data_gl_filtered@other$loc.metrics$TrimmedSequence))

# Filter for loci with trimmed sequences > 100 bp
filtered_loci <- trimmed_lengths > 68
data_gl_filtered_blast <- data_gl_filtered[filtered_loci, ]

# Check how many loci remain
print(paste("Number of loci > 68 bp:", sum(filtered_loci)))

gl2fasta(data_gl_filtered_blast, method=1, outpath = './data' , outfile='test.fasta',verbose=3)

