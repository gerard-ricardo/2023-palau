


library(dplyr)

# 1 Import data -----------------------------------------------------------
# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data1 <- read.excel() # read clipboard from excel
#save(data1, file = file.path("./Rdata", "2023_palau_dart_pos.RData"))
load("./Rdata/2023_palau_dart_pos.RData")

unique(data1$Genotype)
# This will return TRUE for all duplicates, including the first occurrence
duplicate_indices <- duplicated(data1$Genotype) | duplicated(data1$Genotype, fromLast = TRUE)

# To get the indices (row numbers) of these duplicates
duplicate_row_indices <- which(duplicate_indices)
print(duplicate_row_indices)
duplicate_rows <- data1[duplicate_indices, ]
duplicate_rows

which(data1$Genotype == '61')


## join the position ID to simple ID
# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data4_2 <- read.excel() # read clipboard from excel
# data5 = rbind(data3, data4, data4_2)
# data5 = dplyr::arrange(data5, geno) %>%  mutate(Genotype = as.character(geno)) #dplyr - use this. Allows multiple sorts i.e site then orient
#save(data5, file = file.path("./Rdata", "2023_palau_IDs.RData"))
load("./Rdata/2023_palau_IDs.RData")

data6 = left_join(data1, data5, by = 'Genotype')  #joining and keeping left

data7 <-data6[data6$Comment =='larvae',]   #remove numeric treatment level
table(data7$id)

data7[which(data7$id == '5_30'), ]

##########
data2 <-data1[data1$PlateID ==3,]   #remove numeric treatment level


# Define the ranges for spokes, distances, and replicates
spokes <- 1:7
distances <- c(5, 3, 10, 30)
replicates <- 1:2

# Generate all expected combinations
expected_combinations <- expand.grid(spokes = spokes, 
                                     distance = distances, 
                                     replicate = replicates)

# Create a vector of strings as they should appear in the Genotype column
expected_ids <- apply(expected_combinations, 1, function(x) paste(x, collapse = "_"))

# Check which of these are in your actual data
actual_ids <- as.character(data2$Genotype)  # Ensure it's character, if not already
presence_check <- expected_ids %in% actual_ids

# Identify which IDs are missing
(missing_ids <- expected_ids[!presence_check])


which(data2$Genotype == '4_5_1')
