# 2023 palau colony sizes

# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code


# 1 Import data -----------------------------------------------------------
# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data1 <- read.excel() # read clipboard from excel
# save(data1, file = file.path("./Rdata", "2023_palau_col_size.RData"))
load("./Rdata/2023_palau_col_size.RData")


# 2. Labelling -----------------------------------------------
str(data1) # check data type is correct
data1$ID <- as.factor(as.character(data1$ID))



# 3. Wrangling ------------------------------------------------------------
data1
#data1 <- data1[which(data1$notes != 'dup'),]  #remove dup rows
data1 <- data1[which(data1$notes == 'coral'),]  #remove dup rows

duplicated(data1$ID)


data1 <- data1 %>%
  group_by(ID) %>%
  summarise(total_mean_dia = sum(mean_dia, na.rm = TRUE))  %>% data.frame()

mean(data1$total_mean_dia)
median(data1$total_mean_dia)

data2 = data1[grep("c", data1$ID),]   #find rows in column that meet search term
mean(data2$total_mean_dia)
