# 2023 palau colony sizes

# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek3") # set theme in code


# Import data -----------------------------------------------------------
# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data1 <- read.excel() # read clipboard from excel
# save(data1, file = file.path("./Rdata", "2023_palau_col_size.RData"))
load("./Rdata/2023_palau_col_size.RData")


# Labelling -----------------------------------------------
str(data1) # check data type is correct
data1$ID <- as.factor(as.character(data1$ID))


# Wrangling ------------------------------------------------------------
data1
#data1 <- data1[which(data1$notes != 'dup'),]  #remove dup rows
data1 <- data1[which(data1$notes == 'coral'),]  #remove dup rows

duplicated(data1$ID)


data1 <- data1 %>%
  group_by(ID) %>%
  summarise(total_mean_dia = sum(mean_dia, na.rm = TRUE))  %>% data.frame()

mean(data1$total_mean_dia)
median(data1$total_mean_dia)

#centre
data2 = data1[grep("c", data1$ID),]   
mean(data2$total_mean_dia)
median(data2$total_mean_dia)

#spokes
data3 = data1[!grepl("c", data1$ID),]    
mean(data3$total_mean_dia)
median(data3$total_mean_dia)


# density
p1 <- ggplot(data1) +
  geom_density(aes(total_mean_dia), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = total_mean_dia), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
# geom_vline(xintercept = 0, color = "red", lty = 2)+ theme_sleek2()
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.05))
p1 <- p1 + scale_x_continuous(name = "Colony diameter (cm) ")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1

# all together
data1$group <- "overall"
data2$group <- "centre"
data3$group <- "spokes"
combined_data <- bind_rows(data1, data2, data3)

cols <- c("overall" = "#A2CFE3",  # pastel blue
                    "centre" = "#E3A2A8",   # pastel red
                    "spokes" = "#C6A2E3")   # pastel purple

# Density plot with all three distributions
p1 <- ggplot(combined_data, aes(x = total_mean_dia, fill = group, color = group)) +
  geom_density(alpha = 0.3) +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = total_mean_dia), .width = c(.66, .95)) +
  coord_cartesian(ylim = c(0.0, 0.08)) +
  scale_x_continuous(name = "Colony diameter (cm)") +
  scale_y_continuous(name = "Frequency") +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(fill = "Group", color = "Group") +
  theme_sleek3()+
  theme(legend.position = c(0.9, 0.9), legend.text = element_text(size = rel(1), colour = "grey20"))
p1

