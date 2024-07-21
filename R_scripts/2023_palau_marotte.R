#2023 maroote processing

#libraries
library(ggplot2)
library(dplyr)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code



# 1 Import data -----------------------------------------------------------
# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data1 <- read.excel() # read clipboard from excel
# save(data1, file = file.path("./Rdata", "2023_marotte_raw.RData"))
load("./Rdata/2023_marotte_raw.RData") #data1

# Assuming the 'data' variable holds the dataframe with your dataset

# Calibration values from your dataset (replace these with your actual values)
zeropoint_x <- 17
zeropoint_y <- -4076
zeropoint_z <- 10
acccal_x <- 4047
acccal_y <- 4061
acccal_z <- 4053

# Correct the raw acceleration readings
data1$acc_x_corrected <- (data1$acc_x - zeropoint_x) / acccal_x
data1$acc_y_corrected <- (data1$acc_y - zeropoint_y) / acccal_y
data1$acc_z_corrected <- (data1$acc_z - zeropoint_z) / acccal_z

# Calculate tilt angles in degrees
data1$tilt_x <- atan2(data1$acc_x_corrected, data1$acc_z_corrected) * (180 / pi)
data1$tilt_y <- atan2(data1$acc_y_corrected, data1$acc_z_corrected) * (180 / pi)

data1$flow_direction_degrees <- (atan2(data1$acc_y_corrected, data1$acc_x_corrected) * (180 / pi) + 90) %% 360
median(data1$flow_direction_degrees, na.rm = T)

# Create a time sequence for plotting (adjust format according to your actual datetime format)
# Here I'm assuming your time is in the format 'hours:minutes.seconds'
# We need to convert this to a format that R can interpret as a time object
data1$datetime_formatted <- as.POSIXct(strptime(data1$datetime, format="%d/%m/%Y %H:%M:%S"), tz="UTC")
data1$datetime_formatted <- with_tz(data1$datetime_formatted, tzone = "Pacific/Palau")


# Plot tilt over time to observe changes




data_long <- pivot_longer(data1, cols = c("tilt_x", "tilt_y"), names_to = "tilt_direction", values_to = "tilt_value") %>% data.frame()
data_long = data_long[sample(nrow(data_long), 5000), ]


p0 <- ggplot()
p0 <- p0 + geom_point(data = data_long, aes(x = datetime_formatted , y = tilt_value, alpha = 0.8), color = "steelblue", size =  3, position = position_jitter(height = 0.01, width = .01))
p0 <- p0 + labs(
  x = expression(Time),
  y = expression(Tilt))
#p0 <- p0 + scale_y_continuous(limits = c(0, 1))
p0 <- p0 + theme_sleek2()
p0 <- p0 + scale_fill_manual(values = c(c("#3B9AB2", "#78B7C5", "#E1AF00", "#F21A00")), ) # can lookup html  hex colours to get specific eg #FFE6F5
p0 <- p0 + theme(legend.position = "none")
p0 <- p0 + facet_wrap(~tilt_direction, nrow = 2)
p0 <- p0 + guides(color = "none") # this removes excesse legends
p0

#split bytidal cycle
data_long <- data_long %>%
  mutate(tide = if_else(flow_direction_degrees > 200, 'flood', 'ebb'))

data_long %>% group_by(tide) %>% dplyr::summarise(med = median(flow_direction_degrees, na.rm =T)) %>% data.frame()  #basic groups


p0 <- ggplot()
p0 <- p0 + geom_point(data = data_long, aes(x = datetime_formatted , y = flow_direction_degrees, alpha = 0.8), color = "steelblue", size =  3, position = position_jitter(height = 0.01, width = .01))
p0 <- p0 + labs(
  x = expression(Time),
  y = expression(Direction))
#p0 <- p0 + scale_y_continuous(limits = c(0, 1))
p0 <- p0 + theme_sleek2()
p0 <- p0 + scale_fill_manual(values = c(c("#3B9AB2", "#78B7C5", "#E1AF00", "#F21A00")), ) # can lookup html  hex colours to get specific eg #FFE6F5
p0 <- p0 + theme(legend.position = "none")
p0 <- p0 + facet_wrap(~tide, nrow = 2)
p0 <- p0 + guides(color = "none") # this removes excesse legends
p0

# Save the plot if needed
#ggsave("tilt_plot.png", width = 8, height = 4)



# flow rates from maroote--------------------------------------------------------------

# read.excel <- function(header = TRUE, ...) {
#   read.table("clipboard", sep = "\t", header = header, na.strings = c("", "-", "na"), ...)
# }
# data1 <- read.excel() # read clipboard from excel
# data2 <- read.excel() # read clipboard from excel
# data3 = rbind(data1, data2)
# data1 = data3
# 
# save(data1, file = file.path("./Rdata", "2023_palau_flowrates.RData"))
load("./Rdata/2023_palau_flowrates.RData")  #data1

# 2 Labeling and wrangling -----------------------------------------------
str(data1) # check data type is correct
data1$datetime_formatted <- as.POSIXct(strptime(data1$datetime, format="%d/%m/%Y %H:%M:%S"), tz="UTC")

# clip of erros
cutoff_datetime <- as.POSIXct("12/04/2023 17:18:30", format = "%d/%m/%Y %H:%M:%S")
data1 <- data1 %>% filter(datetime_formatted <= cutoff_datetime)


# plotting
data1 = data1[sample(nrow(data1), 5000), ]


p0 <- ggplot()
p0 <- p0 + geom_point(data = data1, aes(x = datetime_formatted , y = speed_m_s , alpha = 0.8), color = "steelblue", size =  3, position = position_jitter(height = 0.01, width = .01))
p0 <- p0 + labs(
  x = expression(Time),
  y = expression(speed))
#p0 <- p0 + scale_y_continuous(limits = c(0, 1))
p0 <- p0 + theme_sleek2()
p0 <- p0 + scale_fill_manual(values = c(c("#3B9AB2", "#78B7C5", "#E1AF00", "#F21A00")), ) # can lookup html  hex colours to get specific eg #FFE6F5
p0 <- p0 + theme(legend.position = "none")
#p0 <- p0 + facet_wrap(~tilt_direction, nrow = 2)
p0 <- p0 + guides(color = "none") # this removes excesse legends
p0


p1 <- ggplot(data1, aes(x = speed_m_s)) +
  geom_density(aes(), alpha = 0.3) +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = speed_m_s), .width = c(.66, .95)) + theme_sleek2()



# flow rates from containers ----------------------------------------------

#x6_30 #   8:55 - 10:17
#5_05 = 8:48 - 9:55  = 67min, 912m
#4_30 = 8:49 - 10:11 = 82min, 978

df1 = data.frame(ID = c('6_30', '5_05', '4_30'), time = c(82, 67, 82), dist = c(942, 912, 978))
df1$time_s = df1$time * 60
df1$speed = df1$dist / df1$time_s
df1
mean(df1$speed )


#percentiles exceding
speeds_above_threshold <- data1$speed_m_s[data1$speed_m_s > 0.21]
percentiles_above_threshold <- quantile(speeds_above_threshold, probs = seq(0, 1, 0.01))
percentiles_above_threshold
# Calculate the percentage of values that exceed the threshold
percentage_above_threshold <- sum(data1$speed_m_s > 0.21) / length(data1$speed_m_s) * 100
percentage_above_threshold  #2.3% of values exceed what we observed on the night


