# 2023 Acropora hyacinthus ------------------------------------------------

# load libraries ----------------------------------------------------------
#install.packages('dartr')
library(dartR)
library(dartR.popgen)
library(PopGenReport)
library(adegenet)
library(tictoc)
library(HardyWeinberg)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggrepel)
library(hierfstat)
library(ape)
library(poppr)
library(pegas)
library(dbscan)
library(sp)
library(rgdal)
library(clustertend)
library(cluster)
library(plotly)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code

# import data -------------------------------------------------------------

#
#meta_2023_acro <- read.csv("./data/Report_DAc24-9371_SNP_2_copy.csv", head = T, skip = 6) # make sure samples are in same order as in data_gl
data_gl <- gl.read.dart(filename = "./data/Report_DAc24-9371_SNP_2_copy.csv", ind.metafile = "./data/2023_palau_meta.csv", topskip = 6)
#Ids for individual metadata does not match the number of ids in the SNP data file. Maybe this is fine if a subset matches.
#ind.metafile ids not matched were:
# [1] "40"  "137" "196"

# recalculate metrics
data_gl <- gl.recalc.metrics(data_gl, v = 3) # recalculate loci metrics

save(data_gl, file = file.path("./Rdata", "2023_Acro_hyac_gl.RData"))
load("./Rdata/2023_Acro_hyac_gl.RData")  #data_gl


data_gl$other$loc.metrics
data_gl$other$ind.metrics
data_gl$other$loc.metrics$coverage <- data_gl$other$loc.metrics$AvgCountRef + data_gl$other$loc.metrics$AvgCountSnp
mean(data_gl$other$loc.metrics$coverage) # AH = 35.8
min((data_gl$other$loc.metrics$coverage)) #  AH = 5
max((data_gl$other$loc.metrics$coverage)) #  AH = 372
sd(data_gl$other$loc.metrics$coverage) / sqrt(1996) #AH = 0.66

# data filtering ----------------------------------------------------------
data_gl_filtered <- data_gl
#following the dartr suggested order (see tut5)
# gl <-gl.filter.secondaries(gl)
# gl <- gl.filter.rdepth(gl)
# gl <- gl.filter.reproducibility(gl)
# gl <- gl.filter.callrate(gl, method=”loc”)
# gl <- gl.filter.callrate(gl, method=”ind”)
# gl <- gl.filter.monomorphs(gl)

## pre-filtering
#indiv = 218  ind, 50405  loc

##secondaries
gl.report.secondaries(data_gl_filtered)
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method="random", verbose = 3) #remove loci fragment that shared SNPs. Only keep 1
#indiv = 218  ind, 34552   loc

#rdepth
gl.report.rdepth(data_gl_filtered)
#platy = 13.8 Generally 10 is considered min
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered,  lower = 10, v = 3) # filter by loci callrate
#indiv = 218 ind, 21260   loc, med depth = 23

##reproducibility 
gl.report.reproducibility(data_gl_filtered )
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.95, v=3) #filter out loci with limited reproducibility
# at 95%: ind = 218, loci = 20276   

# callrate loci (non missing data)
gl.report.callrate(data_gl_filtered, method = "loc") 
# AH = 0.59 (41% missing data)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
#at 70% AH = 218  ind, 7993  loci

#Minor Allele Frequency (MAF) and Coverage Filter:
list.match <- data_gl_filtered$loc.names[
  which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.01 & 
          data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef > 0.01 & 
          data_gl_filtered$other$loc.metrics$coverage > 4)
]
data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)]
#at 70% AH = 218  ind, 5678   loci

## call rate ind (non missing data). low could indicate poor extract or reference genome or contamination.
#individauls
gl.report.callrate(data_gl_filtered, method = "ind") 
# platy = 87%
pre_filt_ind <- data_gl_filtered@ind.names
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.7, v = 3) # filter by ind callrate
filt_ind <- data_gl_filtered@ind.names
(lost_ind <- setdiff(pre_filt_ind, filt_ind))
# at 70% AH = 204  ind, 5678   loci
#just a few larvae lost

data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3) # recalculate loci metrics


##others - not sure if needed
# data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
# not sure if I need HWE filter because remove inbreeding
# data_gl_filtered <- gl.filter.hwe(data_gl_filtered, alpha_val = 0.05, subset = "each", multi_comp_method = 'bonferroni',v=3) #filter out loci that depart from H-W proportions
# list.match <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.05 & data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef > 0.05 & data_gl_filtered$other$loc.metrics$coverage > 5)] #remove loci based on minor allele frequency and low data coverage
# data_gl_filtered <- data_gl_filtered[,match(list.match, data_gl_filtered$loc.names)]#keep only loci in the list above

# population filtering and objects ----------------------------------------


# look into genotype as population
data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "id")
data_gl_filtered


# Convert GENIND OBJECT all indiv
data_genind <- gl2gi(data_gl_filtered)

# Filter out eggs and larvae to keep only adults
adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)

# Convert genind adults only
data_genind_adult <- gl2gi(data_gl_filtered_adult)



