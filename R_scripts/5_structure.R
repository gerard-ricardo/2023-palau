# structure analysis (dart - working) ----------------------------------------------------
# tutorial https://green-striped-gecko.github.io/kioloa/session12.html


##notes
# - seems like you need min K = 1:(K+1) and num.k.rep = 2 to get delta K. So run across 1 higher and lower than you think K is. 


# data wrangling ----------------------------------------------------------
# data_gl_adult_unique <- gl.recalc.metrics(data_gl_adult_unique, v = 3)
# data_gl_adult_unique@other$ind.metrics
# data_gl_adult_unique$pop
# dim(data_gl_adult_unique) # Should return (42, 1340)
# length(indNames(data_gl_adult_unique))  # Should return 21
# colnames(as.matrix(data_gl_adult_unique))  # List all loci names
# length(colnames(as.matrix(data_gl_adult_unique)))  # Count number of loci
# indNames(data_gl_adult_unique)  # List all individuals
# length(indNames(data_gl_adult_unique))  # Count number of individuals
# table(is.na(pop(data_gl_adult_unique)))  # Count how many individuals have NA populations
# indNames(data_gl_adult_unique)[is.na(pop(data_gl_adult_unique))]  # List individuals with missing population data
# 
# 



# structure ---------------------------------------------------------------

#test run
tic("Running structure analysis") # start the timer with a message
struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 2000, numreps = 4000, k.range = 1:5, num.k.rep = 2, 
                                seed = 1, noadmix=FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
toc() 
#No admixture
# tic("Running structure analysis") # start the timer with a message
# struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 2000, numreps = 4000, k.range = 1:5, num.k.rep = 2, 
#                                 seed = 1, noadmix=TRUE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
# toc() 

#46 min burnin = 2000, numreps = 4000, k.range = 1:5, num.k.rep = 2
#save(struct_adult, file = file.path("./Rdata", "struct_adult_test.RData"))

# struct_all = gl.run.structure(data_gl_filtered, verbose = 3, burnin = 2000, numreps = 2000, k.range = 1:5, num.k.rep = 2, 
#                               seed = 1, noadmix=FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
#linux
# tic("Run Structure") # Start the timer
# struct_adult = gl.run.structure(data_gl_filtered_adult, verbose = 3, burnin = 1000, numreps = 1000, k.range = 2:5, num.k.rep = 2, 
#                                 seed = 1, noadmix=FALSE, exec = "/home/gricardo/structure/structure.exe")
# toc() # End the timer (170.7 sec   = is slower)

# formal run
tic("Running structure analysis") #
struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 8000, numreps = 20000, k.range = 2:4, num.k.rep = 2,
                                seed = 1, noadmix = FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
toc()
#seems to work for k.range = 2:4, num.k.rep = 2 (but not other settings)
#save(struct_adult, file = file.path("./Rdata", "struct_adult_1_10.RData"))
load("./Rdata/struct_adult_test.RData")
#seems to vary each time
str(struct_adult)
ev <- gl.evanno(struct_adult, plot.out = TRUE)
#high delta K is best K 
gl.evanno(struct_adult) # Get Evanno's table with likelihood values


#Deviation from Hardy-Windberg
data_gl_adult_unique1 = data_gl_adult_unique
pop(data_gl_adult_unique1) <- factor(rep("MergedPop", length(indNames(data_gl_adult_unique1))))
gl.report.hwe(data_gl_adult_unique1)
gl.filter.hwe(data_gl_adult_unique1, alpha = 0.05, mult.comp.adj = TRUE, mult.comp.adj.method = "fdr")
#no deviation after multiple comaprisons

#ev <- gl.evanno(struct_all, plot.out = TRUE)

qmat <- dartR::gl.plot.structure(struct_adult, K = 3, colors_clusters = list("dodgerblue", "mediumseagreen", "salmon", 'pink'), clumpak = T, save2tmp = T)
qmat <- dartR::gl.plot.structure(struct_adult, K = 2, colors_clusters = list("dodgerblue", "mediumseagreen"), clumpak = T, save2tmp = F, 
                                 plot.out = F)

head(qmat)
p3 = gl.print.reports(1)
gl.list.reports()
#p3$1

# colour seem to rotate. Need ind 1 = red, 12 = green, 4 = blue

# #Melt and append Q matrices
Q_melt <- do.call("rbind", lapply(qmat, reshape2::melt, id.vars = c("Label", "K", "orig.pop", "ord"), variable.name = "Cluster" ))

Q_melt$orig.pop <-
  factor(Q_melt$orig.pop, levels = unique(struct_adult[[1]]$q.mat$orig.pop))

p3 <- ggplot(Q_melt, aes_(x= ~ factor(ord), y = ~value, fill = ~Cluster)) +
  geom_col(color = "black", size = 0.25, width = 1) +
  facet_grid(K ~ orig.pop , scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    breaks = unique(Q_melt$ord), labels = unique(Q_melt$Label), expand = c(0, 0)) +
  scale_fill_manual(values = c("dodgerblue", "mediumseagreen", "salmon")) +
  theme_sleek2() +
  theme(panel.spacing = unit(0, "lines"), panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, angle = 90),
        #axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1),
        #axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank() ,
        legend.position = "none") +
  labs(x = 'Individuals', y = 'K=3')
p3

#save(p3, file = file.path("./Rdata", "structure_plot.RData"))
load("./Rdata/structure_plot.RData") #p3

#create a map showing groupings
gl.map.structure(qmat = qmat, x = data_gl_filtered_adult, K = 3, scalex = 1, scaley = 0.5)
