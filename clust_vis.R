#####################
# Loading Libraries #
#####################


library(foreign)
library(ggplot2)
library(reshape2)
library(readxl)
library(tidyverse)
library(NbClust)
library(FactoMineR)
library(factoextra)
library(data.table)
library(cluster)
library(ClustOfVar)

# TODO: Write the below in functions

########################
# PSD Cluster Analysis #
########################

# Read in data
psd_imp <- read.csv('psd_imp.csv')
psd_imp <- data.frame(psd_imp[,-1], row.names = psd_imp[,1])
names(psd_imp)

# Drop certain columns (redact these)
psd_imp <- psd_imp[, !(names(psd_imp) %in% cbind('X', 'Unnamed..0', 
                                                 'redacted', 
                                                 'redacted', 
                                                 'redacted',
                                                 'redacted', 
                                                 'redacted',
                                                 'redacted',
                                                 'redacted',
                                                 'redacted'))]
colnames(psd_imp)[colSums(is.na(psd_imp)) > 0]

# From FactoMineR: PCA
res_psd <- PCA(psd_imp[1:length(rownames(psd_imp)), 1:length(psd_imp)])
summary(res_psd)

# HCA
tree_psd <- hclustvar(psd_imp)
plot(tree_psd)

# Clustering
psd_impt <- as.data.frame(t(as.data.frame(psd_imp)))
colnames(psd_impt) <- psd_imp$A

kmvar8_psd <- kmeansvar(t(psd_impt), init=8)
summary(kmvar8_psd)    
kmvar9_psd <- kmeansvar(t(psd_impt), init=9)
summary(kmvar9_psd)    
kmvar10_psd <- kmeansvar(t(psd_impt), init=10)
summary(kmvar10_psd)    

write.csv(kmvar8_psd$cluster, 
          file='redacted/kmvar8_psd.csv')
write.csv(kmvar9_psd$cluster, 
          file='redacted/kmvar9_psd.csv')
write.csv(kmvar10_psd$cluster, 
          file='redacted/kmvar10_psd.csv')

# First find the optimal number of clusters
set.seed(31)
fviz_nbclust(psd_impt, kmeans, method='wss') + theme_minimal() + ggtitle("The Elbow Method")
# 2 is optimal?

optim_psd <- NbClust(psd_impt, 
                         distance = "euclidean", 
                         min.nc = 2, 
                         max.nc = 24,
                         method = "kmeans", index="kl") # 2 is optimal
optim_psd # 9 is optimal (for vars)

# Best number of clusters is 9
k_psd <- kmeans(psd_impt, centers = 9, nstart = 25)
fviz_cluster(k_psd,  data=psd_impt)





#########################
# UDUS Cluster Analysis #
#########################

# Read in Data
udus_imp <- read.csv('udus_imp.csv')
udus_imp <- data.frame(udus_imp[,-1], row.names = udus_imp[,1])

# Drop certain columns
udus_imp <- udus_imp[, !(names(udus_imp) %in% cbind('redacted', 'redacted'))]

udus_imp <- udus_imp[, apply(udus_imp, 2, var, na.rm=TRUE) != 0]

# PCA
res_udus <- PCA(udus_imp[1:length(rownames(udus_imp)), 1:length(udus_imp)])
summary(res_udus)

# HCA
tree_udus <- hclustvar(udus_imp)
plot(tree_udus)


# Clustering
udus_impt <- as.data.frame(t(as.data.frame(udus_imp)))
colnames(udus_impt) <- udus_imp$A

kmvar14_udus <- kmeansvar(t(udus_impt), init=14)
summary(kmvar14_udus)    
kmvar15_udus <- kmeansvar(t(udus_impt), init=15)
summary(kmvar15_udus)    
kmvar16_udus <- kmeansvar(t(udus_impt), init=16)
summary(kmvar16_udus)    

write.csv(kmvar14_udus$cluster, 
          file='redacted/kmvar14_udus.csv')
write.csv(kmvar15_udus$cluster, 
          file='redacted/kmvar15_udus.csv')
write.csv(kmvar16_udus$cluster, 
          file='redacted/kmvar16_udus.csv')

set.seed(31)
fviz_nbclust(udus_impt, kmeans, method='wss', k.max=16) + theme_minimal() + ggtitle("The Elbow Method")
# 15 is optimal

optim_udus <- NbClust(udus_impt, 
                          distance = "euclidean", 
                          min.nc = 2, 
                          max.nc = 16,
                          method = "kmeans", index="kl")
optim_udus # 15

# Best number of clusters is 15 for vars
k_udus <- kmeans(udus_impt, centers = 15, nstart = 25)
fviz_cluster(k_udus, data=udus_impt)





########################
# LIH Cluster Analysis #
########################

# Read in Data
lih_imp <- read.csv('lih_imp.csv')
lih_imp <- data.frame(lih_imp[,-1], row.names = lih_imp[,1])

# Drop certain columns
lih_imp <- lih_imp[, !(names(lih_imp) %in% cbind('redacted', 'Unnamed..0'))]

# PCA
res_lih <- PCA(lih_imp[1:length(rownames(lih_imp)), 1:length(lih_imp)])

# HCA
tree_lih <- hclustvar(lih_imp)
plot(tree_lih)

# Clustering
lih_impt <- as.data.frame(t(as.data.frame(lih_imp)))
colnames(lih_impt) <- lih_imp$A

kmvar9_lih <- kmeansvar(t(lih_impt), init=9)
summary(kmvar9_lih)    
kmvar10_lih <- kmeansvar(t(lih_impt), init=10)
summary(kmvar10_lih)    
kmvar11_lih <- kmeansvar(t(lih_impt), init=11)
summary(kmvar11_lih)    

write.csv(kmvar9_lih$cluster, 
          file='redacted/kmvar9_lih.csv')
write.csv(kmvar10_lih$cluster, 
          file='redacted/kmvar10_lih.csv')
write.csv(kmvar11_lih$cluster, 
          file='redacted/kmvar11_lih.csv')


set.seed(31)
fviz_nbclust(lih_impt, kmeans, method='wss', k.max=24) + theme_minimal() + ggtitle("The Elbow Method")
# 10 is optimal

optim_lih <- NbClust(lih_impt, 
                          distance = "euclidean", 
                          min.nc = 2, 
                          max.nc = 24,
                          method = "kmeans", index="kl")

optim_lih

# Best number of clusters is 10
k_lih <- kmeans(lih_impt, centers = 10, nstart = 25)
fviz_cluster(k_lih, data=lih_impt)






###############################################################################
###############################################################################
 
###########################
# LncRNA Cluster Analysis #
###########################

# PSD
psd_lnc <- read.csv('psd_lnc_filt_imp.csv')
names(psd_lnc)
psd_lnc <- psd_lnc[, !(names(psd_lnc) %in% cbind('X', 'Unnamed..0', 
                                                 'Unnamed..0.1', 
                                                 'redacted', 
                                                 'redacted'))]

psd_lnct <- as.data.frame(t(as.data.frame(psd_lnc)))
colnames(psd_lnct) <- psd_lnc$A

kmvar2_psd_lnc <- kmeansvar(t(psd_lnct), init=2)
summary(kmvar2_psd_lnc)    
kmvar3_psd_lnc <- kmeansvar(t(psd_lnct), init=3)
summary(kmvar3_psd_lnc)    
kmvar4_psd_lnc <- kmeansvar(t(psd_lnct), init=4)
summary(kmvar4_psd_lnc)    

# Save results
write.csv(kmvar2_psd_lnc$cluster, 
          file='redacted/kmvar2_psd_lnc.csv')
write.csv(kmvar3_psd_lnc$cluster, 
          file='redacted/kmvar3_psd_lnc.csv')
write.csv(kmvar4_psd_lnc$cluster, 
          file='redacted/kmvar4_psd_lnc.csv')

# First find the optimal number of clusters
set.seed(31)
fviz_nbclust(psd_lnct, kmeans, method='wss', k.max=24) + theme_minimal() + ggtitle("The Elbow Method")

optim_psd_lnc <- NbClust(psd_lnct, 
                     distance = "euclidean", 
                     min.nc = 2, 
                     max.nc = 26,
                     method = "kmeans", index="kl") 
optim_psd_lnc # 2 is optimal





# UDUS
udus_lnc <- read.csv('udus_lnc_filt_imp.csv')
names(udus_lnc)
udus_lnc <- udus_lnc[, !(names(udus_lnc) %in% cbind('X', 'Unnamed..0', 
                                                 'redacted', 
                                                 'redacted'))]

udus_lnct <- as.data.frame(t(as.data.frame(udus_lnc)))
colnames(udus_lnct) <- udus_lnc$A

set.seed(31)
fviz_nbclust(udus_lnct, kmeans, method='wss', k.max=24) + theme_minimal() + ggtitle("The Elbow Method")

optim_udus_lnc <- NbClust(udus_lnct, 
                         distance = "euclidean", 
                         min.nc = 2, 
                         max.nc = 26,
                         method = "kmeans", index="kl") 
optim_udus_lnc # 9 is optimal

kmvar9_udus_lnc <- kmeansvar(t(udus_lnct), init=9)
summary(kmvar9_udus_lnc)    
kmvar10_udus_lnc <- kmeansvar(t(udus_lnct), init=10)
summary(kmvar10_udus_lnc)    
kmvar11_udus_lnc <- kmeansvar(t(udus_lnct), init=11)
summary(kmvar11_udus_lnc)    

# Save results
write.csv(kmvar9_udus_lnc$cluster, 
          file='redacted/kmvar9_udus_lnc.csv')
write.csv(kmvar10_udus_lnc$cluster, 
          file='Credacted/kmvar10_udus_lnc.csv')
write.csv(kmvar11_udus_lnc$cluster, 
          file='redacted/kmvar11_udus_lnc.csv')




# LIH
lih_lnc <- read.csv('lih_lnc_imp.csv')
lih_lnc <- lih_lnc[, !(names(lih_lnc) %in% cbind('X', 'key_0', 
                                                 'redacted'))]

lih_lnc <- lih_lnc[, apply(lih_lnc, 2, var, na.rm=TRUE) != 0]

lih_lnct <- as.data.frame(t(as.data.frame(lih_lnc)))
colnames(lih_lnct) <- lih_lnc$A

set.seed(31)
fviz_nbclust(lih_lnct, kmeans, method = "wss", k.max = 26) + theme_minimal() + ggtitle("The Elbow Method")

optim_lih_lnc <- NbClust(lih_lnct, 
                          distance = "euclidean", 
                          min.nc = 2, 
                          max.nc = 26,
                          method = "kmeans", index="kl") 
optim_lih_lnc # 24 is optimal

kmvar24_lih_lnc <- kmeansvar(t(lih_lnct), init=24)
summary(kmvar24_lih_lnc)    
kmvar26_lih_lnc <- kmeansvar(t(lih_lnct), init=26)
summary(kmvar26_lih_lnc)    
kmvar28_lih_lnc <- kmeansvar(t(lih_lnct), init=28)
summary(kmvar28_lih_lnc)    

# Save results 
write.csv(kmvar24_lih_lnc$cluster, 
          file='redacted/kmvar24_lih_lnc.csv')
write.csv(kmvar26_lih_lnc$cluster, 
          file='redacted/kmvar26_lih_lnc.csv')
write.csv(kmvar28_lih_lnc$cluster, 
          file='redacted/kmvar28_lih_lnc.csv')














