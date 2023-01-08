#####################
# Loading Libraries #
#####################


library(foreign)
library(ggplot2)
library(nnet)
library(reshape2)
library(readxl)
library(tidyverse)
library(MASS)
library(dplyr)

# TODO: Write the below in functions

################
# MVA analysis #
################

# PSD

psd_lnc <- read.csv('psd_lnc.csv')
psd_lnc_pred <- psd_lnc[, !(names(psd_lnc) %in% cbind('X', 'key_0',
                                                 'redacted',
                                                 'redacted'))]

mva_psd <- psd_lnc_pred %>%
    summarize(means=lapply(psd_lnc_pred, function(x) avg=mean(x)), 
              stds=lapply(psd_lnc_pred, function(x) stds=sd(x))) 

mva_psd <- arrange(mva_psd, stds)
mva_psd$means <- as.numeric(mva_psd$means)
mva_psd$stds <- as.numeric(mva_psd$stds)


write.csv(mva_psd, 
          file='redacted/mva_psd.csv')


ggplot(mva_psd, aes(x=sds, y=means)) + geom_point(mapping = aes(x=means, y=stds))


# UDUS

udus_lnc <- read.csv('udus_lnc.csv')
udus_lnc_pred <- udus_lnc[, !(names(udus_lnc) %in% cbind('X', 'key_0',
                                                 'redacted',
                                                 'redacted'))]

mva_udus <- udus_lnc_pred %>%
    summarize(means=lapply(udus_lnc_pred, function(x) avg=mean(x)), 
              stds=lapply(udus_lnc_pred, function(x) stds=sd(x))) 

mva_udus <- arrange(mva_udus, stds)
mva_udus$means <- as.numeric(mva_udus$means)
mva_udus$stds <- as.numeric(mva_udus$stds)


write.csv(mva_udus, 
          file='redacted/mva_udus.csv')

ggplot(mva_udus, aes(x=sds, y=means)) + geom_point(mapping = aes(x=means, y=stds))



# LIH

lih_lnc <- read.csv('lih_lnc.csv')
lih_lnc <- lih_lnc[, !(names(lih_lnc) %in% cbind('X', 'key_0',
                                                      'redacted'))]


mva_lih <- lih_lnc_pred %>%
    summarize(means=lapply(lih_lnc_pred, function(x) avg=mean(x)), 
              stds=lapply(lih_lnc_pred, function(x) stds=sd(x))) 

mva_lih <- arrange(mva_lih, stds)
mva_lih$means <- as.numeric(mva_lih$means)
mva_lih$stds <- as.numeric(mva_lih$stds)

mva_lih

write.csv(mva_lih, 
          file='redacted/mva_lih.csv')


ggplot(mva_lih, aes(x=sds, y=means)) + geom_point(mapping = aes(x=means, y=stds))





variation <- function(df, name) {
    mva <- df %>%
        summarize(means=lapply(df, function(x) avg=mean(x)), 
                  stds=lapply(df, function(x) stds=sd(x))) 
    
    mva <- arrange(mva, stds)
    mva$means <- as.numeric(mva$means)
    mva$stds <- as.numeric(mva$stds)
    
    write.csv(mva, 
              file=paste('redacted',
                         name))
    
}