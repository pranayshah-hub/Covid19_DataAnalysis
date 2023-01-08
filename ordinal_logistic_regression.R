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


################
# Loading Data #
################

# PSD Clinical Variables

psd_cc <- read.csv('PSD_COV_cc.csv')
psd_cc$covid_19_severity <- relevel(factor(psd_cc$covid_19_severity), ref='0')

drops <- c('results_date')
psd_cc <- psd_cc[, !(names(psd_cc) %in% drops)]

# PSD Clinical Variables After Mean Imputation

psd_cc_imp <- read.csv('psd_imp.csv')
psd_cc_imp$covid_19_severity <- relevel(factor(psd_cc_imp$covid_19_severity), ref='0')

psd_cc_imp <- psd_cc_imp[, !(names(psd_cc_imp) %in% drops)]

# PSD lncRNA

psd_lnc <- read.csv('psd_lnc_filt.csv')

drops <- c('key_0', 'X') #, 'SEQ1663', 'SEQ1466', 'SEQ0832', 'SEQ0930',
#'SEQ1174', 'SEQ2542', 'SEQ0424', 'SEQ1498', 'SEQ0310',
#'SEQ2864', 'SEQ2451', 'SEQ1629', 'SEQ0455', 'SEQ1727', 'SEQ2268', 
#'SEQ0312', 'SEQ0875', 'SEQ0530', 'SEQ1448', 'SEQ2811', 'SEQ2705',
#'SEQ2508', 'SEQ0045', 'SEQ2096', 'SEQ2609', 'SEQ2725', 'SEQ2575',
#'SEQ1163', 'SEQ1803')
psd_lnc <- psd_lnc[, !(names(psd_lnc) %in% drops)]

psd_lnc$covid_19_severity <- relevel(factor(psd_lnc$covid_19_severity), ref=1)


# PSD lncRNA After Mean Imputation

psd_lnc_imp <- read.csv('psd_lnc_filt_imp.csv')
drops <- c('key_0', 'X')#, 'SEQ1663', 'SEQ1466', 'SEQ0832', 'SEQ0930',
#'SEQ1174', 'SEQ2542', 'SEQ0424', 'SEQ1498', 'SEQ0310',
#'SEQ2864', 'SEQ2451', 'SEQ1629', 'SEQ0455', 'SEQ1727', 'SEQ2268', 
#'SEQ0312', 'SEQ0875', 'SEQ0530', 'SEQ1448')
psd_lnc_imp <- psd_lnc_imp[, !(names(psd_lnc_imp) %in% drops)]


#########
# Setup #
#########

reg <- list()
ctable <- list()
ps <- list() # data.frame(matrix(vector(), 5, 2))
or_ci <- list()


#############
# Functions #
#############

# Using the standard confint function doesn't work
polr_reg_psd <- function(data, reg, ctable, ps, or) {
    # Iterate over columns
    for (i in 2:length(colnames(data))) {
        if ((colnames(data)[i] != 'redacted')
            & (colnames(data)[i] != 'redacted')) {
            print(colnames(data)[i]) 
            
            # Fit model
            tryCatch({
            reg[[i]] <- polr(reformulate(colnames(data)[i],
                                             'redacted'),
                                 data=data, Hess=TRUE)
            print('Done')
            
            # Store coefficients, p-values and odds-ratios
            ctable[[i]] <- coef(summary(reg[[i]]))
            p <- pnorm(abs(ctable[[i]][ , "t value"]), 
                       lower.tail = FALSE) * 2
            ps[[i]] <- p
            ctable[[i]] <- cbind(ctable[[i]], "p value" = ps[[i]])
            or_ci[[i]] <- exp(cbind(OR=coef(reg[[i]]),
                                        confint(reg[[i]])))
            }, error=function(e){})
            
    # Return ctable
    return(ctable)
        }
    }
}


######################
# Fitting polr Model #
######################

# PSD Clinical
ctable_psd <- polr_reg_psd(psd_cc, reg, ctable, ps, or_ci)
ctable_psd <- do.call(rbind, ctable_psd)

# PSD Clinical After Mean Imputation
ctable_psd_imp <- polr_reg_psd(psd_cc_imp, reg, ctable, ps, or_ci)
ctable_psd_imp <- do.call(rbind, ctable_psd_imp)

# PSD LncRNA
ctable_psd_lnc <- polr_reg_psd(psd_lnc_filt, reg, ctable, ps, or_ci)
ctable_psd_lnc <- do.call(rbind, ctable_psd_lnc)

# PSD LncRNA After Mean Imputation
ctable_psd_lnc_imp <- polr_reg_psd(psd_lnc_filt_imp, reg, ctable, ps, or_ci)
ctable_psd_lnc_imp <- do.call(rbind, ctable_psd_lnc_imp)


##################
# saving Results #
##################

write.table(ctable_psd, 
            file='redacted/ctable_psd.csv', 
            sep=',')

write.table(ctable_psd_imp, 
            file='redacted/ctable_psd_imp.csv', 
            sep=',')

write.table(ctable_psd_lnc, 
            file='redacted/ctable_psd_lnc_sev.csv', 
            sep=',')

write.table(ctable_psd_lnc_imp, 
            file='redacted/ctable_psd_lnc_imp_sev.csv', 
            sep=',')

###############################################################################
###############################################################################


################
# Loading Data #
################

# LIH Ordinal Logistic Regression: CLINICAL

lih_c <- read.csv('LIH_PRD_c_clean.csv')

# Check levels
print(factor(lih_c[,'redacted']))
lih_c$NIH_classification_or <- factor(lih_c$NIH_classification_or, levels=c('redacted',
                                                                            'redacted',
                                                                            'redacted',
                                                                            'redacted'))

# Drop temp_mhyn for now
drops_lih <- cbind('redacted')
lih_c <- lih_c[, !(names(lih_c) %in% drops_lih)]


# LIH Ordinal Logistic Regression: CLINICAL After Mean Imputation

lih_c_imp <- read.csv('lih_imp.csv')

# Check levels
print(factor(lih_c_imp[,'redacted']))

lih_c_imp$NIH_classification_or <- factor(lih_c_imp$NIH_classification_or, levels=c('redacted',
                                                                                    'redacted',
                                                                                    'redacted',
                                                                                    'redacted'))

# Drop temp_mhyn for now
drops_lih_imp <- cbind('Unnamed..0', 'redacted')
lih_c_imp <- lih_c_imp[, !(names(lih_c_imp) %in% drops_lih_imp)]


# LIH LncRNA

lih_lnc <- read.csv('lih_lnc_filt.csv')

lih_lnc$NIH_classification_or <- factor(lih_lnc$NIH_classification_or, levels=c('redacted',
                                                                                'redacted',
                                                                                'redacted',
                                                                                'redacted'))
lih_lnc <- lih_lnc[, !(names(lih_lnc) %in% cbind('X', 'key_0'))]


# LIH lncRNA After Mean Imputation

lih_lnc_imp <- read.csv('lih_lnc_filt_imp.csv')

lih_lnc_imp$NIH_classification_or <- factor(lih_lnc_imp$NIH_classification_or, levels=c('redacted',
                                                                                        'redacted',
                                                                                        'redacted',
                                                                                        'redacted'))
lih_lnc_imp <- lih_lnc_imp[, !(names(lih_lnc_imp) %in% cbind('X', 'key_0'))]


#########
# Setup #
#########

reg <- list()
ctable <- list()
ps <- list() #data.frame(matrix(vector(), 5, 2))
or_ci <- list()


#############
# Functions #
#############

polr_reg_lih <- function(data, reg, ctable, ps, or, starting) { 
    for (i in starting:length(colnames(lih_c))) {
        if ((colnames(lih_c)[i] != 'redacted')
            & (colnames(lih_c)[i] == 'redacted')) {
            next
        }
        
        else if ((colnames(lih_c)[i] != 'redacted')
                 & (colnames(lih_c)[i] != 'redacted')) {
            print(colnames(lih_c)[i])
                
            # Fit model
            tryCatch({
            reg_lih[[i]] <- polr(reformulate(colnames(lih_c)[i],
                                             'redacted'),
                                 data=lih_c, Hess=TRUE)
            print('Done')
            
            # Store coefficients, p-values and odds-ratios
            ctable_lih[[i]] <- coef(summary(reg_lih[[i]]))
            p <- pnorm(abs(ctable_lih[[i]][, "t value"]), 
                       lower.tail = FALSE) * 2
            ps[[i]] <- p
            ctable_lih[[i]] <- cbind(ctable_lih[[i]], "p value" = ps[[i]])
            or_ci_lih[[i]] <- exp(cbind(OR=coef(reg_lih[[i]]),
                                        confint(reg_lih[[i]])))
            }, error=function(e){})
        }
    }
    # Return ctable
    return(ctable)
}


######################
# Fitting polr Model #
######################

ctable_lih <- polr_reg_lih(lih_c, reg, ctable, ps, or_ci, starting=2)
ctable_lih <- do.call(rbind, ctable_lih)

ctable_lih_imp <- polr_reg_lih(lih_c_imp, reg, ctable, ps, or_ci, starting=55)
ctable_lih_imp <- do.call(rbind, ctable_lih_imp)

ctable_lih_lnc <- polr_reg_lih(lih_lnc, reg, ctable, ps, or_ci, starting=1)
ctable_lih_lnc <- do.call(rbind, ctable_lih_lnc)

ctable_lih_lnc_imp <- polr_reg_lih(lih_lnc_imp, reg, ctable, ps, or_ci, starting=1)
ctable_lih_lnc_imp <- do.call(rbind, ctable_lih_lnc_imp)


##########
# Saving #
##########

write.table(ctable_lih, 
            file='redacted/ctable_lih.csv', 
            sep=',')

write.table(ctable_lih_imp, 
            file='redacted/ctable_lih_imp.csv', 
            sep=',')

write.table(ctable_lih_lnc, 
            file='redacted/ctable_lih_lnc.csv', 
            sep=',')

write.table(ctable_lih_lnc_imp, 
            file='redacted/ctable_lih_lnc_imp.csv', 
            sep=',')




