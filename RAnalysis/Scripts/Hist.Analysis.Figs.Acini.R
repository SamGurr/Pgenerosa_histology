#Title: Pgenerosa histology (Acini)
#Project: FFAR
#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200729
#See Readme file for details

rm(list=ls())

# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(tidyr)
library(ggpubr)
library(car)

#set working directory---------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_histology/RAnalysis/Data") #set working

# UPLOAD DATA------------------------------------------------------------------------------------------
dat<-read.csv("Master_summary_male_acini.csv", header=T, sep=",", na.string="NA", as.is=T) 
dat # view data
dat$Area = as.numeric(dat$Area)

# PIVOT THE TABLE TO CALC RELATIVE VALUES--------------------------------------------------------------
dat2 <- dat %>% dplyr::select(-c('Meas_num','Label','hue','saturation','brightness')) %>% 
  tidyr::pivot_wider(names_from=type, values_from=Area)

# CALC RELATIVE VALUES---------------------------------------------------------------------------------
dat2$TOTAL_AREA <- (dat2$cytes_zoa + dat2$lumen)
dat2$perc_zoa <- (dat2$zoa/dat2$TOTAL_AREA)*100 # percent area of spermatozoa
dat2$perc_cytes <- ((dat2$cytes_zoa - dat2$zoa)/dat2$TOTAL_AREA)*100 # percent area of spermatocytes
dat2$perc_lumen <- (dat2$lumen/dat2$TOTAL_AREA)*100 # percent area of lumen
dat2$zoa_cyte_ratio <- (dat2$zoa/(dat2$cytes_zoa - dat2$zoa)) # ratio of spermatozoa : spermatocytes
dat2$cytes <- (dat2$cytes_zoa - dat2$zoa) # area of spermatocytes

# CALC THE MEAN FOR EACH SAMPLE-------------------------------------------------------------------------
Means_Table <- dat2 %>%
  dplyr::select(-'acini_segment') %>% # remove unecessary columns
  group_by(ID,Treatment,Date) %>%
  summarize(
            mean_zoa= mean(zoa, na.rm = TRUE),
            mean_cytes = mean(cytes, na.rm = TRUE),
            mean_lumen= mean(lumen, na.rm = TRUE),
            mean_perc_zoa= mean(perc_zoa, na.rm = TRUE),
            mean_perc_cytes = mean(perc_cytes, na.rm = TRUE),
            mean_perc_lumen= mean(perc_lumen, na.rm = TRUE),
            mean_zoa_cyte_ratio = mean(zoa_cyte_ratio, na.rm = TRUE),
            num_acini=n())
Means_Table # view the table

# PIVOT LONGER-----------------------------------------------------------------------------------------
Means_Table_long <- Means_Table %>% 
  tidyr::pivot_longer(cols = c(4:10), names_to='scoring_metric', values_to='means')

#  PLOTS----------------------------------------------------------------------------------------
Means_Table_long$Date_Treat <- paste((substr(Means_Table_long$Date,5,8)),Means_Table_long$Treatment, sep ="_")

plot <- ggboxplot(Means_Table_long, x = "Date_Treat", y = "means",  fill = "Treatment",
                     palette = c("#FC4E07","#00AFBB"), add = "none")
plot2 <- plot %>% ggadd(shape ="Treatment",fill = "white") %>% ggadd("jitter", size = 4,shape ="Treatment",fill = "white") +
  facet_wrap( ~ scoring_metric, ncol=2, scales = "free") + theme_classic()
plot2

# STATISTICS------------------------------------------------------------------------------------------
Means_Table$Treatment <- as.factor(Means_Table$Treatment)
Means_Table$Date <- as.factor(Means_Table$Date)
par(mfrow=c(1,4)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots

# Two-Way ANOVA (treat×time)

## cytes test---------------------------------
mod_cytes_TIME  <- aov(mean_cytes~Treatment*Date, data = Means_Table)
anova(mod_cytes_TIME) # p = 0.9391
shapiro.test(residuals(mod_cytes_TIME)) #  normal residuals p-value = 7.261e-06
leveneTest(mod_cytes_TIME) # p = 0.6463
hist(residuals(mod_cytes_TIME)) #plot histogram of residuals
boxplot(residuals(mod_cytes_TIME)) #plot boxplot of residuals
plot(fitted(mod_cytes_TIME),residuals(mod_cytes_TIME))
qqnorm(residuals(mod_cytes_TIME)) # qqplot

## zoa test---------------------------------
mod_zoa_TIME  <- aov(mean_zoa~Treatment*Date, data = Means_Table)
anova(mod_zoa_TIME) # p = 0.01051 * (time)
shapiro.test(residuals(mod_zoa_TIME)) #  normal residuals p-value = 0.08541
leveneTest(mod_zoa_TIME) # p = 0.1719
hist(residuals(mod_zoa_TIME)) #plot histogram of residuals
boxplot(residuals(mod_zoa_TIME)) #plot boxplot of residuals
plot(fitted(mod_zoa_TIME),residuals(mod_zoa_TIME))
qqnorm(residuals(mod_zoa_TIME)) # qqplot

TukeyHSD(mod_zoa_TIME)

## lumen test---------------------------------
mod_lumen_TIME  <- aov(mean_lumen~Treatment*Date, data = Means_Table)
anova(mod_lumen_TIME) # p = 0.06423 . (treat) ; 0.09255 . (time)
shapiro.test(residuals(mod_lumen_TIME)) #  normal residuals p-value = 0.2685
leveneTest(mod_lumen_TIME) # p = 0.1227
hist(residuals(mod_lumen_TIME)) #plot histogram of residuals
boxplot(residuals(mod_lumen_TIME)) #plot boxplot of residuals
plot(fitted(mod_lumen_TIME),residuals(mod_lumen_TIME))
qqnorm(residuals(mod_lumen_TIME)) # qqplot

TukeyHSD(mod_lumen_TIME)

# One-Way anova
## cytes test---------------------------------
mod_cytes  <- aov(mean_lumen~Treatment*Date, data = Means_Table)
anova(mod_cytes) # p = 0.9391
shapiro.test(residuals(mod_cytes)) #  normal residuals p-value = 0.3263
leveneTest(mod_cytes) # p = 0.6256
hist(residuals(mod_cytes)) #plot histogram of residuals
boxplot(residuals(mod_cytes)) #plot boxplot of residuals
plot(fitted(mod_cytes),residuals(mod_cytes))
qqnorm(residuals(mod_zoa)) # qqplot


### zoa test---------------------------------
anova(mod_zoa) # p = 0.494
shapiro.test(residuals(mod_zoa)) #  normal residuals p-value = 0.2042
leveneTest(mod_zoa) # p = 0.5209
hist(residuals(mod_zoa)) #plot histogram of residuals
boxplot(residuals(mod_zoa)) #plot boxplot of residuals
plot(fitted(mod_zoa),residuals(mod_zoa))
qqnorm(residuals(mod_zoa)) # qqplot

## cytes test---------------------------------
mod_cytes  <- aov(mean_perc_cytes~Treatment, data = Means_Table)
anova(mod_cytes) # p = 0.9391
shapiro.test(residuals(mod_cytes)) #  normal residuals p-value = 0.3263
leveneTest(mod_cytes) # p = 0.6256
hist(residuals(mod_cytes)) #plot histogram of residuals
boxplot(residuals(mod_cytes)) #plot boxplot of residuals
plot(fitted(mod_cytes),residuals(mod_cytes))
qqnorm(residuals(mod_zoa)) # qqplot

## lumen test---------------------------------
mod_lumen  <- aov(mean_perc_lumen~Treatment, data = Means_Table)
anova(mod_lumen) # p = 0.1876
shapiro.test(residuals(mod_lumen)) #  normal residuals p-value = 0.4007
leveneTest(mod_lumen) # p = 0.6256
hist(residuals(mod_lumen)) #plot histogram of residuals
boxplot(residuals(mod_lumen)) #plot boxplot of residuals
plot(fitted(mod_lumen),residuals(mod_lumen))
qqnorm(residuals(mod_lumen)) # qqplot

## zoa : cytes test---------------------------
mod_zoa_cytes <- aov(mean_zoa_cyte_ratio~Treatment, data = Means_Table)
anova(mod_zoa_cytes) # p = 0.8147
shapiro.test(residuals(mod_zoa_cytes)) #  normal residuals p-value = 0.4334
leveneTest(mod_zoa_cytes) # p = 0.4916
hist(residuals(mod_zoa_cytes)) #plot histogram of residuals
boxplot(residuals(mod_zoa_cytes)) #plot boxplot of residuals
plot(fitted(mod_zoa_cytes),residuals(mod_zoa_cytes))
qqnorm(residuals(mod_zoa_cytes)) # qqplot
