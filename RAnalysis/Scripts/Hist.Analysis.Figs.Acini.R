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
            mean_total_acini_area = mean(total_area, na.rm = TRUE),
            num_acini=n())
Means_Table # view the table

# PIVOT LONGER-----------------------------------------------------------------------------------------
Means_Table_long <- Means_Table %>% 
  tidyr::pivot_longer(cols = c(4:7), names_to='scoring_metric', values_to='means')

#  PLOTS----------------------------------------------------------------------------------------
# ocrete variable Date_Treat
Means_Table_long$Date_Treat <- paste((substr(Means_Table_long$Date,5,8)),Means_Table_long$Treatment, sep ="_")
list(Means_Table_long$Date_Treat)
# order the levels of the factor date_treatment to order correly in the plot
Means_Table_long$Date_Treat <- factor(Means_Table_long$Date_Treat,levels = c("0123_Ambient", "0123_Low", "0221_Ambient", "0221_Low"))



# plots
plot <- ggboxplot(Means_Table_long, x = "Date_Treat", y = "means",  fill = "Treatment",
                     palette = c("#00AFBB", "#FC4E07"), add = "none")
plot2 <- plot %>% ggadd(shape ="Treatment",fill = "white") %>% ggadd("jitter", size = 4,shape ="Treatment",fill = "white") +
  facet_wrap( ~ scoring_metric, ncol=2, scales = "free") + theme_classic()
plot2

# STATISTICS------------------------------------------------------------------------------------------
Means_Table$Treatment <- as.factor(Means_Table$Treatment)
Means_Table$Date <- as.factor(Means_Table$Date)
par(mfrow=c(1,4)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots

# Two-Way ANOVA (treat×time)
typeof(Means_Table$Date) # currenty an integer - must change to a character
Means_Table$Date <- as.character(Means_Table$Date)
## total acini area test---------------------------------
mod_total_acini_TIME  <- aov(mean_total_acini_area~Treatment*Date, data = Means_Table)
anova(mod_total_acini_TIME) # p = 0.9391
shapiro.test(residuals(mod_total_acini_TIME)) #  normal residuals p-value = 7.261e-06
leveneTest(mod_total_acini_TIME) # p = 0.6463
hist(residuals(mod_total_acini_TIME)) #plot histogram of residuals
boxplot(residuals(mod_total_acini_TIME)) #plot boxplot of residuals
plot(fitted(mod_total_acini_TIME),residuals(mod_total_acini_TIME))
qqnorm(residuals(mod_total_acini_TIME)) # qqplot

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