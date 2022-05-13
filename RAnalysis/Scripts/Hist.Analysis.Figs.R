#Title: Pgenerosa histology
#Project: FFAR
#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200716
#See Readme file for details

rm(list=ls())

install.packages('ggpubr')
# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(tidyr)
library(ggpubr)

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_histology/RAnalysis") #set working

# upload data
dat<-read.csv("Data/20190123_male_scoring.csv", header=T, sep=",", na.string="NA", as.is=T) 
dat # view data

#  TABLE 
table <- dat %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(VALUE = area_.um.3.[target_characteristic == "Spematozoa+Spermatocytes"] / area_.um.3.[target_characteristic == "Total Area"]) %>% 
  select(c("Treatment", "ID","target_characteristic", "area_.um.3.", "VALUE")) %>% 
  dplyr::filter(target_characteristic %in% ("Spematozoa+Spermatocytes"))

#  PLOT
plot <- ggviolin(table, x = "Treatment", y = "VALUE",  ylab = "% area spermatozoa and spermatocytes (µm^3)",  fill = "Treatment",
                          palette = c("#00AFBB","#FC4E07"), add = "none", title = "Gonad histology (Male; samples from 20190123)")
plot_box_violin <- plot %>% ggadd("boxplot",shape ="Treatment",fill = "white") %>% ggadd("jitter", size = 4,shape ="Treatment",fill = "white")
plot_box_violin

#  STAT
table$Treatment.num <- as.factor(table$Treatment)
table$Treatment.num <- as.numeric(table$Treatment.num)
t.test(table$VALUE~table$Treatment.num) # not different! t = 0.47357, df = 6.4737, p-value = 0.6514

#  SAVE
ggsave(file="Output/Acini_area_plot.pdf", plot_box_violin, width = 12, height = 8, units = c("in")) # save plot
