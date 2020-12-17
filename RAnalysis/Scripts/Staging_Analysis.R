#Title: Pgenerosa histology staging analysis
#Project: FFAR
#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200716
#See Readme file for details

rm(list=ls())

# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(tidyr)
library(ggpubr)
library(viridis)
library(hrbrthemes)

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_histology/RAnalysis/") #set working

# upload data
data<-read.csv("Data/Staging/Hist_Staging_Kaitlyn.csv", header=T, sep=",", na.string="NA", as.is=T) 
data # view data

# RE: 'Stage' and 'Staging_number' columns
# Immature: very early active (1),
# early active (2)) 
# Mature: late active (3)
# ripe (4))
# Spent (5)

data_pull <- data %>% dplyr::select(c('Tank','Treatment','Date', 'Sex', 'Geoduck_ID', 'Staging_number', 'Stage_ID'))
data_pull$Staging_number <- as.character(data_pull$Staging_number)
# calcaulate the proportions of staging number by date and treatment 
prop_ALL <-
  data_pull %>% 
  group_by(Date, Treatment, Sex, Stage_ID) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
prop_ALL$prop <- prop_ALL$freq*100

prop_Jan23_Female <- prop_ALL %>% dplyr::filter('123' %in% Date) %>%  dplyr::filter('F' %in% Sex)
prop_Jan23_Male <- prop_ALL %>% dplyr::filter('123' %in% Date) %>%  dplyr::filter('M' %in% Sex)
prop_Feb21_Female <- prop_ALL %>% dplyr::filter('221' %in% Date) %>%  dplyr::filter('F' %in% Sex)
prop_Feb21_Male <- prop_ALL %>% dplyr::filter('221' %in% Date) %>%  dplyr::filter('M' %in% Sex)

### separate datasets for six, date, and treatment 
# female
prop_Jan23_Female_Amb <- prop_Jan23_Female %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Jan23_Female_Low <- prop_Jan23_Female %>%
  dplyr::filter('Low' %in% Treatment)

prop_Feb21_Female_Amb <- prop_Feb21_Female %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Feb21_Female_Low <- prop_Feb21_Female %>%
  dplyr::filter('Low' %in% Treatment)
# Male
prop_Jan23_Male_Amb <- prop_Jan23_Male %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Jan23_Male_Low <- prop_Jan23_Male %>%
  dplyr::filter('Low' %in% Treatment)

prop_Feb21_Male_Amb <- prop_Feb21_Male %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Feb21_Male_Low <- prop_Feb21_Male %>%
  dplyr::filter('Low' %in% Treatment)


# STACKED BAR PLOT 

# Stacked
prop_ALL$Date_Treatment <- paste(prop_ALL$Date, prop_ALL$Treatment, sep ="_") # new column merges date and treat for stacked bar plots
prop_female <- prop_ALL %>% dplyr::filter(Sex %in% 'F') # new dataset to plot JUST female histology
prop_male <- prop_ALL %>% dplyr::filter(Sex %in% 'M')  # new dataset to plot JUST male histology

StackedBar_FEMALE <- ggplot(prop_female, aes(fill=Stage_ID, y= n, x=Date_Treatment)) + 
                      geom_bar(position="stack", stat="identity") +
                      scale_fill_viridis(discrete = T) +
                      ggtitle("Female histology staging") +
                      theme_ipsum() +
                      xlab("Date_Treatment") +
                      ylab("number of samples")
StackedBar_FEMALE# view plot

StackedBar_MALE <- ggplot(prop_male, aes(fill=Stage_ID, y= n, x=Date_Treatment)) + 
                      geom_bar(position="stack", stat="identity") +
                      scale_fill_viridis(discrete = T) +
                      ggtitle("Male histology staging") +
                      theme_ipsum() +
                      xlab("Date_Treatment") +
                      ylab("number of samples")
StackedBar_MALE# view plot

# grid plots
StackedBar_Plots <- grid.arrange(StackedBar_FEMALE, StackedBar_MALE, ncol =2, nrow = 1)
StackedBar_Plots # view plot
#  SAVE
ggsave(file="Stacked_barplots_Staging.pdf", StackedBar_Plots, width = 24, height = 12, units = c("in")) 


# DONUT PLOTS
F_amb_123 <- ggdonutchart(prop_Jan23_Female_Amb, "prop", label = "Stage_ID",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Ambient_0123", 
             fill = "Stage_ID", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00"))

F_low_123 <- ggdonutchart(prop_Jan23_Female_Low, "prop", label = "Stage_ID",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Low_0123", 
             fill = "Stage_ID", color = "white", palette = c("#56B4E9", "#0072B2"))
             
F_amb_221 <- ggdonutchart(prop_Feb21_Female_Amb, "prop", label = "Stage_ID",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Ambient_0221", 
             fill = "Stage_ID", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))

F_low_221 <- ggdonutchart(prop_Feb21_Female_Low, "prop", label = "Stage_ID",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Low_0221", 
             fill = "Stage_ID", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00"))


M_amb_123 <- ggdonutchart(prop_Jan23_Male_Amb, "prop", label = "Stage_ID",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Ambient_0123", 
                          fill = "Stage_ID", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))

M_low_123 <- ggdonutchart(prop_Jan23_Male_Low, "prop", label = "Stage_ID",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Low_0123", 
                          fill = "Stage_ID", color = "white", palette = c("#56B4E9","#CC79A7"))

M_amb_221 <- ggdonutchart(prop_Feb21_Male_Amb, "prop", label = "Stage_ID",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Ambient_0221",                
                          fill = "Stage_ID", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))

M_low_221 <- ggdonutchart(prop_Feb21_Male_Low, "prop", label = "Stage_ID",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Low_0221", 
                          fill = "Stage_ID", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))


ggarrange(M_amb_123, M_low_123,F_amb_123, F_low_123,
          M_amb_221, M_low_221, F_amb_221, F_low_221,  ncol = 4, nrow = 2)
