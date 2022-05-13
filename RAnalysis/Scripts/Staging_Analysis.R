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
library(gridExtra)
library(rstatix)
#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/Github_repositories/paper-GeoduckReproDev_pH/data/histology/RAnalysis") #set working

# upload data
histology<-read.csv("Data/Staging/Hist_Staging_Trigg_Crandall.csv", header=T, sep=",", na.string="NA", as.is=T) 
histology # view data
histology$Stage <- ifelse(is.na(histology$Stage_Crandall), histology$Stage_Trigg, histology$Stage_Crandall) # call Grace's Stages (Stage_Crandall) - however when NA (photos analyzed by Shelly only) replace with Stage_Trigg

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







prop_Jan23_Female # view for the row numbers in Low and Ambient
Fem_Jan23 <-  as.table(rbind(c(0,60,40), c(66.7, 33.3, 0)))
dimnames(Fem_Jan23) <- list(
  group = c("Ambient", "Low"),
  stage = c("stage5", "stage6", "stage7") )
prop_test(Fem_Jan23) # significantly different! 




prop_Feb21_Female # view for the row numbers in Low and Ambient
Fem_Feb21 <- as.table(rbind(c(prop_Feb21_Female$prop[c(1:2)]), c(prop_Feb21_Female$prop[c(3:4)])))
dimnames(Fem_Feb21) <- list(
  group = c("Low", "Ambient"),
  stage = c("stage6", "stage7") )
prop_test(Fem_Feb21) # not significantly different!


res <- prop.test(x = c(6, 2), 
                 n = c(8, 3))
res # p-value = 1

res <- prop.test(x = c(2, 1), 
                 n = c(8, 3))
res # p-value = 1



prop_Jan23_Male
Male_Jan23 <- as.table(rbind(c(14.3,85.7,0), c(40,0,60)))
dimnames(Male_Jan23) <- list(
  group = c("Low", "Ambient"),
  stage = c("stage1_2", "stage3_4","stage5")
)
prop_test(Male_Jan23) # significant

res <- prop.test(x = c(6, 3), 
                 n = c(7, 5))
res # p-value = 0.7353 # stages 3 and 4 






prop_Feb21_Male
Male_Feb21 <- as.table(rbind(c(0,0,33.3, 66.7), c(22.2,11.1,44.4,22.2)))
dimnames(Male_Feb21) <- list(
  group = c("Ambient", "Low"),
  stage = c("stage2", "stage3","stage4", "stage5")
)
prop_test(Male_Feb21) # significant!


res <- prop.test(x = c(1, 4), 
                 n = c(3, 9))
res # p-value = 0.7353 # stages 3 and 4 
















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
