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
library(gridExtra)
#set working directory---------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/Github_repositories/Pgenerosa_histology/RAnalysis") #set working


######################################################################### #
############################## DATA PREP ################################ #
######################################################################### #

# UPLOAD DATA------------------------------------------------------------------------------------------
staging<-read.csv("Data/Staging/Hist_Staging_Kaitlyn.csv", header=T, sep=",", na.string="NA", as.is=T) 
staging
Male_hist<-read.csv("Data/Master_summary_male_acini.csv", header=T, sep=",", na.string="NA", as.is=T) 
Male_hist # view data
Male_hist$Area = as.numeric(Male_hist$Area)
# PIVOT THE TABLE TO CALC RELATIVE VALUES--------------------------------------------------------------
Male_hist_2 <- Male_hist %>% dplyr::select(-c('Meas_num','Label','hue','saturation','brightness')) %>% 
  tidyr::pivot_wider(names_from=type, values_from=Area)
Male_hist_2 <- Male_hist_2 %>% dplyr::select(-'total_area')
# CALC RELATIVE VALUES---------------------------------------------------------------------------------
Male_hist_2$TotalArea <- (Male_hist_2$lumen + Male_hist_2$zoa + (Male_hist_2$cytes_zoa - Male_hist_2$zoa))
Male_hist_2$perc_zoa <- (Male_hist_2$zoa/Male_hist_2$TotalArea)*100 # percent area of spermatozoa
Male_hist_2$perc_cytes <- ((Male_hist_2$cytes_zoa - Male_hist_2$zoa)/Male_hist_2$TotalArea)*100 # percent area of spermatocytes
Male_hist_2$perc_lumen <- (Male_hist_2$lumen/Male_hist_2$TotalArea)*100 # percent area of lumen
Male_hist_2$cytes <- (Male_hist_2$cytes_zoa - Male_hist_2$zoa) # area of spermatocytes
# convert to characer
typeof(Male_hist_2$ID) 
typeof(Male_hist_2$Date) 
Male_hist_2$Date <- as.character(Male_hist_2$Date)


######################################################################### #
##################      TWO-WAY ANOVA     ############################### #
######################################################################### #
# ZOA TEST
Zoa_2factorialAOV <- aov(zoa ~ Treatment * Date, data = Male_hist_2)
summary(Zoa_2factorialAOV) # date is a sig effect in raw data (not transformed!)
TukeyHSD(Zoa_2factorialAOV, which = "Date") # 02/21 > 01/23; increased with time
plot(Zoa_2factorialAOV, 1) # TEST ASSUMPTIONS     # 1. Homogeneity of variances (Levenes test) of the model data
leveneTest(zoa ~ Treatment * Date, data = Male_hist_2) # DOES NOT PASS homogeneity of variance; leveneTest is from the 'car' package
plot(Zoa_2factorialAOV, 2)  # TEST ASSUMPTIONS    # 2. Normality - QQ plot (quantile quantile) of the model residuals
Zoa_aov_residuals <- residuals(object = Zoa_2factorialAOV) # Extract the residuals
shapiro.test(x = Zoa_aov_residuals) # Run Shapiro-Wilk test; NOT NORMAL
hist(Male_hist_2$zoa)
####  ZOA TEST transformed #### #
Zoa_2factorialAOV.transform <- aov(log(zoa+1) ~ Treatment * Date, data = Male_hist_2)
summary(Zoa_2factorialAOV.transform) # still effect of Date in transformation
TukeyHSD(Zoa_2factorialAOV.transform, which = "Date") # 02/21 > 01/23; increased with time
plot(Zoa_2factorialAOV.transform, 1) # TEST ASSUMPTIONS     # 1. Homogeneity of variances (Levenes test) of the model data
leveneTest(sqrt(zoa) ~ Treatment * Date, data = Male_hist_2) # DOES NOT PASS homogeneity of variance; leveneTest is from the 'car' package
plot(Zoa_2factorialAOV.transform, 2)  # TEST ASSUMPTIONS    # 2. Normality - QQ plot (quantile quantile) of the model residuals
Zoa_aov_residuals.transform <- residuals(object = Zoa_2factorialAOV.transform) # Extract the residuals
shapiro.test(x = Zoa_aov_residuals.transform) # Run Shapiro-Wilk test; NORMAL
hist(log(Male_hist_2$zoa)+1)
# LUMEN TEST
lumen_2factorialAOV <- aov(lumen ~ Treatment * Date, data = Male_hist_2)
summary(lumen_2factorialAOV) # treatment and data are sig diff in raw data (not transformed!)
TukeyHSD(lumen_2factorialAOV)
TukeyHSD(lumen_2factorialAOV, which = "Treatment") # Low < Ambient
TukeyHSD(lumen_2factorialAOV, which = "Date") # 02/21 > 01/23; increased with time
plot(lumen_2factorialAOV, 1) # TEST ASSUMPTIONS     # 1. Homogeneity of variances (Levenes test) of the model data
leveneTest(lumen ~ Treatment * Date, data = Male_hist_2) # DOES NOT PASS homogeneity of varianceE; leveneTest is from the 'car' package
plot(lumen_2factorialAOV, 2)  # TEST ASSUMPTIONS    # 2. Normality - QQ plot (quantile quantile) of the model residuals
lumen_aov_residuals <- residuals(object = lumen_2factorialAOV) # Extract the residuals
shapiro.test(x = lumen_aov_residuals) # Run Shapiro-Wilk test; NOT NORMAL
hist(Male_hist_2$lumen)
####  lumen TEST transformed #### #
lumen_2factorialAOV.transform <- aov(sqrt(lumen) ~ Treatment * Date, data = Male_hist_2)
summary(lumen_2factorialAOV.transform) # same effects as raw data; no marginal interaction
TukeyHSD(lumen_2factorialAOV.transform, which = "Treatment") # Low < Ambient
TukeyHSD(lumen_2factorialAOV.transform, which = "Date") # 02/21 > 01/23; increased with time
plot(lumen_2factorialAOV.transform, 1) # TEST ASSUMPTIONS     # 1. Homogeneity of variances (Levenes test) of the model data
leveneTest(sqrt(lumen) ~ Treatment * Date, data = Male_hist_2) # DOES NOT PASS homogeneity of varianceE; leveneTest is from the 'car' package
plot(lumen_2factorialAOV.transform, 2)  # TEST ASSUMPTIONS    # 2. Normality - QQ plot (quantile quantile) of the model residuals
lumen_aov_residuals.transform <- residuals(object = lumen_2factorialAOV.transform) # Extract the residuals
shapiro.test(x = lumen_aov_residuals.transform) # Run Shapiro-Wilk test; NORMAL
hist(sqrt(Male_hist_2$lumen))
hist(log(Male_hist_2$lumen))
# CYTES TEST
cytes_2factorialAOV <- aov(cytes ~ Treatment * Date, data = Male_hist_2)
summary(cytes_2factorialAOV) # Date has a marginal diff (not transformed!)
TukeyHSD(cytes_2factorialAOV, which = "Date") # 02/21 > 01/23; increased with time
plot(cytes_2factorialAOV, 1) # TEST ASSUMPTIONS     # 1. Homogeneity of variances (Levenes test) of the model data
leveneTest(cytes ~ Treatment * Date, data = Male_hist_2) # PASSED homogeneity of variance; 0.3056
plot(cytes_2factorialAOV, 2)  # TEST ASSUMPTIONS    # 2. Normality - QQ plot (quantile quantile) of the model residuals
cytes_aov_residuals <- residuals(object = cytes_2factorialAOV) # Extract the residuals
shapiro.test(x = cytes_aov_residuals) # Run Shapiro-Wilk test; NOT NORMAL
####  cytes TEST transformed #### #
cytes_2factorialAOV.transform <- aov(sqrt(cytes) ~ Treatment * Date, data = Male_hist_2)
summary(cytes_2factorialAOV.transform) # additional effect of treatment!
TukeyHSD(cytes_2factorialAOV.transform, which = "Treatment") # Low < Ambient
plot(cytes_2factorialAOV.transform, 1) # TEST ASSUMPTIONS     # 1. Homogeneity of variances (Levenes test) of the model data
leveneTest(sqrt(cytes) ~ Treatment * Date, data = Male_hist_2) # PASS homogeneity of variance; 0.2546
plot(cytes_2factorialAOV.transform, 2)  # TEST ASSUMPTIONS    # 2. Normality - QQ plot (quantile quantile) of the model residuals
cytes_aov_residuals.transform <- residuals(object = cytes_2factorialAOV.transform) # Extract the residuals
shapiro.test(x = cytes_aov_residuals.transform) # Run Shapiro-Wilk test; NON NORMAL - 0.00325






######################################################################### #
##################      TWO-WAY ANOVA PROPORTION     ############################### #
######################################################################### #
par(mfrow=c(2,2)) # 2 by 2 grid for plotting
# PERC ZOA TEST
PERCZoa_2factorialAOV <- aov(perc_zoa ~ Treatment * Date, data = Male_hist_2)
summary(PERCZoa_2factorialAOV) # date is a sig effect in raw data (not transformed!)
TukeyHSD(PERCZoa_2factorialAOV, which = "Date") # 02/21 > 01/23; increased with time
plot(PERCZoa_2factorialAOV) # residuals look normal  - no influenctial outliers 
leveneTest(perc_zoa ~ Treatment * Date, data = Male_hist_2) # PASS
PERCZoa_aov_residuals <- residuals(object = PERCZoa_2factorialAOV) # Extract the residuals
shapiro.test(x = PERCZoa_aov_residuals) # Run Shapiro-Wilk test; PASS
hist(Male_hist_2$perc_zoa) # right (positive) skew

# PERC LUMEN TEST
PERClumen_2factorialAOV <- aov(perc_lumen ~ Treatment * Date, data = Male_hist_2)
summary(PERClumen_2factorialAOV) # effect of treatment
TukeyHSD(PERClumen_2factorialAOV, which = "Treatment")  # Low < Ambient -3.451896 % area
plot(PERClumen_2factorialAOV) # residuals look normal  - no influenctial outliers 
leveneTest(perc_lumen ~ Treatment * Date, data = Male_hist_2) # PASSES
PERClumen_aov_residuals <- residuals(object = PERClumen_2factorialAOV) # Extract the residuals
shapiro.test(x = PERClumen_aov_residuals) # Run Shapiro-Wilk test; NOT NORMAL
hist(Male_hist_2$perc_lumen) # right (positive) skew

# PERC CYTES TEST
PERCcytes_2factorialAOV <- aov(perc_cytes ~ Treatment * Date, data = Male_hist_2)
summary(PERCcytes_2factorialAOV) # Date has a significant effect
TukeyHSD(PERCcytes_2factorialAOV, which = "Date") # 02/21 < 01/23; cytes DECREASED with time -11.29722%
plot(PERCcytes_2factorialAOV) # residuals look normal  - no influenctial outliers 
leveneTest(perc_cytes ~ Treatment * Date, data = Male_hist_2) # PASSED homogeneity of variance
PERCcytes_aov_residuals <- residuals(object = PERCcytes_2factorialAOV) # Extract the residuals
shapiro.test(x = PERCcytes_aov_residuals) # Run Shapiro-Wilk test;  PASS, normal
hist(Male_hist_2$perc_cytes) # right (positive) skew

######################################################################### #
##################      BOX PLOT MODEL    ############################### #
######################################################################### #
Male_hist_2$Treatment <- ordered(Male_hist_2$Treatment, levels = c("Ambient", "Low"))
cb_Colors <- c("#4E84C4", "#D16103")
# Zoa
zoa_plot <- ggplot(Male_hist_2, aes(x = Date, y = zoa, colour = Treatment)) +
  theme_bw() +
  geom_boxplot(outlier.size = 1, position=position_dodge(0.9), fill = "white") +
  geom_point(pch = 21, size = 3, fill = "white", position = position_jitterdodge(0.2)) +
  scale_color_manual(values=cb_Colors)
zoa_plot # view the plot
zoa_plot_transformed <- ggplot(Male_hist_2, aes(x = Date, y = sqrt(zoa), colour= Treatment)) +
  theme_bw() +
  geom_boxplot(outlier.size = 1, position=position_dodge(0.9), fill = "white") +
  geom_point(pch = 21, size = 3, fill = "white", position = position_jitterdodge(0.2)) +
  scale_color_manual(values=cb_Colors)
zoa_plot_transformed # view the plot
# lumen
lumen_plot <- ggplot(Male_hist_2, aes(x = Date, y = lumen, colour = Treatment)) +
  theme_bw() +
  geom_boxplot(outlier.size = 1, position=position_dodge(0.9), fill = "white") +
  geom_point(pch = 21, size = 3, fill = "white", position = position_jitterdodge(0.2)) +
  scale_color_manual(values=cb_Colors)
lumen_plot # view the plot
lumen_plot_transformed <- ggplot(Male_hist_2, aes(x = Date, y = sqrt(lumen), colour= Treatment)) +
  theme_bw() +
  geom_boxplot(outlier.size = 1, position=position_dodge(0.9), fill = "white") +
  geom_point(pch = 21, size = 3, fill = "white", position = position_jitterdodge(0.2)) +
  scale_color_manual(values=cb_Colors)
lumen_plot_transformed # view the plot
# cytes
cytes_plot <- ggplot(Male_hist_2, aes(x = Date, y = cytes, colour = Treatment)) +
  theme_bw() +
  geom_boxplot(outlier.size = 1, position=position_dodge(0.9), fill = "white") +
  geom_point(pch = 21, size = 3, fill = "white", position = position_jitterdodge(0.2)) +
  scale_color_manual(values=cb_Colors)
cytes_plot # view the plot
cytes_plot_transformed <- ggplot(Male_hist_2, aes(x = Date, y = log(cytes), colour= Treatment)) +
  theme_bw() +
  geom_boxplot(outlier.size = 1, position=position_dodge(0.9), fill = "white") +
  geom_point(pch = 21, size = 3, fill = "white", position = position_jitterdodge(0.2)) +
  scale_color_manual(values=cb_Colors)
cytes_plot_transformed # view the plot

# grid plots
ALL_PLOTS <- grid.arrange(zoa_plot, lumen_plot, cytes_plot,
                          zoa_plot_transformed, lumen_plot_transformed, cytes_plot_transformed, ncol =3, nrow = 2)
ALL_PLOTS # view plot
#  SAVE
ggsave(file="Output/Grid_plot.pdf", ALL_PLOTS, width = 12, height = 8, units = c("in")) 


######################################################################### #
###LINEAR REG WITH STAGE AS INDEP VAR  ############################### #
######################################################################### #
Male_hist_2$Geoduck_ID <- Male_hist_2$ID
staging_coltrim        <- staging %>% dplyr::select(c('Geoduck_ID','Geoduck_ID','Stage_ID','Staging_number'))
MaleHistStage_merge    <- merge(Male_hist_2,staging_coltrim,by="Geoduck_ID") %>% 
                                dplyr::mutate(ID = gsub('.*_00', '', Geoduck_ID)) %>% 
                                dplyr::mutate(ID = gsub('.*_0', '', Geoduck_ID))

MaleHistStage_221     <- MaleHistStage_merge %>% dplyr::filter(Date %in% '20190221')
MaleHistStage_123     <- MaleHistStage_merge %>% dplyr::filter(Date %in% '20190123')

MaleHistStage_low     <- MaleHistStage_merge %>% dplyr::filter(Treatment %in% 'Low')
MaleHistStage_ambient <- MaleHistStage_merge %>% dplyr::filter(Treatment %in% 'Ambient')

PERClumen_staging <- ggplot(MaleHistStage_merge, aes(factor(Staging_number), perc_lumen)) +
                        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
                        # geom_point(size = 3, position = position_jitterdodge(jitter.width = 0.1)) +
                        #geom_jitter(position = position_jitter()) +
                        geom_text(aes(label=ID, color = Treatment), position = position_jitter(height = .25, width = .25), size = 4) +
                         scale_color_manual(values = c("#00BFC4","#F8766D")) +
                        # scale_color_manual(values = c("#0072B2", "#D55E00")) + # colorblindness palette blue and orange
                        theme_classic() +
                        labs(y=expression("Lumen (% area total acini)"), x=expression("Staging ID")) +
                        geom_smooth(method = "lm", se=T, color="grey25", alpha = 0.2, aes(group=1)) +
                        theme(legend.position = "none")
PERClumen_staging2 <-PERClumen_staging + theme(text = element_text(size = 15))# view plot
PERClumen_staging2 # ambient is blue and elevated pCO2 is orange


PERCzoa_staging <- ggplot(MaleHistStage_merge, aes(as.factor(Staging_number), perc_zoa)) +
                    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
                    # geom_point(size = 3, position = position_jitterdodge(jitter.width = 0.1)) +
                    #geom_jitter(position = position_jitter()) +
                    geom_text(aes(label=ID, color = Treatment), position = position_jitter(height = .25, width = .25), size = 4) +
                     scale_color_manual(values = c("#00BFC4","#F8766D")) +
                    #scale_color_manual(values = c("#0072B2", "#D55E00")) + # colorblindness palette blue and orange
                    theme_classic() +
                    labs(y=expression("Spermatozoa (% area total acini)"), x=expression("Staging ID")) +
                    geom_smooth(method = "lm", se=T, color="grey25", alpha = 0.2, aes(group=1)) +
                    theme(legend.position = "none")
                  PERCzoa_staging2 <- PERCzoa_staging + theme(text = element_text(size = 15))# view plot
PERCzoa_staging2



PERCcytes_staging <- ggplot(MaleHistStage_merge, aes(factor(Staging_number), perc_cytes)) +
                      geom_boxplot(alpha = 0.5, outlier.shape = NA) +
                      # geom_point(size = 3, position = position_jitterdodge(jitter.width = 0.1)) +
                      #geom_jitter(position = position_jitter()) +
                      geom_text(aes(label=ID, color = Treatment), position = position_jitter(height = .25, width = .25), size = 4) +
                       scale_color_manual(values = c("#00BFC4","#F8766D")) +
                      #scale_color_manual(values = c("#0072B2", "#D55E00")) + # colorblindness palette blue and orange
                      theme_classic() +
                      labs(y=expression("Spermatocytes (% area total acini)"), x=expression("Staging ID")) +
                      geom_smooth(method = "lm", se=T, color="grey25", alpha = 0.2, aes(group=1)) +
                      theme(legend.position = "none")
PERCcytes_staging2 <-PERCcytes_staging + theme(text = element_text(size = 15))  # view plot
PERCcytes_staging2


# grid plots
staging_hist_plots <- grid.arrange(PERClumen_staging2, PERCzoa_staging2, PERCcytes_staging2, ncol =3, nrow = 1)
staging_hist_plots # view plot
#  SAVE
ggsave(file="Output/StagingHist_regression_plot.pdf", staging_hist_plots, width = 14, height = 7, units = c("in")) 

### same plots but for Low and Ambient treatment

PERClumen_staging_LOW <- ggplot(MaleHistStage_low, aes(factor(Staging_number), perc_lumen)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Lumen_Low") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))
PERClumen_staging_LOW2 <-PERClumen_staging_LOW + theme(text = element_text(size = 20))# view plot

PERClumen_staging_AMBIENT <- ggplot(MaleHistStage_ambient, aes(factor(Staging_number), perc_lumen)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Lumen_Ambient") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))
PERClumen_staging_AMBIENT2 <-PERClumen_staging_AMBIENT + theme(text = element_text(size = 20))# view plot

PERCzoa_staging_LOW <- ggplot(MaleHistStage_low, aes(factor(Staging_number), perc_zoa)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Zoa_Low") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))
PERCzoa_staging_LOW2 <- PERCzoa_staging_LOW + theme(text = element_text(size = 20))# view plot

PERCzoa_staging_AMBIENT <- ggplot(MaleHistStage_ambient, aes(factor(Staging_number), perc_zoa)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Zoa_Ambient") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))
PERCzoa_staging_AMBIENT2 <- PERCzoa_staging_AMBIENT + theme(text = element_text(size = 20))# view plot

PERCcytes_staging_LOW <- ggplot(MaleHistStage_low, aes(factor(Staging_number), perc_cytes)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Cytes_Low") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))
PERCcytes_staging_LOW2 <-PERCcytes_staging_LOW + theme(text = element_text(size = 20))  # view plot

PERCcytes_staging_AMBIENT <- ggplot(MaleHistStage_ambient, aes(factor(Staging_number), perc_cytes)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Cytes_Ambient") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))
PERCcytes_staging_AMBIENT2 <-PERCcytes_staging_AMBIENT + theme(text = element_text(size = 20))  # view plot

# grid plots
staging_hist_plots_TREATMENTS <- grid.arrange(PERClumen_staging_AMBIENT2, PERCzoa_staging_AMBIENT2, PERCcytes_staging_AMBIENT2, 
                                              PERClumen_staging_LOW2, PERCzoa_staging_LOW2, PERCcytes_staging_LOW2, ncol =3, nrow = 2)
staging_hist_plots_TREATMENTS # view plot
#  SAVE
ggsave(file="Output/StagingHist_regressionTREATMENT_plot.pdf", staging_hist_plots_TREATMENTS, width = 12, height = 8, units = c("in")) 


######################################################################### #
################## PROPORTIONS BAR CHART  ############################### #
######################################################################### #

# call subset for plotting and pivot longer with tidyr
colnames(Male_hist_2)
Male_hist_plots <- Male_hist_2[,c(1:3,9:11)]
Male_hist_plots_long <- Male_hist_plots %>% 
  tidyr::pivot_longer(cols = c(4:6), names_to='prop_metric', values_to='value')
# change labels for plotting with gsub
Male_hist_plots_long$Date <- gsub("20190123", "72",Male_hist_plots_long$Date)
Male_hist_plots_long$Date <- gsub("20190221", "93 + 8 day recovery",Male_hist_plots_long$Date)
Male_hist_plots_long$prop_metric <- gsub("perc_cytes", "spermatocytes",Male_hist_plots_long$prop_metric)
Male_hist_plots_long$prop_metric <- gsub("perc_zoa", "spermatozoa",Male_hist_plots_long$prop_metric)
Male_hist_plots_long$prop_metric <- gsub("perc_lumen", "lumen",Male_hist_plots_long$prop_metric)
#convert metric to ordered factor
Male_hist_plots_long_ordered <- Male_hist_plots_long[order(Male_hist_plots_long$prop_metric ),]
Male_hist_plots_long_ordered$prop_metric <- factor(Male_hist_plots_long_ordered$prop_metric, levels = rev(c(rep("lumen",1), rep("spermatocytes",1), rep("spermatozoa",1))))
Male_hist_plots_long_ordered_2 <- Male_hist_plots_long_ordered%>%
  group_by(Date, Treatment, prop_metric) %>% 
  dplyr::summarise(mean = mean(value), sd =sd(value)) %>% 
  mutate(y_pos = cumsum(mean))
# stacked proportion plot
pd <- position_dodge2(width = 0.2)
Male_hist_plots_long_ordered_2$prop_metric <- as.character(Male_hist_plots_long_ordered_2$prop_metric)

proportion_Figure <- ggplot(Male_hist_plots_long_ordered_2, aes(x = Treatment, y = mean, fill = Treatment, alpha = prop_metric)) + 
  geom_bar(stat = "identity", width=0.7) + 
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  #scale_fill_manual(values = c("#0072B2", "#D55E00")) + # colorblindness palette blue and orange
  geom_errorbar(aes(ymax = y_pos + sd, ymin=y_pos - sd), stat = "identity", width = 0.1, alpha = 0.7, position = pd) + 
  facet_wrap(~Date)+ 
  scale_alpha_manual(values=c(seq(0.3,1, length.out = 3))) + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(from = 0, to = 100, by = 10)) + 
  ylim(0,110)+
  ylab("mean proportion of total acini area (%)") + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
proportion_Figure

#  SAVE
ggsave(file="Proportion_MaleHist_plot.pdf", proportion_Figure, width = 12, height = 8, units = c("in")) 

################################################## #
################################################## #
######   Contingency Table  & Chi-sq test  ####### #
################################################## #

# load libraries for chisq analysis and figures
library(reshape)
library(reshape2)
library(gplots)
library(corrplot)
# Run with ID combining both treatment and date
Male_hist_cont_table <- Male_hist_2  %>%  dplyr::select(Date, Treatment, perc_cytes, perc_lumen, perc_zoa)
Male_hist_cont_table

Male_hist_cont_table <-
  Male_hist_2 %>%
  group_by(Treatment,Date) %>%
  summarize(
    prop_zoa    = mean(perc_zoa, na.rm = TRUE),
    prop_cytes  = mean(perc_cytes, na.rm = TRUE),
    prop_lumen  = mean(perc_lumen, na.rm = TRUE))
Male_hist_cont_table # view the table

Male_contingencyTable <- Male_hist_cont_table %>% pivot_longer(
  Male_hist_cont_table,
  cols = starts_with("prop"),
  names_to = "hist_meas",
  values_to = "value")
Male_contingencyTable$ID <- paste(substr(Male_contingencyTable$Treatment,1,1), substr(Male_contingencyTable$Date,5,8), sep='_')
Male_contingencyTable2 <- Male_contingencyTable %>% dplyr::select(c(ID,value,hist_meas))
Male_contingencyTable2 <- Male_contingencyTable2[,-1]

Male_contingencyTable2 <- cast(Male_contingencyTable2, hist_meas~ID) # cast to move categorical treatments as columns
MATRIX_Male_contingencyTable2 <- as.table(as.matrix(Male_contingencyTable2)) # 1. convert the data as a table

balloonplot(t(MATRIX_Male_contingencyTable2), main ="ID", xlab ="", ylab="", label = FALSE, show.margins = FALSE) # 2. Graph

chisq <- chisq.test(MATRIX_Male_contingencyTable2) # Chi squared test
chisq # note: df = num row -1 * num col -1 - Chi Sq test is NOT significant
round(chisq$observed,2) # Observed rounded to 2 sig figs
round(chisq$expected,2) # Expected  rounded to 2 sig figs
round(chisq$residuals, 3) # Pearson residuals can be easily extracted from the output of the function chisq.test():
corrplot(chisq$residuals, is.cor = FALSE) # visualize Pearson residuals using the package corrplot:
contrib <- 100*chisq$residuals^2/chisq$statistic # Contibution in percentage (%)
round(contrib, 3)
corrplot(contrib, is.cor = FALSE) # Visualize the contribution
# Although chi-square was not significant - 
# observed patterns reinforce results from the two-way ANOVA 
# showing a decrease in spermatocytes with time and an increase 
# in spermatozoa with time. This 'seesaw' effect suggests that 
# male P. generosa are able to continue reproductive development 
# regardless of pCO2 treatment

# NOTE: since we do see an affect of time in ANOVA - run with time and without treatment chi-squared
# Run with ID just as time
Male_hist_cont_table3 <-
  Male_hist_2 %>%
  group_by(Date) %>%
  summarize(
    prop_zoa    = mean(perc_zoa, na.rm = TRUE),
    prop_cytes  = mean(perc_cytes, na.rm = TRUE),
    prop_lumen  = mean(perc_lumen, na.rm = TRUE))
Male_hist_cont_table3 <- Male_hist_cont_table3 %>% pivot_longer(
  Male_hist_cont_table,
  cols = starts_with("prop"),
  names_to = "hist_meas",
  values_to = "value")
Male_hist_cont_table3 # view the table
Male_hist_cont_table3$Date <- as.character(substr(Male_hist_cont_table3$Date,5,8))
Male_hist_cont_table3 <- cast(Male_hist_cont_table3, hist_meas~Date) # cast to move categorical treatments as columns
MATRIX_Male_hist_cont_table3 <- as.table(as.matrix(Male_hist_cont_table3)) # 1. convert the data as a table

balloonplot(t(MATRIX_Male_hist_cont_table3), main ="Date", xlab ="", ylab="", label = FALSE, show.margins = FALSE) # 2. Graph
chisq.date <- chisq.test(MATRIX_Male_hist_cont_table3) # Chi squared test
chisq.date # note: df = num row -1 * num col -1 - Chi Sq test is NOT significant
round(chisq.date$observed,2) # Observed rounded to 2 sig figs
round(chisq.date$expected,2) # Expected  rounded to 2 sig figs
round(chisq.date$residuals, 3) # Pearson residuals can be easily extracted from the output of the function chisq.test():
corrplot(chisq.date$residuals, is.cor = FALSE) # visualize Pearson residuals using the package corrplot:
contrib.date <- 100*chisq.date$residuals^2/chisq.date$statistic # Contibution in percentage (%)
round(contrib.date, 3)
corrplot(contrib.date, is.cor = FALSE) # Visualize the contribution


# contingency table with staging ID and treatment
MaleHistStage_merge # dataset
Male_hist_cont_table_STAGE <- MaleHistStage_merge  %>%  dplyr::select(Date, Treatment, Stage_ID, perc_cytes, perc_lumen, perc_zoa)
Male_hist_cont_table_STAGE$ID <- paste(substr(Male_hist_cont_table_STAGE$Treatment, 1,3), Male_hist_cont_table_STAGE$Stage_ID, sep="_")
Male_hist_cont_table_STAGE <- Male_hist_cont_table_STAGE %>% dplyr::select(ID, perc_cytes, perc_lumen, perc_zoa)
Male_hist_cont_table_STAGE <-
  Male_hist_cont_table_STAGE %>%
  group_by(ID) %>%
  summarize(
    prop_zoa    = mean(perc_zoa, na.rm = TRUE),
    prop_cytes  = mean(perc_cytes, na.rm = TRUE),
    prop_lumen  = mean(perc_lumen, na.rm = TRUE))
Male_hist_cont_table_STAGE <- Male_hist_cont_table_STAGE %>% pivot_longer(
  Male_hist_cont_table_STAGE,
  cols = starts_with("prop"),
  names_to = "hist_meas",
  values_to = "value")
Male_hist_cont_table_STAGE # view the table

Cont_table_STAGE <- cast(Male_hist_cont_table_STAGE, hist_meas~ID) # cast to move categorical treatments as columns
MATRIX_Cont_table_STAGE <- as.table(as.matrix(Cont_table_STAGE)) # 1. convert the data as a table

balloonplot(t(MATRIX_Cont_table_STAGE), main ="ID", xlab ="", ylab="", label = FALSE, show.margins = FALSE) # 2. Graph
chisq.STAGE <- chisq.test(MATRIX_Cont_table_STAGE) # Chi squared test
chisq.STAGE # note: df = num row -1 * num col -1 - Chi Sq test is NOT significant
round(chisq.STAGE$observed,2) # Observed rounded to 2 sig figs
round(chisq.STAGE$expected,2) # Expected  rounded to 2 sig figs
round(chisq.STAGE$residuals, 3) # Pearson residuals can be easily extracted from the output of the function chisq.test():
corrplot(chisq.STAGE$residuals, is.cor = FALSE) # visualize Pearson residuals using the package corrplot:
contrib.STAGE<- 100*chisq.STAGE$residuals^2/chisq.STAGE$statistic # Contibution in percentage (%)
round(contrib.STAGE, 3)
corrplot(contrib.STAGE, is.cor = FALSE) # Visualize the contribution

################################################## #
################################################## #
######   CONVERT TO MEANS TO TACKLE  ############# #
################################################## #
########MIGHT BE THE WRONG APPROACH ############## #

# CALC THE MEAN FOR EACH SAMPLE-------------------------------------------------------------------------
Male_hist_2$ID2 <- as.numeric(Male_hist_2$ID)
Means_Table <- Male_hist_2 %>%
  dplyr::select(-'acini_segment') %>% # remove unecessary columns
  group_by(ID,Treatment,Date) %>%
  summarize(
            mean_zoa    = mean(zoa, na.rm = TRUE),
            mean_cytes  = mean(cytes, na.rm = TRUE),
            mean_lumen  = mean(lumen, na.rm = TRUE),
            mean_total_acini_area = mean(total_area, na.rm = TRUE),
            mean_AREA   = mean(TOTAL_AREA, na.rm = TRUE),
            prop_zoa    = mean(perc_zoa, na.rm = TRUE),
            prop_cytes  = mean(perc_cytes, na.rm = TRUE),
            prop_lumen  = mean(perc_lumen, na.rm = TRUE),
            prop_zoa    = mean(perc_zoa, na.rm = TRUE),
            num_acini   = n())
Means_Table # view the table

# PIVOT LONGER-----------------------------------------------------------------------------------------
Means_Table_long <- Means_Table %>% 
  tidyr::pivot_longer(cols = c(4:11), names_to='scoring_metric', values_to='means')

#  PLOTS----------------------------------------------------------------------------------------
# ALL PLOTS (using Means_Table_long)
# order variable Date_Treat
Means_Table_long$Date_Treat <- paste((substr(Means_Table_long$Date,5,8)),Means_Table_long$Treatment, sep ="_")
list(Means_Table_long$Date_Treat)
# order the levels of the factor date_treatment to order correly in the plot
Means_Table_long$Date_Treat <- factor(Means_Table_long$Date_Treat,levels = c("0123_Ambient", "0123_Low", "0221_Ambient", "0221_Low"))
# plot
plot2 <- plot %>% ggadd(shape ="Treatment",fill = "white") %>% ggadd("jitter", size = 4,shape ="Treatment",fill = "white") +
  facet_wrap( ~ scoring_metric, ncol=2, scales = "free") + theme_classic()
plot2

# INDIVIIDUAL PLOTS (using Means_Table)
# MEAN AREA OF ACINI
Means_Table$Date <- as.character(Means_Table$Date)
plot_area <- ggplot(Means_Table, aes(x=Date, y=mean_AREA, color=Treatment, fill - white)) +
              theme_bw() +
              geom_boxplot() +
              ylab("acini area (µm)") +
              scale_x_discrete(name ="exposure time (days)", labels=c("72","93 + 8 day recovery")) +
              annotate("text", size = 12, x=0.55, y=1300, label= "A")
pA <- plot_area %>% ggadd("jitter", size = 3, fill = "white") + theme(legend.position="none")
# plot proportion spermatozoa 
plot_zoa_proportion <- ggplot(Means_Table, aes(x=Date, y=prop_zoa, color=Treatment, fill - white)) +
  theme_bw() +
  geom_boxplot() +
  ylab("spermatozoa area proportion (%)") +
  scale_x_discrete(name ="exposure time (days)", labels=c("72","93 + 8 day recovery")) +
  annotate("text", size = 12, x=0.55, y=45, label= "B")
pB <- plot_zoa_proportion %>% ggadd("jitter", size = 3, fill = "white") + theme(legend.position="none")
# plot proportion spermatocytes 
plot_cytes_proportion <- ggplot(Means_Table, aes(x=Date, y=prop_cytes, color=Treatment, fill - white)) +
  theme_bw() +
  geom_boxplot() +
  ylab("spermatocytes area proportion (%)") +
  scale_x_discrete(name ="exposure time (days)", labels=c("72","93 + 8 day recovery")) +
  annotate("text", size = 12, x=0.55, y=95, label= "C")
pC <- plot_cytes_proportion %>% ggadd("jitter", size = 3, fill = "white") + theme(legend.position="none")
# plot proportion lumen
plot_lumen_proportion <- ggplot(Means_Table, aes(x=Date, y=prop_lumen, color=Treatment, fill - white)) +
  theme_bw() +
  geom_boxplot() +
  ylab("lumen area proportion (%)") +
  scale_x_discrete(name ="exposure time (days)", labels=c("72","93 + 8 day recovery")) +
  annotate("text", size = 12, x=0.55, y=22, label= "D")
pD <- plot_lumen_proportion %>% ggadd("jitter", size = 3, fill = "white") + theme(legend.position="none")
# gird plots
Final_Plot_Grid <- grid.arrange(pA, pB,pC, pD, ncol =2, nrow = 2)
Final_Plot_Grid # view plot

# Fig Grid 2
Means_Table[11,]
MeansTable_2 <- Means_Table[-c(11), ]
# plot  spermatozoa 
plot_zoa_area <- ggplot(MeansTable_2, aes(x=Date, y=mean_zoa, color=Treatment, fill - white)) +
  theme_bw() +
  geom_boxplot() +
  ylab("mean area spermatozoa (µm^3)") +
  scale_x_discrete(name ="exposure time (days)", labels=c("72","93 + 8 day recovery")) +
  annotate("text", size = 12, x=0.55, y=320, label= "B")
pB2 <- plot_zoa_area %>% ggadd("jitter", size = 3, fill = "white") + theme(legend.position="none")
# plot  spermatocytes 
plot_cytes_area <- ggplot(MeansTable_2, aes(x=Date, y=mean_cytes, color=Treatment, fill - white)) +
  theme_bw() +
  geom_boxplot() +
  ylab("mean area spermatocytes (µm^3)") +
  scale_x_discrete(name ="exposure time (days)", labels=c("72","93 + 8 day recovery")) +
  annotate("text", size = 12, x=0.55, y=975, label= "C")
pC2 <- plot_cytes_area %>% ggadd("jitter", size = 3, fill = "white") + theme(legend.position="none")
# plot  lumen
plot_lumen_area <- ggplot(MeansTable_2, aes(x=Date, y=mean_lumen, color=Treatment, fill - white)) +
  theme_bw() +
  geom_boxplot() +
  ylab("mean area lumen (µm^3)") +
  scale_x_discrete(name ="exposure time (days)", labels=c("72","93 + 8 day recovery")) +
  annotate("text", size = 12, x=0.55, y=200, label= "D")
pD2 <- plot_lumen_area %>% ggadd("jitter", size = 3, fill = "white") + theme(legend.position="none")
# gird plots
Final_Plot_Grid2 <- grid.arrange(pA, pB2,pC2, pD2, ncol =2, nrow = 2)
Final_Plot_Grid2 # view plot

#  SAVE
ggsave(file="Grid_plot.pdf", Final_Plot_Grid2, width = 12, height = 8, units = c("in")) 

################################################################# #
################################################################# #
# STATISTICS WITH MEAN DATA------------------------------------------------------------------------------------------
################################################################# #
################################################################# #
################################################################# #

Means_Table$Treatment <- as.factor(Means_Table$Treatment)
Means_Table$Date <- as.factor(Means_Table$Date)
par(mfrow=c(2,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots

# Two-Way ANOVA (treat×time)
typeof(Means_Table$Date) # currenty an integer - must change to a character
Means_Table$Date <- as.character(Means_Table$Date)


## TOTAL_AREA test (lumen + cytes + zoa) ---------------------------------
mod_total_area_TIME  <- aov(mean_AREA~Treatment*Date, data = Means_Table)
anova(mod_total_area_TIME) 
shapiro.test(residuals(mod_total_area_TIME)) #  normal residuals p-value = 0.001496 (driven by outlier?)
leveneTest(mod_total_area_TIME) # p = 0.3766
hist(residuals(mod_total_area_TIME)) #plot histogram of residuals
boxplot(residuals(mod_total_area_TIME)) #plot boxplot of residuals
plot(fitted(mod_total_area_TIME),residuals(mod_total_area_TIME))
qqnorm(residuals(mod_total_area_TIME)) # qqplot

## PROPORTION zoa test  ---------------------------------
mod_prop_zoa_TIME  <- aov(prop_zoa~Treatment*Date, data = Means_Table)
anova(mod_prop_zoa_TIME) 
shapiro.test(residuals(mod_prop_zoa_TIME)) #  normal residuals p-value = 0.2917
leveneTest(mod_prop_zoa_TIME) # p = 0.7
hist(residuals(mod_prop_zoa_TIME)) #plot histogram of residuals
boxplot(residuals(mod_prop_zoa_TIME)) #plot boxplot of residuals
plot(fitted(mod_prop_zoa_TIME),residuals(mod_prop_zoa_TIME))
qqnorm(residuals(mod_prop_zoa_TIME)) # qqplot

## PROPORTION cytes test  ---------------------------------
mod_prop_cytes_TIME  <- aov(prop_cytes~Treatment*Date, data = Means_Table)
anova(mod_prop_cytes_TIME) 
shapiro.test(residuals(mod_prop_cytes_TIME)) #  normal residuals p-value = 0.1608
leveneTest(mod_prop_cytes_TIME) # p = 0.5509
hist(residuals(mod_prop_cytes_TIME)) #plot histogram of residuals
boxplot(residuals(mod_prop_cytes_TIME)) #plot boxplot of residuals
plot(fitted(mod_prop_cytes_TIME),residuals(mod_prop_cytes_TIME))
qqnorm(residuals(mod_prop_cytes_TIME)) # qqplot

## PROPORTION lumen test  ---------------------------------
mod_prop_lumen_TIME  <- aov(prop_lumen~Treatment*Date, data = Means_Table)
anova(mod_prop_lumen_TIME) 
shapiro.test(residuals(mod_prop_lumen_TIME)) #  normal residuals p-value = 7.261e-06
leveneTest(mod_prop_lumen_TIME) # p = 0.6463
hist(residuals(mod_prop_lumen_TIME)) #plot histogram of residuals
boxplot(residuals(mod_prop_lumen_TIME)) #plot boxplot of residuals
plot(fitted(mod_prop_lumen_TIME),residuals(mod_prop_lumen_TIME))
qqnorm(residuals(mod_prop_lumen_TIME)) # qqplot

## total acini area test 2 (total of whole area + white space) ---------------------------------
mod_total_acini_TIME  <- aov(mean_total_acini_area~Treatment*Date, data = Means_Table)
anova(mod_total_acini_TIME) # p = 0.9391
shapiro.test(residuals(mod_total_acini_TIME)) #  normal residuals p-value = 7.261e-06
leveneTest(mod_total_acini_TIME) # p = 0.6463
hist(residuals(mod_total_acini_TIME)) #plot histogram of residuals
boxplot(residuals(mod_total_acini_TIME)) #plot boxplot of residuals
plot(fitted(mod_total_acini_TIME),residuals(mod_total_acini_TIME))
qqnorm(residuals(mod_total_acini_TIME)) # qqplot

## cytes test---------------------------------
mod_cytes_TIME  <- aov(mean_cytes~Treatment*Date, data = MeansTable_2)
anova(mod_cytes_TIME) # p = 0.9391
shapiro.test(residuals(mod_cytes_TIME)) #  normal residuals p-value = 7.261e-06
leveneTest(mod_cytes_TIME) # p = 0.6463
hist(residuals(mod_cytes_TIME)) #plot histogram of residuals
boxplot(residuals(mod_cytes_TIME)) #plot boxplot of residuals
plot(fitted(mod_cytes_TIME),residuals(mod_cytes_TIME))
qqnorm(residuals(mod_cytes_TIME)) # qqplot

## zoa test---------------------------------
mod_zoa_TIME  <- aov(mean_zoa~Treatment*Date, data = MeansTable_2)
anova(mod_zoa_TIME) # p = 0.01051 * (time)
shapiro.test(residuals(mod_zoa_TIME)) #  normal residuals p-value = 0.08541
leveneTest(mod_zoa_TIME) # p = 0.1719
hist(residuals(mod_zoa_TIME)) #plot histogram of residuals
boxplot(residuals(mod_zoa_TIME)) #plot boxplot of residuals
plot(fitted(mod_zoa_TIME),residuals(mod_zoa_TIME))
qqnorm(residuals(mod_zoa_TIME)) # qqplot

TukeyHSD(mod_zoa_TIME)

## lumen test---------------------------------
mod_lumen_TIME  <- aov(mean_lumen~Treatment*Date, data = MeansTable_2)
anova(mod_lumen_TIME) # p = 0.06423 . (treat) ; 0.09255 . (time)
shapiro.test(residuals(mod_lumen_TIME)) #  normal residuals p-value = 0.2685
leveneTest(mod_lumen_TIME) # p = 0.1227
hist(residuals(mod_lumen_TIME)) #plot histogram of residuals
boxplot(residuals(mod_lumen_TIME)) #plot boxplot of residuals
plot(fitted(mod_lumen_TIME),residuals(mod_lumen_TIME))
qqnorm(residuals(mod_lumen_TIME)) # qqplot

TukeyHSD(mod_lumen_TIME)

