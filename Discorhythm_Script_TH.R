# This script repeats DiscoRhythm analysis for four different conditions of the RNA-seq PAC2 dataset

library(tidyverse)
library(DiscoRhythm)
library(VennDiagram)
options(max.print = 1000000)

#################################################################################

SDWT <- read.table("SDWTREPS.txt", header=TRUE, sep = "\t")

colnames(SDWT) <- c("IDs", "CT3_1_A", "CT3_2_B", "CT3_3_C", "CT9_4_D", "CT9_5_E", "CT9_6_F", "CT15_7_G", "CT15_8_H", "CT15_9_I", "CT21_10_J", "CT21_11_K", "CT21_12_L", "CT27_13_M", "CT27_14_N", "CT27_15_O", "CT33_16_P", "CT33_17_Q", "CT33_18_R", "CT39_19_A", "CT39_20_B", "CT39_21_C", "CT45_22_D", "CT45_23_E", "CT45_24_F")

se <- discoDFtoSE(SDWT)

discoODAres <- discoODAs(se, period=24, method="JTK", ncores=1, circular_t=FALSE)

x <- discoODAres[["JTK"]]
x <- rownames_to_column(x,"name")
y <- select(x,name,pvalue,period)
z <- filter(y, pvalue < 0.05, period == 24)
s <- select(z,name,pvalue)
write_delim(s, file = "trueSDWT")

discoBatch(indata=se, report="discoSDWT.html", ncores=1, main_per=24, timeType="linear", cor_threshold=3, cor_method="pearson", cor_threshType="sd", pca_threshold=3, pca_scale=TRUE, pca_pcToCut=paste0("PC",seq_len(4)), aov_method="None", aov_pcut=0.05, aov_Fcut=0, avg_method="Median", osc_method="JTK", osc_period=24)

#################################################################################

LDWT <- read.table("LDWTREPS.txt", header=TRUE, sep = "\t")

colnames(LDWT) <- c("IDs", "CT3_1_A", "CT3_2_B", "CT3_3_C", "CT9_4_D", "CT9_5_E", "CT9_6_F", "CT15_7_G", "CT15_8_H", "CT15_9_I", "CT21_10_J", "CT21_11_K", "CT21_12_L", "CT27_13_M", "CT27_14_N", "CT27_15_O", "CT33_16_P", "CT33_17_Q", "CT33_18_R", "CT39_19_A", "CT39_20_B", "CT39_21_C", "CT45_22_D", "CT45_23_E", "CT45_24_F")

se <- discoDFtoSE(LDWT)

discoODAres <- discoODAs(se, period=24, method="JTK", ncores=1, circular_t=FALSE)

x <- discoODAres[["JTK"]]
x <- rownames_to_column(x,"name")
y <- select(x,name,pvalue,period)
z <- filter(y, pvalue < 0.05, period == 24)
s <- select(z,name,pvalue)
write_delim(s, file = "trueLDWT")

discoBatch(indata=se, report="discoLDWT.html", ncores=1, main_per=24, timeType="linear", cor_threshold=3, cor_method="pearson", cor_threshType="sd", pca_threshold=3, pca_scale=TRUE, pca_pcToCut=paste0("PC",seq_len(4)), aov_method="None", aov_pcut=0.05, aov_Fcut=0, avg_method="Median", osc_method="JTK", osc_period=24)

#################################################################################

SDCLOCK <- read.table("SDCLOCKREPS.txt", header=TRUE, sep = "\t")

colnames(SDCLOCK) <- c("IDs", "CT3_1_A", "CT3_2_B", "CT3_3_C", "CT9_4_D", "CT9_5_E", "CT9_6_F", "CT15_7_G", "CT15_8_H", "CT15_9_I", "CT21_10_J", "CT21_11_K", "CT21_12_L", "CT27_13_M", "CT27_14_N", "CT27_15_O", "CT33_16_P", "CT33_17_Q", "CT33_18_R", "CT39_19_A", "CT39_20_B", "CT39_21_C", "CT45_22_D", "CT45_23_E", "CT45_24_F")

se <- discoDFtoSE(SDCLOCK)

discoODAres <- discoODAs(se, period=24, method="JTK", ncores=1, circular_t=FALSE)

x <- discoODAres[["JTK"]]
x <- rownames_to_column(x,"name")
y <- select(x,name,pvalue,period)
z <- filter(y, pvalue < 0.05, period == 24)
s <- select(z,name,pvalue)
write_delim(s, file = "trueSDCLOCK")

discoBatch(indata=se, report="discoSDCLOCK.html", ncores=1, main_per=24, timeType="linear", cor_threshold=3, cor_method="pearson", cor_threshType="sd", pca_threshold=3, pca_scale=TRUE, pca_pcToCut=paste0("PC",seq_len(4)), aov_method="None", aov_pcut=0.05, aov_Fcut=0, avg_method="Median", osc_method="JTK", osc_period=24)

#################################################################################

LDCLOCK <- read.table("LDCLOCKREPS.txt", header=TRUE, sep = "\t")

colnames(LDCLOCK) <- c("IDs", "CT3_1_A", "CT3_2_B", "CT3_3_C", "CT9_4_D", "CT9_5_E", "CT9_6_F", "CT15_7_G", "CT15_8_H", "CT15_9_I", "CT21_10_J", "CT21_11_K", "CT21_12_L", "CT27_13_M", "CT27_14_N", "CT27_15_O", "CT33_16_P", "CT33_17_Q", "CT33_18_R", "CT39_19_A", "CT39_20_B", "CT39_21_C", "CT45_22_D", "CT45_23_E", "CT45_24_F")

se <- discoDFtoSE(LDCLOCK)

discoODAres <- discoODAs(se, period=24, method="JTK", ncores=1, circular_t=FALSE)

x <- discoODAres[["JTK"]]
x <- rownames_to_column(x,"name")
y <- select(x,name,pvalue,period)
z <- filter(y, pvalue < 0.05, period == 24)
s <- select(z,name,pvalue)
write_delim(s, file = "trueLDCLOCK")

discoBatch(indata=se, report="discoLDCLOCK.html", ncores=1, main_per=24, timeType="linear", cor_threshold=3, cor_method="pearson", cor_threshType="sd", pca_threshold=3, pca_scale=TRUE, pca_pcToCut=paste0("PC",seq_len(4)), aov_method="None", aov_pcut=0.05, aov_Fcut=0, avg_method="Median", osc_method="JTK", osc_period=24)

#################################################################################

 trueSDWT <- read.csv("~/ZebrafishJTKcycle/trueSDWT", sep="")
   View(trueSDWT)
 trueSDCLOCK <- read.csv("~/ZebrafishJTKcycle/trueSDCLOCK", sep="")
   View(trueSDCLOCK)
 trueLDWT <- read.csv("~/ZebrafishJTKcycle/trueLDWT", sep="")
   View(trueLDWT)
 trueLDCLOCK <- read.csv("~/ZebrafishJTKcycle/trueLDCLOCK", sep="")
   View(trueLDCLOCK)
z <- unique(trueLDCLOCK$name)
x <- unique(trueSDCLOCK$name)
y <- unique(trueLDWT$name)
s <- unique(trueSDWT$name)
#################################################################################
