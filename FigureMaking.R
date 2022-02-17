# E Reichert, 2022
# Plot figures for NICD COVID-19 lag analysis manuscript

#load necessary libraries
library(lubridate)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)
library(ggridges)

# Here, we first run a script that reads in and cleans surveillance data (CaseAdmitAvg) 
# containing daily COVID-19 case counts, hospital admissions with COVID-19, and within-hospital
# COVID-19 attributable deaths (3 separate files), disaggregated by age category. 
source("Documents/2021 CCDD research/NICDdatacleaning.R")

# Now run case to hospitalization lag analysis file
source("Documents/2021 CCDD research/HospLag.R")
# Now run case/hosp to death lag analysis file
source("Documents/2021 CCDD research/DeathLag.R")

############################
#Fig. 1 Age-Stratified Epi Curves
############################

## The code requires a dataset (here, CaseAdmitAvg) with Date (Date), Age Category (agecat), 
## and 07da average variables for confirmed COVID-19 cases (cases_07da), hospitalizations 
## (hosp_07da), and deaths (deaths_07da). Dates used to represent specific VOC-dominated
## epidemic waves in this script are specific to South Africa.

#Plot deaths (07da average) over time
Deaths_curve <- ggplot(CaseAdmitAvg) +
  geom_line(aes(x = Date, y = deaths_07da, color = agecat), size = 0.7) +
  #add vertical lines demarcating t0 and t1 for each VOC dominated wave
  geom_vline(xintercept = CaseAdmitAvg$Date[27], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[222], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[237], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[376], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[419], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[593], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[608], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[690], linetype = "dashed") +
  theme_classic() + scale_color_brewer(palette = "Paired") + ylab("7d Average Deaths") +
  #add labels for each VOC dominated wave
  annotate(geom="text", x=CaseAdmitAvg$Date[27]+days(30), y=250, label="D614G",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[237]+days(30), y=250, label="Beta",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[419]+days(30), y=250, label="Delta",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[608]+days(32), y=250, label="Omicron",
           color="black", size = 3) +
  ggtitle("C.")
Deaths_curve <- Deaths_curve + theme(legend.position = "bottom") + labs(col = "Age Category")
Deaths_curve

#plot hospital admissions (07da average) over time
Hosp_curve <- ggplot(CaseAdmitAvg) +
  geom_line(aes(x = Date, y = hosp_07da, color = agecat), size = 0.7) +
  #add vertical lines demarcating t0 and t1 for each VOC dominated wave
  geom_vline(xintercept = CaseAdmitAvg$Date[27], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[222], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[237], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[376], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[419], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[593], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[608], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[690], linetype = "dashed") +
  theme_classic() + scale_color_brewer(palette = "Paired") + ylab("7d Average Hospitalizations") +
  #add labels for each VOC dominated wave
  annotate(geom="text", x=CaseAdmitAvg$Date[27]+days(28), y=900, label="D614G",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[237]+days(28), y=900, label="Beta",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[419]+days(28), y=900, label="Delta",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[608]+days(30), y=900, label="Omicron",
           color="black", size = 3) +
  ggtitle("B.")
Hosp_curve <- Hosp_curve + theme(legend.position = "none")
Hosp_curve

Case_curve <- ggplot(CaseAdmitAvg) +
  geom_line(aes(x = Date, y = cases_07da, color = agecat), size = 0.7) +
  #add vertical lines demarcating t0 and t1 for each VOC dominated wave
  geom_vline(xintercept = CaseAdmitAvg$Date[27], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[222], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[237], linetype = "dashed") +
  geom_vline(xintercept = CaseAdmitAvg$Date[376], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[419], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[593], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[608], linetype = "dashed") +  
  geom_vline(xintercept = CaseAdmitAvg$Date[690], linetype = "dashed") +
  theme_classic() + scale_color_brewer(palette = "Paired") + ylab("7d Average Cases") +
  #add labels for each VOC dominated wave
  annotate(geom="text", x=CaseAdmitAvg$Date[27]+days(28), y=6000, label="D614G",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[237]+days(28), y=6000, label="Beta",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[419]+days(28), y=6000, label="Delta",
           color="black", size = 3) +
  annotate(geom="text", x=CaseAdmitAvg$Date[608]+days(30), y=6000, label="Omicron",
           color="black", size = 3) +
  ggtitle("A.")
Case_curve <- Case_curve + theme(legend.position= "none") 
Case_curve

#create multi-faceted plot with all 3 visualized
grid.arrange(Case_curve, Hosp_curve, Deaths_curve, nrow = 3, heights = c(1,1,1.3))

############################
#Fig. 2 Age-Stratified Distribution of COVID-19 Infections, for each VOC
############################

#Step 1 - for each VOC's dataset, turn it into long format so that each record
# represents 1 case
WT_long <- WT_deaths[rep(1:nrow(WT_deaths), WT_deaths$CaseCount), ] %>%
  dplyr::select(Date, agecat)

#create ridgeplot of cases over time, stratified by age category
WT_ridge <- ggplot(WT_long, aes(x = Date, y = agecat)) +
  geom_density_ridges(aes(fill = agecat), bandwidth = 230000, color = "white") +
  ggtitle("D614G: Apr 2020-Oct 2020") + theme_classic() +
  geom_vline(xintercept = WT_deaths$Date[107], color = "black", linetype = "dashed") +
  scale_fill_brewer(palette = "Paired") + ylab("Age Category")
WT_ridge <- WT_ridge + theme(legend.position = "none")
#WT - peak case # on 2020-07-20, added vertical line indicator

Beta_long <- Beta_deaths[rep(1:nrow(Beta_deaths), Beta_deaths$CaseCount), ] %>%
  dplyr::select(Date, agecat)

Beta_ridge <- ggplot(Beta_long, aes(x = Date, y = agecat)) +
  geom_density_ridges(aes(fill = agecat), bandwidth = 230000, color = "white") +
  ggtitle("Beta: Nov 2020-Mar 2021") + theme_classic() +
  geom_vline(xintercept = Beta_deaths$Date[65], color = "black", linetype = "dashed") +
  scale_fill_brewer(palette = "Paired") + ylab("Age Category")
Beta_ridge <- Beta_ridge + theme(legend.position = "none")
#beta - peak case # on 2021-01-04, added vertical line indicator

Delta_long <- Delta_deaths[rep(1:nrow(Delta_deaths), Delta_deaths$CaseCount), ] %>%
  dplyr::select(Date, agecat)

Delta_ridge <- ggplot(Delta_long, aes(x = Date, y = agecat)) +
  geom_density_ridges(aes(fill = agecat), bandwidth = 230000, color = "white") +
  ggtitle("Delta: May 2021-Oct 2021") + theme_classic() +
  geom_vline(xintercept = Delta_deaths$Date[590], color = "black", linetype = "dashed") +
  scale_fill_brewer(palette = "Paired") + ylab("Age Category")
Delta_ridge <- Delta_ridge + theme(legend.position = "none")
#delta - peak case # on 2021-07-05, added vertical line indicator

Omi_long <- Omi_deaths[rep(1:nrow(Omi_deaths), Omi_deaths$CaseCount), ] %>%
  dplyr::select(Date, agecat)

Omi_ridge <- ggplot(Omi_long, aes(x = Date, y = agecat)) +
  geom_density_ridges(aes(fill = agecat), bandwidth = 230000, color = "white") +
  ggtitle("Omicron: Nov 2021-Jan 2022") + theme_classic() +
  geom_vline(xintercept = Omi_deaths$Date[197], color = "black", linetype = "dashed") +
  scale_fill_brewer(palette = "Paired") + ylab("Age Category")
Omi_ridge <- Omi_ridge + theme(legend.position  = "none")
#omicron - peak case # on 2021-12-13, added vertical line indicator

#create aggregate plot
grid.arrange(WT_ridge, Beta_ridge, Delta_ridge, Omi_ridge, nrow = 1)

############################
#Fig. 3 Heatmap of Distance Correlations for all associations
# Case-Hosp, Hosp-Death, and Case-Death, by VOC and Age Category
############################

VOC_corr_CasetoHosp <- VOC_corr_CasetoHosp %>%
  mutate(Type = rep("A. Cases to Hospitalizations")) 
VOC_corr_CasetoDeath <- VOC_corr_CasetoDeath %>%
  mutate(Type = rep("C. Cases to Deaths"))
VOC_corr_HosptoDeath <- VOC_corr_HosptoDeath %>%
  mutate(Type = rep("B. Hospitalizations to Deaths"))
VOC_corr <- rbind(VOC_corr_CasetoHosp, VOC_corr_CasetoDeath, VOC_corr_HosptoDeath)

VOC_corr_max_CasetoHosp <- VOC_corr_max_CasetoHosp %>%
  mutate(Type = rep("A. Cases to Hospitalizations"))
VOC_corr_max_CasetoDeath <- VOC_corr_max_CasetoDeath %>%
  mutate(Type = rep("C. Cases to Deaths"))
VOC_corr_max_HosptoDeath <- VOC_corr_max_HosptoDeath %>%
  mutate(Type = "B. Hospitalizations to Deaths")
VOC_corr_max <- rbind(VOC_corr_max_CasetoHosp, VOC_corr_max_HosptoDeath, VOC_corr_max_CasetoDeath)

VOC_corr_90_CasetoHosp <- VOC_corr_90_CasetoHosp %>%
  mutate(Type = rep("A. Cases to Hospitalizations"))
VOC_corr_90_CasetoDeath <- VOC_corr_90_CasetoDeath %>%
  mutate(Type = rep("C. Cases to Deaths"))
VOC_corr_90_HosptoDeath <- VOC_corr_90_HosptoDeath %>%
  mutate(Type = rep("B. Hospitalizations to Deaths"))
VOC_corr_90 <- rbind(VOC_corr_90_CasetoHosp, VOC_corr_90_CasetoDeath, VOC_corr_90_HosptoDeath)

P_all <- ggplot(VOC_corr, aes(x = LagTime, y = agecat)) +
  geom_tile(aes(fill = Corr)) + 
  geom_point(data=VOC_corr_max, size=1, fill=NA, colour="black",
             aes(x = LagTime, y = agecat)) +
  geom_line(data = VOC_corr_90, aes(x = LagTime, y = agecat), colour = "black", size = 0.4, linetype = "dotted") +
  scale_fill_viridis(option = "mako") +
  xlab("Lag Interval (Days)") + 
  ylab("Age Category") + 
  facet_grid(VOC ~ Type, scales = "free")  + theme_classic() + labs(fill = expression(paste("Distance \nCorrelation")))
P_all <- P_all + theme(legend.position = "bottom", text = element_text(size = 12))

