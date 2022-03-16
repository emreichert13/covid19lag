#E Reichert 2022
# COVID-19 Case to Death and Hospitalization to Death Lag Analysis 
# using NICD Surveillance Data

#load necessary libraries
library(readr)
library(readxl)
library(lubridate)
library(ggplot2)
library(zoo)
library(tidyr)
library(stringr)
library(viridis)
library(energy)
library(dplyr)
library(EnvStats)
library(ggridges)

# Here, we first run a script that reads in and cleans surveillance data (CaseAdmitAvg) 
# containing daily COVID-19 case counts, hospital admissions with COVID-19, and within-hospital
# COVID-19 attributable deaths (3 separate files), disaggregated by age category. 
source("Documents/2021 CCDD research/NICDdatacleaning.R")

# Now run case to hospitalization lag analysis file
source("Documents/2021 CCDD research/HospLag.R")

###################################################
# Cases to Deaths
###################################################

# Function to create lead time variables for deaths
# note, each death_lead[x] variable contains 
# COVID-19 death counts for x days in the future
CreateLeads2 <- function(x) {
  x %>%
    mutate(Death_lead1 = lead(DeathCount, 1)) %>%
    mutate(Death_lead2 = lead(DeathCount, 2)) %>%
    mutate(Death_lead3 = lead(DeathCount, 3)) %>%
    mutate(Death_lead4 = lead(DeathCount, 4)) %>%
    mutate(Death_lead5 = lead(DeathCount, 5)) %>%
    mutate(Death_lead6 = lead(DeathCount, 6)) %>%
    mutate(Death_lead7 = lead(DeathCount, 7)) %>%
    mutate(Death_lead8 = lead(DeathCount, 8)) %>%
    mutate(Death_lead9 = lead(DeathCount, 9)) %>%
    mutate(Death_lead10 = lead(DeathCount, 10)) %>%
    mutate(Death_lead11 = lead(DeathCount, 11)) %>%
    mutate(Death_lead12 = lead(DeathCount, 12)) %>%
    mutate(Death_lead13 = lead(DeathCount, 13)) %>%
    mutate(Death_lead14 = lead(DeathCount, 14)) %>%
    mutate(Death_lead15 = lead(DeathCount, 15)) %>%
    mutate(Death_lead16 = lead(DeathCount, 16)) %>%
    mutate(Death_lead17 = lead(DeathCount, 17)) %>%
    mutate(Death_lead18 = lead(DeathCount, 18)) %>%
    mutate(Death_lead19 = lead(DeathCount, 19)) %>%
    mutate(Death_lead20 = lead(DeathCount, 20)) %>%
    mutate(Death_lead21 = lead(DeathCount, 21)) %>%
    mutate(Death_lead22 = lead(DeathCount, 22)) %>%
    mutate(Death_lead23 = lead(DeathCount, 23)) %>%
    mutate(Death_lead24 = lead(DeathCount, 24)) %>%
    mutate(Death_lead25 = lead(DeathCount, 25))
}

#Apply death leads function to each VOC dataset
WT_deaths <- WT %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads2() %>% 
  dplyr::ungroup() 

Beta_deaths <- BETA %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads2() %>% 
  dplyr::ungroup() 

Delta_deaths <- DELTA %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads2() %>% 
  dplyr::ungroup()

Omi_deaths <- OMI %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads2() %>% 
  dplyr::ungroup() %>%
  filter(!is.na(cases_07da))

# Create correlation function - assess distance correlation of 
# COVID-19 cases (7da average) with # of 
# COVID-19 deaths x days in the future, where x is 1-25 days, by age cat

CreateCorr2 <- function(x) {
  x %>%
    dplyr::group_by(agecat) %>%
    summarise(corr_Deathlead1 = dcor(cases_07da[!is.na(Death_lead1)], Death_lead1[!is.na(Death_lead1)]),
              corr_Deathlead2 = dcor(cases_07da[!is.na(Death_lead2)], Death_lead2[!is.na(Death_lead2)]),
              corr_Deathlead3 = dcor(cases_07da[!is.na(Death_lead3)], Death_lead3[!is.na(Death_lead3)]),
              corr_Deathlead4 = dcor(cases_07da[!is.na(Death_lead4)], Death_lead4[!is.na(Death_lead4)]),
              corr_Deathlead5 = dcor(cases_07da[!is.na(Death_lead5)], Death_lead5[!is.na(Death_lead5)]),
              corr_Deathlead6 = dcor(cases_07da[!is.na(Death_lead6)], Death_lead6[!is.na(Death_lead6)]),
              corr_Deathlead7 = dcor(cases_07da[!is.na(Death_lead7)], Death_lead7[!is.na(Death_lead7)]),
              corr_Deathlead8 = dcor(cases_07da[!is.na(Death_lead8)], Death_lead8[!is.na(Death_lead8)]),
              corr_Deathlead9 = dcor(cases_07da[!is.na(Death_lead9)], Death_lead9[!is.na(Death_lead9)]),
              corr_Deathlead10 = dcor(cases_07da[!is.na(Death_lead10)], Death_lead10[!is.na(Death_lead10)]),
              corr_Deathlead11 = dcor(cases_07da[!is.na(Death_lead11)], Death_lead11[!is.na(Death_lead11)]),
              corr_Deathlead12 = dcor(cases_07da[!is.na(Death_lead12)], Death_lead12[!is.na(Death_lead12)]),
              corr_Deathlead13 = dcor(cases_07da[!is.na(Death_lead13)], Death_lead13[!is.na(Death_lead13)]),
              corr_Deathlead14 = dcor(cases_07da[!is.na(Death_lead14)], Death_lead14[!is.na(Death_lead14)]),
              corr_Deathlead15 = dcor(cases_07da[!is.na(Death_lead15)], Death_lead15[!is.na(Death_lead15)]),
              corr_Deathlead16 = dcor(cases_07da[!is.na(Death_lead16)], Death_lead16[!is.na(Death_lead16)]),
              corr_Deathlead17 = dcor(cases_07da[!is.na(Death_lead17)], Death_lead17[!is.na(Death_lead17)]),
              corr_Deathlead18 = dcor(cases_07da[!is.na(Death_lead18)], Death_lead18[!is.na(Death_lead18)]),
              corr_Deathlead19 = dcor(cases_07da[!is.na(Death_lead19)], Death_lead19[!is.na(Death_lead19)]),
              corr_Deathlead20 = dcor(cases_07da[!is.na(Death_lead20)], Death_lead20[!is.na(Death_lead20)]),
              corr_Deathlead21 = dcor(cases_07da[!is.na(Death_lead21)], Death_lead21[!is.na(Death_lead21)]),
              corr_Deathlead22 = dcor(cases_07da[!is.na(Death_lead22)], Death_lead22[!is.na(Death_lead22)]),
              corr_Deathlead23 = dcor(cases_07da[!is.na(Death_lead23)], Death_lead23[!is.na(Death_lead23)]),
              corr_Deathlead24 = dcor(cases_07da[!is.na(Death_lead24)], Death_lead24[!is.na(Death_lead24)]),
              corr_Deathlead25 = dcor(cases_07da[!is.na(Death_lead25)], Death_lead25[!is.na(Death_lead25)]))
}


#calculate correlation for each VOC at each lag time
WT_corr <- CreateCorr2(WT_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("D614G"))

Beta_corr <- CreateCorr2(Beta_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Beta"))

Delta_corr <- CreateCorr2(Delta_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Delta"))

Omi_corr <- CreateCorr2(Omi_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Omicron"))

#aggregate results for each VOC-dominated wave
VOC_corr_CasetoDeath <- rbind(WT_corr, Beta_corr, Delta_corr, Omi_corr)
VOC_corr_CasetoDeath$VOC <- factor(VOC_corr_CasetoDeath$VOC, levels=c('D614G','Beta','Delta','Omicron'))

#determine optimal lag time for each VOC, age category strata
#by looking at maximum value of distance correlation within each group
VOC_corr_max_CasetoDeath <- VOC_corr_CasetoDeath %>%
  group_by(agecat, VOC) %>%
  summarise(LagTime = LagTime[which.max(Corr)], Corr = Corr[which.max(Corr)]) %>%
  mutate(VOC_num = ifelse(VOC == "D614G", 1,
                          ifelse("Beta",2,
                                 ifelse(VOC == "Delta",3,
                                        ifelse(VOC == "Omicron", 4, NA)))))

#subset of correlation data for which dCor >= 0.90 (for visualization)
VOC_corr_90_CasetoDeath <- VOC_corr_CasetoDeath %>% filter(Corr >= .90)

#heatmap of distance correlations by lag time for Cases-Deaths
P3 <- ggplot(VOC_corr_CasetoDeath, aes(x = LagTime, y = agecat)) +
  geom_tile(aes(fill = Corr)) + 
  geom_point(data=VOC_corr_max_CasetoDeath, size=1, fill=NA, colour="black",
            aes(x = LagTime, y = agecat)) +
  geom_line(data = VOC_corr_90_CasetoDeath, aes(x = LagTime, y = agecat), colour = "black", size = 0.4, linetype = "dotted") +
  scale_fill_viridis(option = "mako") +
  xlab("Lag Interval (Days)") + 
  ylab("Age Category") + 
  facet_grid(VOC ~ .)  + theme_classic() + labs(title = "C. Cases to Deaths", fill = expression(paste("Distance \nCorrelation")))
P3 <- P3 + theme(plot.title = element_text(size = 12), legend.position = "bottom")

#create CFR as function of lag time
CalcDeathRate <- function(x) {
  x %>%
    dplyr::group_by(agecat, Date) %>%
    summarise(rate_deathlead1 = Death_lead1/cases_07da,
              rate_deathlead2 = Death_lead2/cases_07da,
              rate_deathlead3 = Death_lead3/cases_07da,
              rate_deathlead4 = Death_lead4/cases_07da,
              rate_deathlead5 = Death_lead5/cases_07da,
              rate_deathlead6 = Death_lead6/cases_07da,
              rate_deathlead7 = Death_lead7/cases_07da,
              rate_deathlead8 = Death_lead8/cases_07da,
              rate_deathlead9 = Death_lead9/cases_07da,
              rate_deathlead10 = Death_lead10/cases_07da,
              rate_deathlead11 = Death_lead11/cases_07da,
              rate_deathlead12 = Death_lead12/cases_07da,
              rate_deathlead13 = Death_lead13/cases_07da,
              rate_deathlead14 = Death_lead14/cases_07da,
              rate_deathlead15 = Death_lead15/cases_07da,
              rate_deathlead16 = Death_lead16/cases_07da,
              rate_deathlead17 = Death_lead17/cases_07da,
              rate_deathlead18 = Death_lead18/cases_07da,
              rate_deathlead19 = Death_lead19/cases_07da,
              rate_deathlead20 = Death_lead20/cases_07da,
              rate_deathlead21 = Death_lead21/cases_07da,
              rate_deathlead22 = Death_lead22/cases_07da,
              rate_deathlead23 = Death_lead23/cases_07da,
              rate_deathlead24 = Death_lead24/cases_07da,
              rate_deathlead25 = Death_lead25/cases_07da)
}

#apply CFR function to each VOC dataset
WT_rates <- CalcDeathRate(WT_deaths) %>%
  mutate(VOC = rep("D614G"))

Beta_rates <- CalcDeathRate(Beta_deaths) %>%
  mutate(VOC = rep("Beta"))

Delta_rates <- CalcDeathRate(Delta_deaths) %>%
  mutate(VOC = rep("Delta"))

Omi_rates <- CalcDeathRate(Omi_deaths) %>%
  mutate(VOC = rep("Omicron"))

#aggregate results for each VOC-dominated wave
VOC_rates <- rbind(WT_rates, Beta_rates, Delta_rates, Omi_rates) %>%
  gather(., key = LagTime, value = DeathRate, -agecat, -Date, -VOC) %>%
  mutate(LagTime = parse_number(LagTime),
         VOC = factor(VOC, levels=c('D614G','Beta','Delta','Omicron')))

#correct -Inf values
VOC_rates$DeathRate[VOC_rates$DeathRate < 0.0001] <- 0.0001

#for each lag time, age cat, and VOC, calculate Mean of Log10-transformed CFR
VOC_rates <- VOC_rates %>%
  group_by(VOC, agecat, LagTime) %>%
  filter(!is.na(DeathRate)) %>%
  summarise(AvgDeathRate = geoMean(DeathRate, na.rm = T))

VOC_corr_max_CasetoDeath <- left_join(VOC_corr_max_CasetoDeath, VOC_rates, by = c("VOC", "agecat", "LagTime"))

#plot CFR by lag time, for each age strata and VOC
P5 <- ggplot(VOC_rates) +
  facet_wrap(~agecat) + theme_light() +
  geom_line(aes(x = LagTime, y = log10(AvgDeathRate), color = VOC), size = 0.8) + 
  geom_point(data=VOC_corr_max_CasetoDeath, size=2.5,
             aes(x = LagTime, y = log10(AvgDeathRate), color = VOC)) + 
  scale_color_manual(values = c("turquoise3", "royalblue4", "darkorange2", "red4"))+ 
  xlab("Lag Interval (days)") +
  ylab("Geomtric Mean CFR (log10 scale)") + 
  scale_y_continuous(breaks = c(-3, -2,-1), labels = c("0.1%", "1%", "10%")) +
  ggtitle("B. Case Fatality Ratio by Lag Interval")
P5 <- P5 + theme(plot.title = element_text(size = 12))
P5

VOC_CFR_summary <- VOCcases_age %>%
  left_join(., VOC_corr_max_CasetoDeath, by = c("VOC", "agecat")) %>%
  mutate(Weighted_CFR = Prop * AvgDeathRate)

VOC_CFR_summary %>%
  group_by(VOC) %>%
  summarise(CFR = sum(Weighted_CFR))

#CHR and CFR aggregate plot
grid.arrange(P4, P5, nrow = 1, widths = c(1, 1.3))

###################################################
# Hospital Admissions to Deaths
###################################################

#Using same death lead variables, now calculate distance correlation of COVID-19
# hospital admissions (daily raw count) with COVID-19 deaths X days in the future
# , where x is 1-25 days, by age category
CreateCorr3 <- function(x) {
  x %>%
    dplyr::group_by(agecat) %>%
    summarise(corr_Deathlead1 = dcor(HospCount[!is.na(Death_lead1)], Death_lead1[!is.na(Death_lead1)]),
              corr_Deathlead2 = dcor(HospCount[!is.na(Death_lead2)], Death_lead2[!is.na(Death_lead2)]),
              corr_Deathlead3 = dcor(HospCount[!is.na(Death_lead3)], Death_lead3[!is.na(Death_lead3)]),
              corr_Deathlead4 = dcor(HospCount[!is.na(Death_lead4)], Death_lead4[!is.na(Death_lead4)]),
              corr_Deathlead5 = dcor(HospCount[!is.na(Death_lead5)], Death_lead5[!is.na(Death_lead5)]),
              corr_Deathlead6 = dcor(HospCount[!is.na(Death_lead6)], Death_lead6[!is.na(Death_lead6)]),
              corr_Deathlead7 = dcor(HospCount[!is.na(Death_lead7)], Death_lead7[!is.na(Death_lead7)]),
              corr_Deathlead8 = dcor(HospCount[!is.na(Death_lead8)], Death_lead8[!is.na(Death_lead8)]),
              corr_Deathlead9 = dcor(HospCount[!is.na(Death_lead9)], Death_lead9[!is.na(Death_lead9)]),
              corr_Deathlead10 = dcor(HospCount[!is.na(Death_lead10)], Death_lead10[!is.na(Death_lead10)]),
              corr_Deathlead11 = dcor(HospCount[!is.na(Death_lead11)], Death_lead11[!is.na(Death_lead11)]),
              corr_Deathlead12 = dcor(HospCount[!is.na(Death_lead12)], Death_lead12[!is.na(Death_lead12)]),
              corr_Deathlead13 = dcor(HospCount[!is.na(Death_lead13)], Death_lead13[!is.na(Death_lead13)]),
              corr_Deathlead14 = dcor(HospCount[!is.na(Death_lead14)], Death_lead14[!is.na(Death_lead14)]),
              corr_Deathlead15 = dcor(HospCount[!is.na(Death_lead15)], Death_lead15[!is.na(Death_lead15)]),
              corr_Deathlead16 = dcor(HospCount[!is.na(Death_lead16)], Death_lead16[!is.na(Death_lead16)]),
              corr_Deathlead17 = dcor(HospCount[!is.na(Death_lead17)], Death_lead17[!is.na(Death_lead17)]),
              corr_Deathlead18 = dcor(HospCount[!is.na(Death_lead18)], Death_lead18[!is.na(Death_lead18)]),
              corr_Deathlead19 = dcor(HospCount[!is.na(Death_lead19)], Death_lead19[!is.na(Death_lead19)]),
              corr_Deathlead20 = dcor(HospCount[!is.na(Death_lead20)], Death_lead20[!is.na(Death_lead20)]),
              corr_Deathlead21 = dcor(HospCount[!is.na(Death_lead21)], Death_lead21[!is.na(Death_lead21)]),
              corr_Deathlead22 = dcor(HospCount[!is.na(Death_lead22)], Death_lead22[!is.na(Death_lead22)]),
              corr_Deathlead23 = dcor(HospCount[!is.na(Death_lead23)], Death_lead23[!is.na(Death_lead23)]),
              corr_Deathlead24 = dcor(HospCount[!is.na(Death_lead24)], Death_lead24[!is.na(Death_lead24)]),
              corr_Deathlead25 = dcor(HospCount[!is.na(Death_lead25)], Death_lead25[!is.na(Death_lead25)]))
}

#calculate correlation for each VOC at each lag time
WT_corr <- CreateCorr3(WT_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("D614G"))  

Beta_corr <- CreateCorr3(Beta_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Beta"))

Delta_corr <- CreateCorr3(Delta_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Delta"))

Omi_corr <- CreateCorr3(Omi_deaths) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Omicron"))

#aggregate results for each VOC-dominated wave
VOC_corr_HosptoDeath <- rbind(WT_corr, Beta_corr, Delta_corr, Omi_corr)
VOC_corr_HosptoDeath$VOC <- factor(VOC_corr_HosptoDeath$VOC, levels=c('D614G','Beta','Delta','Omicron'))

#determine optimal lag time for each VOC, age category strata
#by looking at maximum value of distance correlation within each group
VOC_corr_max_HosptoDeath <- VOC_corr_HosptoDeath %>%
  group_by(agecat, VOC) %>%
  summarise(LagTime = LagTime[which.max(Corr)], Corr = Corr[which.max(Corr)]) %>%
  mutate(VOC_num = ifelse(VOC == "D614G", 1,
                          ifelse("Beta",2,
                                 ifelse(VOC == "Delta",3,
                                        ifelse(VOC == "Omicron", 4, NA)))))

#subset of correlation data for which dCor >= 0.90 (for visualization)
VOC_corr_90_HosptoDeath <- VOC_corr_HosptoDeath %>% filter(Corr >= .90)

#heatmap of distance correlations by lag time for Hosp-Deaths
P2 <- ggplot(VOC_corr_HosptoDeath, aes(x = LagTime, y = agecat)) +
  geom_tile(aes(fill = Corr)) + 
  geom_point(data=VOC_corr_max_HosptoDeath, size=1, fill=NA, colour="black",
             aes(x = LagTime, y = agecat)) +
  geom_line(data = VOC_corr_90_HosptoDeath, aes(x = LagTime, y = agecat), colour = "black", size = 0.4, linetype = "dotted") +
  scale_fill_viridis(option = "mako") +
  xlab("Lag Interval (Days)") + 
  ylab("Age Category") + 
  facet_grid(VOC ~ .)  + theme_classic() + labs(title = "B. Hospitalizations to Deaths", fill = expression(paste("Distance \nCorrelation")))
P2 <- P2 + theme(plot.title = element_text(size = 12), legend.position = 'bottom')

#create in-hospital CFR as function of lag time
CalcDeathRate2 <- function(x) {
  x %>%
    dplyr::group_by(agecat, Date) %>%
    summarise(rate_deathlead1 = Death_lead1/HospCount,
              rate_deathlead2 = Death_lead2/HospCount,
              rate_deathlead3 = Death_lead3/HospCount,
              rate_deathlead4 = Death_lead4/HospCount,
              rate_deathlead5 = Death_lead5/HospCount,
              rate_deathlead6 = Death_lead6/HospCount,
              rate_deathlead7 = Death_lead7/HospCount,
              rate_deathlead8 = Death_lead8/HospCount,
              rate_deathlead9 = Death_lead9/HospCount,
              rate_deathlead10 = Death_lead10/HospCount,
              rate_deathlead11 = Death_lead11/HospCount,
              rate_deathlead12 = Death_lead12/HospCount,
              rate_deathlead13 = Death_lead13/HospCount,
              rate_deathlead14 = Death_lead14/HospCount,
              rate_deathlead15 = Death_lead15/HospCount,
              rate_deathlead16 = Death_lead16/HospCount,
              rate_deathlead17 = Death_lead17/HospCount,
              rate_deathlead18 = Death_lead18/HospCount,
              rate_deathlead19 = Death_lead19/HospCount,
              rate_deathlead20 = Death_lead20/HospCount,
              rate_deathlead21 = Death_lead21/HospCount,
              rate_deathlead22 = Death_lead22/HospCount,
              rate_deathlead23 = Death_lead23/HospCount,
              rate_deathlead24 = Death_lead24/HospCount,
              rate_deathlead25 = Death_lead25/HospCount)
}

#apply within-hosp CFR function to each VOC dataset
WT_rates <- CalcDeathRate2(WT_deaths) %>%
  mutate(VOC = rep("D614G"))

Beta_rates <- CalcDeathRate2(Beta_deaths) %>%
  mutate(VOC = rep("Beta"))

Delta_rates <- CalcDeathRate2(Delta_deaths) %>%
  mutate(VOC = rep("Delta"))

Omi_rates <- CalcDeathRate2(Omi_deaths) %>%
  mutate(VOC = rep("Omicron"))

#aggregate results for each VOC-dominated wave
VOC_rates <- rbind(WT_rates, Beta_rates, Delta_rates, Omi_rates) %>%
  gather(., key = LagTime, value = HospDeathRate, -agecat, -Date, -VOC) %>%
  mutate(LagTime = parse_number(LagTime),
         VOC = factor(VOC, levels=c('D614G','Beta','Delta','Omicron')))

#correct -Inf values
VOC_rates$HospDeathRate[VOC_rates$HospDeathRate < 0.001] <- 0.001

#for each lag time, age cat, and VOC, calculate Mean of Log10-transformed CFR
VOC_rates <- VOC_rates %>%
  group_by(VOC, agecat, LagTime) %>%
  summarise(AvgHospDeathRate = geoMean(HospDeathRate, na.rm = T))

VOC_corr_max_HosptoDeath <- left_join(VOC_corr_max_HosptoDeath, VOC_rates, by = c("VOC", "agecat", "LagTime"))

#plot within-hosp CFR by lag time, for each age strata and VOC
P6 <- ggplot(VOC_rates) +
  facet_wrap(~agecat) + theme_light() +
  geom_line(aes(x = LagTime, y = AvgHospDeathRate, color = VOC), size = 0.8) + 
  geom_point(data=VOC_corr_max_HosptoDeath, size=2.5,
             aes(x = LagTime, y = AvgHospDeathRate, color = VOC)) + 
  scale_color_manual(values = c("turquoise3", "royalblue4", "darkorange2", "red4"))+ 
  xlab("Lag Interval (days)") +
  ylab("Mean of Log10 CFR") + 
  #scale_y_continuous(breaks = c(-3, -2,-1), labels = c("0.1%", "1%", "10%")) +
  ggtitle("B. Case Fatality Ratio by Lag Interval")
P6 <- P6 + theme(plot.title = element_text(size = 12))
P6

#calculate age group composition of hospitalizations for each VOC
WT_byage2 <- WT %>%
  group_by(agecat) %>%
  summarise(N = sum(HospCount), Prop = N/106568) %>%
  mutate(VOC = rep("D614G"))

Beta_byage2 <- BETA %>%
  group_by(agecat) %>%
  summarise(N = sum(HospCount), Prop = N/137397) %>%
  mutate(VOC = rep("Beta"))

Delta_byage2 <- DELTA %>%
  group_by(agecat) %>%
  summarise(N = sum(HospCount), Prop = N/172446) %>%
  mutate(VOC = rep("Delta"))

Omi_byage2 <- OMI %>%
  group_by(agecat) %>%
  summarise(N = sum(HospCount), Prop = N/60250) %>%
  mutate(VOC = rep("Omicron"))

VOChosp_age <- rbind(WT_byage2, Beta_byage2, Delta_byage2, Omi_byage2)

VOC_CFRHosp_summary <- VOChosp_age %>%
  left_join(., VOC_corr_max_HosptoDeath, by = c("VOC", "agecat")) %>%
  mutate(Weighted_CFR_Hosp = Prop * AvgHospDeathRate)

VOC_CFRHosp_summary %>%
  group_by(VOC) %>%
  summarise(CFRHosp = sum(Weighted_CFR_Hosp))