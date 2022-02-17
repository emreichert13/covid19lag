#E Reichert 2022
# COVID-19 Case to Hospitalization Lag Analysis using NICD Surveillance Data

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

# Here, we first run a script that reads in and cleans surveillance data (CaseAdmitAvg) 
# containing daily COVID-19 case counts, hospital admissions with COVID-19, and within-hospital
# COVID-19 attributable deaths (3 separate files), disaggregated by age category. 
source("Documents/2021 CCDD research/NICDdatacleaning.R")

#create 7 day avg for COVID-19 cases, hospital admissions, and deaths
# Note: 7d avg hospital and death variables will only be used for epi curve data visualization
CaseAdmitAvg <- CaseAdmit %>%
  dplyr::arrange(desc(agecat)) %>% 
  dplyr::group_by(agecat) %>% 
  dplyr::mutate(cases_07da = zoo::rollmean(CaseCount, k = 7, fill = NA),
                hosp_07da = zoo::rollmean(HospCount, k = 7, fill = NA),
                deaths_07da = zoo::rollmean(DeathCount, k = 7, fill = NA)) %>% 
  dplyr::ungroup()

#visualize selected indicator over time, by agecat
CaseAdmitAvg %>%
  ggplot() + 
  geom_line(aes(x = Date, y = log10(cases_07da)), col = "blue") + 
  geom_line(aes(x = Date, y = log10(hosp_07da)), col = "red") + 
  facet_wrap(~agecat)


#Set t0 and t1 for each VOC-dominated epidemic wave
# note - dates specific to South Africa here
WT_t0 <- ymd("2020-04-05")
beta_t0 <- ymd("2020-11-01")
delta_t0 <- ymd("2021-05-02")
omicron_t0 <- ymd("2021-11-07")

WT_t1 <- ymd("2020-10-17")
beta_t1 <- ymd("2021-03-20")
delta_t1 <- ymd("2021-10-23")
omicron_t1 <- ymd("2022-01-28")

WT <- CaseAdmitAvg %>%
  filter(Date >= WT_t0 & Date <= WT_t1)

BETA <- CaseAdmitAvg %>%
  filter(Date >= beta_t0 & Date <= beta_t1)

DELTA <- CaseAdmitAvg %>%
  filter(Date >= delta_t0 & Date <= delta_t1)

OMI <- CaseAdmitAvg %>%
  filter(Date >= omicron_t0 & Date <= omicron_t1)

# Function to create lead time variables for hospitalizations
# note, each hosp_lead[x] variable contains hospital admission counts for x days in the future
CreateLeads <- function(x) {
  x %>%
    mutate(hosp_lead1 = lead(HospCount, 1)) %>%
    mutate(hosp_lead2 = lead(HospCount, 2)) %>%
    mutate(hosp_lead3 = lead(HospCount, 3)) %>%
    mutate(hosp_lead4 = lead(HospCount, 4)) %>%
    mutate(hosp_lead5 = lead(HospCount, 5)) %>%
    mutate(hosp_lead6 = lead(HospCount, 6)) %>%
    mutate(hosp_lead7 = lead(HospCount, 7)) %>%
    mutate(hosp_lead8 = lead(HospCount, 8)) %>%
    mutate(hosp_lead9 = lead(HospCount, 9)) %>%
    mutate(hosp_lead10 = lead(HospCount, 10)) %>%
    mutate(hosp_lead11 = lead(HospCount, 11)) %>%
    mutate(hosp_lead12 = lead(HospCount, 12)) %>%
    mutate(hosp_lead13 = lead(HospCount, 13)) %>%
    mutate(hosp_lead14 = lead(HospCount, 14)) %>%
    mutate(hosp_lead15 = lead(HospCount, 15)) %>%
    mutate(hosp_lead16 = lead(HospCount, 16)) %>%
    mutate(hosp_lead17 = lead(HospCount, 17)) %>%
    mutate(hosp_lead18 = lead(HospCount, 18)) %>%
    mutate(hosp_lead19 = lead(HospCount, 19)) %>%
    mutate(hosp_lead20 = lead(HospCount, 20))
    
}

#Now apply hospital leads function to each VOC dataset
WT_lags <- WT %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads() %>% 
  dplyr::ungroup() 

Beta_lags <- BETA %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads() %>% 
  dplyr::ungroup()

Delta_lags <- DELTA %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads() %>% 
  dplyr::ungroup()

Omi_lags <- OMI %>%
  dplyr::arrange(agecat, Date) %>% 
  dplyr::group_by(agecat) %>% 
  CreateLeads() %>% 
  dplyr::ungroup() %>%
  filter(!is.na(cases_07da))

# Create correlation function - assess distance correlation of 
# COVID-19 cases (7da average) with # of 
# COVID-19 hospitalizations x days in the future, where x is 1-20 days, by age cat

CreateCorr <- function(x) {
  x %>%
    dplyr::group_by(agecat) %>%
    summarise(corr_hosplead1 = dcor(cases_07da[!is.na(hosp_lead1)], hosp_lead1[!is.na(hosp_lead1)]),
              corr_hosplead2 = dcor(cases_07da[!is.na(hosp_lead2)], hosp_lead2[!is.na(hosp_lead2)]),
              corr_hosplead3 = dcor(cases_07da[!is.na(hosp_lead3)], hosp_lead3[!is.na(hosp_lead3)]),
              corr_hosplead4 = dcor(cases_07da[!is.na(hosp_lead4)], hosp_lead4[!is.na(hosp_lead4)]),
              corr_hosplead5 = dcor(cases_07da[!is.na(hosp_lead5)], hosp_lead5[!is.na(hosp_lead5)]),
              corr_hosplead6 = dcor(cases_07da[!is.na(hosp_lead6)], hosp_lead6[!is.na(hosp_lead6)]),
              corr_hosplead7 = dcor(cases_07da[!is.na(hosp_lead7)], hosp_lead7[!is.na(hosp_lead7)]),
              corr_hosplead8 = dcor(cases_07da[!is.na(hosp_lead8)], hosp_lead8[!is.na(hosp_lead8)]),
              corr_hosplead9 = dcor(cases_07da[!is.na(hosp_lead9)], hosp_lead9[!is.na(hosp_lead9)]),
              corr_hosplead10 = dcor(cases_07da[!is.na(hosp_lead10)], hosp_lead10[!is.na(hosp_lead10)]),
              corr_hosplead11 = dcor(cases_07da[!is.na(hosp_lead11)], hosp_lead11[!is.na(hosp_lead11)]),
              corr_hosplead12 = dcor(cases_07da[!is.na(hosp_lead12)], hosp_lead12[!is.na(hosp_lead12)]),
              corr_hosplead13 = dcor(cases_07da[!is.na(hosp_lead13)], hosp_lead13[!is.na(hosp_lead13)]),
              corr_hosplead14 = dcor(cases_07da[!is.na(hosp_lead14)], hosp_lead14[!is.na(hosp_lead14)]),
              corr_hosplead15 = dcor(cases_07da[!is.na(hosp_lead15)], hosp_lead15[!is.na(hosp_lead15)]),
              corr_hosplead16 = dcor(cases_07da[!is.na(hosp_lead16)], hosp_lead16[!is.na(hosp_lead16)]),
              corr_hosplead17 = dcor(cases_07da[!is.na(hosp_lead17)], hosp_lead17[!is.na(hosp_lead17)]),
              corr_hosplead18 = dcor(cases_07da[!is.na(hosp_lead18)], hosp_lead18[!is.na(hosp_lead18)]),
              corr_hosplead19 = dcor(cases_07da[!is.na(hosp_lead19)], hosp_lead19[!is.na(hosp_lead19)]),
              corr_hosplead20 = dcor(cases_07da[!is.na(hosp_lead20)], hosp_lead20[!is.na(hosp_lead20)]))
}

#apply correlation function for each VOC dataset
WT_corr <- CreateCorr(WT_lags) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("D614G"))

Beta_corr <- CreateCorr(Beta_lags) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Beta"))

Delta_corr <- CreateCorr(Delta_lags) %>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Delta"))

Omi_corr <- CreateCorr(Omi_lags)%>%
  gather(., key = LagTime, value = Corr, -agecat) %>%
  mutate(LagTime = parse_number(LagTime))%>%
  mutate(VOC = rep("Omicron"))

#aggregate results for each VOC-dominated wave
VOC_corr_CasetoHosp <- rbind(WT_corr, Beta_corr, Delta_corr, Omi_corr)
VOC_corr_CasetoHosp$VOC = factor(VOC_corr_CasetoHosp$VOC, levels=c('D614G','Beta','Delta','Omicron'))

#determine optimal lag time for each VOC, age category strata
#by looking at maximum value of distance correlation within each group
VOC_corr_max_CasetoHosp <- VOC_corr_CasetoHosp %>%
  group_by(agecat, VOC) %>%
  summarise(LagTime = LagTime[which.max(Corr)], Corr = Corr[which.max(Corr)]) %>%
  mutate(VOC_num = ifelse(VOC == "D614G", 1,
                          ifelse("Beta",2,
                                 ifelse(VOC == "Delta",3,
                                        ifelse(VOC == "Omicron", 4, NA)))))

#subset of correlation data for which dCor >= 0.90 (for visualization)
VOC_corr_90_CasetoHosp <- VOC_corr_CasetoHosp %>% filter(Corr >= .90)

#heatmap of distance correlations by lag time for Cases-Hospitalizations
P1 <- ggplot(VOC_corr_CasetoHosp, aes(x = LagTime, y = agecat)) +
  geom_tile(aes(fill = Corr)) + 
  geom_point(data=VOC_corr_max_CasetoHosp, size=1, fill=NA, colour="black",
            aes(x = LagTime, y = agecat)) +
  geom_line(data = VOC_corr_90_CasetoHosp, aes(x = LagTime, y = agecat), colour = "black", size = 0.4, linetype = "dotted") +  scale_fill_viridis(option = "mako") + 
  ggtitle("A. Cases to Hospitalizations") +
  theme_classic() + xlab("Lag Interval (Days)") + 
  ylab("Age Category") +
  facet_grid(VOC ~ .) + labs(fill = expression(paste("Distance \nCorrelation")))
P1 <- P1 + theme(plot.title = element_text(size = 12), legend.position = "bottom")
P1

#create CHR as function of lag time
CalcHospRate <- function(x) {
  x %>%
    dplyr::group_by(agecat, Date) %>%
    summarise(rate_hosplead1 = hosp_lead1/cases_07da,
              rate_hosplead2 = hosp_lead2/cases_07da,
              rate_hosplead3 = hosp_lead3/cases_07da,
              rate_hosplead4 = hosp_lead4/cases_07da,
              rate_hosplead5 = hosp_lead5/cases_07da,
              rate_hosplead6 = hosp_lead6/cases_07da,
              rate_hosplead7 = hosp_lead7/cases_07da,
              rate_hosplead8 = hosp_lead8/cases_07da,
              rate_hosplead9 = hosp_lead9/cases_07da,
              rate_hosplead10 = hosp_lead10/cases_07da,
              rate_hosplead11 = hosp_lead11/cases_07da,
              rate_hosplead12 = hosp_lead12/cases_07da,
              rate_hosplead13 = hosp_lead13/cases_07da,
              rate_hosplead14 = hosp_lead14/cases_07da,
              rate_hosplead15 = hosp_lead15/cases_07da,
              rate_hosplead16 = hosp_lead16/cases_07da,
              rate_hosplead17 = hosp_lead17/cases_07da,
              rate_hosplead18 = hosp_lead18/cases_07da,
              rate_hosplead19 = hosp_lead19/cases_07da,
              rate_hosplead20 = hosp_lead20/cases_07da)
}

#apply hospitalization rate function to each VOC dataset
WT_rates <- CalcHospRate(WT_lags) %>%
  mutate(VOC = rep("D614G"))

Beta_rates <- CalcHospRate(Beta_lags) %>%
  mutate(VOC = rep("Beta"))

Delta_rates <- CalcHospRate(Delta_lags) %>%
  mutate(VOC = rep("Delta"))

Omi_rates <- CalcHospRate(Omi_lags) %>%
  mutate(VOC = rep("Omicron"))

#aggregate results for each VOC-dominated wave
VOC_rates <- rbind(WT_rates, Beta_rates, Delta_rates, Omi_rates) %>%
  gather(., key = LagTime, value = HospRate, -agecat, -Date, -VOC) %>%
  mutate(LagTime = parse_number(LagTime),
         VOC = factor(VOC, levels=c('D614G','Beta','Delta','Omicron')))
#correct -Inf values
VOC_rates$HospRate[VOC_rates$HospRate < 0.01] <- 0.01

#for each lag time, age cat, and VOC, calculate Mean of Log10-transformed CHR
VOC_rates <- VOC_rates %>%
  group_by(VOC, agecat, LagTime) %>%
  summarise(AvgHospRate = mean(log10(HospRate), na.rm = T))

VOC_corr_max_CasetoHosp <- left_join(VOC_corr_max_CasetoHosp, VOC_rates, by = c("VOC", "agecat", "LagTime"))

#plot CHR by lag time, for each age strata and VOC
P4 <- ggplot(VOC_rates) +
  facet_wrap(~agecat) + theme_light() +
  geom_line(aes(x = LagTime, y = AvgHospRate, color = VOC, fill = VOC), size = 0.8) + 
  geom_point(data=VOC_corr_max_CasetoHosp, size=2.5,
             aes(x= LagTime, y = AvgHospRate, color = VOC)) + 
  scale_color_manual(values = c("turquoise3", "royalblue4", "darkorange2", "red4"))+ 
  xlab("Lag Interval (days)") +
  ylab("Mean of Log10 CHR") + scale_y_continuous(breaks = c(-1.30103, -1, -0.69897, -0.39794,-0.09691001), labels = c("5%", "10%", "20%", "40%", "80%")) +
  ggtitle("A. Case Hospitalization Ratio by Lag Interval")
P4 <- P4 + theme(legend.position = 'none', plot.title = element_text(size = 12))
P4

