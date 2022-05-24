## SENSITIVITY ANALYSIS PT II ##


#calculate CHR, CFR, and within-hospital CFR using 'conventional' method
# i.e., simply dividing cumulative number of outcomes by number of cases from t0-t1
# doesn't account for lag interval, assumes all outcomes have accumulated
wt_totals <- WT %>%
  group_by(agecat) %>%
  summarise(Cases = sum(CaseCount),
            Hosp = sum(HospCount),
            Deaths = sum(DeathCount),
            CHR = Hosp/Cases,
            CFR = Deaths/Cases,
            CFRHosp = Deaths/Hosp)

beta_totals <- BETA %>%
  group_by(agecat) %>%
  summarise(Cases = sum(CaseCount),
            Hosp = sum(HospCount),
            Deaths = sum(DeathCount),
            CHR = Hosp/Cases,
            CFR = Deaths/Cases,
            CFRHosp = Deaths/Hosp)

delta_totals <- DELTA %>%
  group_by(agecat) %>%
  summarise(Cases = sum(CaseCount),
            Hosp = sum(HospCount),
            Deaths = sum(DeathCount),
            CHR = Hosp/Cases,
            CFR = Deaths/Cases,
            CFRHosp = Deaths/Hosp)

omi_totals <- OMI %>%
  group_by(agecat) %>%
  summarise(Cases = sum(CaseCount),
            Hosp = sum(HospCount),
            Deaths = sum(DeathCount),
            CHR = Hosp/Cases,
            CFR = Deaths/Cases,
            CFRHosp = Deaths/Hosp)

#read in results from sensitivity analysis using Proposed Method
omicron_sens <- read_csv("Documents/2021 CCDD research/Omicronsens.csv") %>%
  select(range = range, agecat = age, CFR = log10CFR) %>%
  mutate(Method = "Proposed Method")

#calculate using conventional method within-epidemic for Omicron VOC
omicron_t0 <- ymd("2021-11-07")
dataranges <- seq(5,80,5)
OMI_res <- data.frame()
for (i in seq(5,80,5)) {
  OMI_i <- OMI %>%
    filter(Date <= omicron_t0+i) %>%
    mutate(range = i, Method = "Conventional Method")
  OMI_res <- rbind(OMI_res, OMI_i)
}

OMI_res <- OMI_res %>%
  group_by(agecat, range, Method) %>%
  summarise(CFR = sum(DeathCount)/sum(CaseCount))

#Compare the two methods
OMI_sensitivity <- rbind(OMI_res, omicron_sens)

p1 <- ggplot() +
  geom_line(aes(x=range, y=(CFR), group=agecat, color=agecat), data = OMI_sensitivity, size = 0.7) +
  theme_classic() +
  labs(x="Days of Available Data after t0", y="CFR",color="Age") +
  scale_color_brewer(palette = "Paired") +
  facet_wrap(~Method) +
  geom_hline(aes(yintercept = (CFR), group = agecat, color = agecat), data = omi_totals, linetype = "dashed") +
  ggtitle("Omicron VOC")

#read in results from sensitivity analysis using Proposed Method
delta_sens <- read_csv("Documents/2021 CCDD research/L10CFRbyRange_delta.csv") %>%
  gather(., key = "agecat", value = "CFR", -range) %>%
  mutate(Method = "Proposed Method") 

#calculate using conventional method within-epidemic for Delta VOC
delta_t0 <- ymd("2021-05-02")
dataranges <- seq(5,175,5)
DELTA_res <- data.frame()
for (i in seq(5,175,5)) {
  DELTA_i <- DELTA %>%
    filter(Date <= delta_t0+i) %>%
    mutate(range = i, Method = "Conventional Method")
  DELTA_res <- rbind(DELTA_res, DELTA_i)
}

DELTA_res <- DELTA_res %>%
  group_by(agecat, range, Method) %>%
  summarise(CFR = sum(DeathCount)/sum(CaseCount))

#Compare the two methods
DELTA_sensitivity <- rbind(DELTA_res, delta_sens)

p2 <- ggplot() +
  geom_line(aes(x=range, y=(CFR), group=agecat, color=agecat), data = DELTA_sensitivity, size = 0.7) +
  theme_classic() +
  labs(x="Days of Available Data after t0", y="CFR",color="Age") +
  scale_color_brewer(palette = "Paired") +
  facet_wrap(~Method) +
  geom_hline(aes(yintercept = (CFR), group = agecat, color = agecat), data = delta_totals, linetype = "dashed") +
  ggtitle("Delta VOC")

pdf("Documents/2021 CCDD research/FigureS2.pdf", width = 9, height = 7)
grid.arrange(p1, p2, nrow = 2)
dev.off()

