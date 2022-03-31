# Plots - sensitivity analysis

library(tidyverse)

##### read in results data

OmiRes <- read_csv("results/ResBind_omi.csv")
DeltaRes <- read_csv("results/ResBind_delta.csv")

OmiCFRs <- read_csv("results/L10CFRbyRange_omi.csv")
DeltaCFRs <- read_csv("results/L10CFRbyRange_delta.csv")

##### Supplementary 1a

OmiRes <- OmiRes |> mutate(VOC="Omicron")
DeltaRes <- DeltaRes |> mutate(VOC="Delta")

s1a <- bind_rows(OmiRes, DeltaRes) |> 
  ggplot(aes(x=range, y=BestLag, group=agecat, color=agecat)) +
  geom_line() +
  theme_classic() +
  labs(x="Days of Available Data After Wave Start",
       y="Optimal Lag (Days)",
       color="Age") +
  scale_color_brewer(palette = "Paired") +
  facet_wrap(~ VOC, scales = "free_x", nrow = 2) +
  ggtitle("A. Lag Intervals Calculated with Limited Data")
s1a
ggsave("figures/S1a.pdf")

##### Supplementary 1b

OmiCFRs <- OmiCFRs |> 
  pivot_longer(!range, values_to = "cfr", 
               names_to = "age") |> 
  mutate(VOC="Omicron")

DeltaCFRs <- DeltaCFRs |> 
  pivot_longer(!range, values_to = "cfr", 
               names_to = "age") |> 
  mutate(VOC="Delta")

s1b <- bind_rows(OmiCFRs, DeltaCFRs) |> 
  ggplot(aes(x=range, y=cfr, 
             group=age, color=age)) +
  geom_line() +
  theme_classic() +
  labs(x="Days of Available Data After Wave Start", 
       y="Geometric Mean CFR (log10 scale)",
       color="Age") +
  scale_color_brewer(palette = "Paired") +
  facet_wrap(~ VOC, scales = "free_x", nrow = 2) +
  ggtitle("B. CFRs Calculated with Limited Data")
s1b
ggsave("figures/S1b.pdf")




