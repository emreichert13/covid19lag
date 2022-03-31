# sensitivity analysis with limited data - omicron

library(tidyverse)
library(lubridate)
library(viridis)
library(zoo)
library(energy)

##### load data
cases <- read_csv("data/cases_by_age/E_Reichert_01022022.csv")

deaths <- readxl::read_xlsx("data/hosp_by_age/DATCOV.01022022.xlsx",
                            sheet = "InHospitalDeaths")

##### cleaning

## cases
cases <- cases |> 
  
  # parsing dates
  mutate(DateCollect = parse_date_time(speccollectiondate, orders = c('dmy', 'mdy')),
         DateReceive = parse_date_time(specreceiveddate, orders = c('dmy', 'mdy')),
         DateReport = parse_date_time(specreportdate, orders = c('dmy', 'mdy'))) |> 
  
  mutate(DateCollect = as_date(DateCollect),
         DateReceive = as_date(DateReceive),
         DateReport = as_date(DateReport)) |> 
  
  # erroneous ages present
  rename(Age = patage) |> 
  filter(Age <= 110 & Age >= 0) |> 
  
  # no cases in SA prior to 2020-03-01
  filter(DateCollect >= "2020-03-01" & DateCollect <= "2022-02-01") |> 
  
  # did not do any filtering by date report, date receive
  
  # create age categories
  mutate(agecat = ifelse(Age <= 17, "0-17",
                         ifelse(Age > 17 & Age <= 29, "18-29",
                                ifelse(Age > 29 & Age <= 39, "30-39",
                                       ifelse(Age > 39 & Age <= 49, "40-49",
                                              ifelse(Age > 49 & Age <= 64, "50-64",
                                                     ifelse(Age > 64 & Age <= 74, "65-74",
                                                            ifelse(Age > 74, "75+", NA)))))))) |> 
  
  mutate(province = as_factor(province),
         district = as_factor(district),
         subdistrict = as_factor(subdistrict),
         agecat = as_factor(agecat)) |>
  
  mutate(agecat = factor(agecat,
                         levels = c("0-17", "18-29", "30-39", "40-49", 
                                    "50-64", "65-74", "75+"))) |> 
  
  # including DateReport for sensitivity analysis
  select(DateCollect, DateReport, Age, agecat, province, district, subdistrict) |>
  
  arrange(DateCollect) |> 
  
  rename(Date = DateCollect)

## deaths
deaths <- deaths |> 
  rename(Date = `Death Date`) |> 
  mutate(Date = as_date(Date)) |> 
  
  filter(Age <= 110) |> 
  
  # no cases in SA prior to 2020-03-01
  filter(Date >= "2020-03-01" & Date <= "2022-02-01") |> 
  
  rename(province = Province) |> 
  mutate(province = as_factor(province)) |> 
  
  mutate(agecat = ifelse(Age <= 17, "0-17",
                         ifelse(Age > 17 & Age <= 29, "18-29",
                                ifelse(Age > 29 & Age <= 39, "30-39",
                                       ifelse(Age > 39 & Age <= 49, "40-49",
                                              ifelse(Age > 49 & Age <= 64, "50-64",
                                                     ifelse(Age > 64 & Age <= 74, "65-74",
                                                            ifelse(Age > 74, "75+", NA)))))))) |>
  
  mutate(agecat = factor(agecat,
                         levels = c("0-17", "18-29", "30-39", "40-49", 
                                    "50-64", "65-74", "75+"))) |> 
  
  arrange(Date)


##### wave bounds

omicron_t0 <- ymd("2021-11-07")
omicron_tt <- ymd("2022-01-28")


##### creating limited data lists

# ymd("2022-01-28") - ymd("2021-11-07")
# full analysis utilizes 82 days

# here we see how accuracy increases as time progresses
# from 7d to 82d of data, in 5d increments

# seq(7,82,5)

DataList <- list()

for (i in seq(7,82,5)) {
  
  nam <- paste("CaseDeath", i, sep = "_")
  
  c <- cases |> 
    
    filter(Date >= omicron_t0-3) |> # 3d before so that rolling ave starts at t0
    filter(DateReport >= omicron_t0-3 & DateReport <= omicron_t0+i) |> 
    
    group_by(Date, agecat) |>
    summarise(CaseCount = n())
  
  d <- deaths |> 
    
    filter(Date >= omicron_t0-3 & Date <= omicron_t0+i) |> # 3d before so that rolling ave starts at t0
    # no reporting date for deaths  
    group_by(Date, agecat) |>
    summarise(DeathCount = sum(Deaths))
  
  DataList[[nam]] <- full_join(c, d, by = c("Date", "agecat")) |> 
    mutate(CaseCount = replace_na(CaseCount, 0),
           DeathCount = replace_na(DeathCount, 0))
  
}

# dim(DataList$CaseDeath_7)
# dim(DataList$CaseDeath_82)

# For each restricted data set, we have 4 cols: Date, agecat, CaseCount, DeathCount with raw data
# This simulates what was known at each time point
# tibbles range from 77 to 602 rows

rm(c,d,i,nam)

agecats <- levels(DataList$CaseDeath_7$agecat)
dataranges <- seq(7,82,5)
list <- list()

##### for each data range, we generate daily aggregate case and death counts, age stratified

for (i in seq(7,82,5)) {
  
  rangename <- paste("CaseDeath", i, sep = "_")
  
  for (j in agecats) {
    
    d <- DataList[[rangename]] |> ungroup() |> filter(agecat==j) |>
      group_by(Date) |>
      summarise(AggCaseCount = sum(CaseCount), AggDeathCount = sum(DeathCount)) |>
      mutate(AggCaseCountR7 = rollmean(AggCaseCount, k=7, fill = NA)) #|>
    
    # filter(Date >= t0 & Date <= tt)
    
    age <- paste("Age", j, sep = "_")
    age <- gsub("-","_", age)
    age <- gsub("\\+","_", age)
    
    range <- paste("Range",i,sep = "")
    
    list[[range]][[age]] <- d
  }
  
}

rm(d,i,j,range,rangename,age)

# now we have a list of lists - a list for each data range containing a tibble for each agecat

ranges <- names(list)

#####  create death lead variables

leads <- list()

# load functions CreateLeadsD and DropMissingCases

CreateLeadsD <- function(x) {
  x |>
    mutate(death_lead1 = lead(AggDeathCount, 1)) |>
    mutate(death_lead2 = lead(AggDeathCount, 2)) |>
    mutate(death_lead3 = lead(AggDeathCount, 3)) |>
    mutate(death_lead4 = lead(AggDeathCount, 4)) |>
    mutate(death_lead5 = lead(AggDeathCount, 5)) |>
    mutate(death_lead6 = lead(AggDeathCount, 6)) |>
    mutate(death_lead7 = lead(AggDeathCount, 7)) |>
    mutate(death_lead8 = lead(AggDeathCount, 8)) |>
    mutate(death_lead9 = lead(AggDeathCount, 9)) |>
    mutate(death_lead10 = lead(AggDeathCount, 10)) |>
    mutate(death_lead11 = lead(AggDeathCount, 11)) |>
    mutate(death_lead12 = lead(AggDeathCount, 12)) |>
    mutate(death_lead13 = lead(AggDeathCount, 13)) |>
    mutate(death_lead14 = lead(AggDeathCount, 14)) |> 
    mutate(death_lead15 = lead(AggDeathCount, 15)) |> 
    mutate(death_lead16 = lead(AggDeathCount, 16)) |> 
    mutate(death_lead17 = lead(AggDeathCount, 17)) |> 
    mutate(death_lead18 = lead(AggDeathCount, 18)) |> 
    mutate(death_lead19 = lead(AggDeathCount, 19)) |> 
    mutate(death_lead20 = lead(AggDeathCount, 20)) |> 
    mutate(death_lead21 = lead(AggDeathCount, 21)) |> 
    mutate(death_lead22 = lead(AggDeathCount, 22)) |> 
    mutate(death_lead23 = lead(AggDeathCount, 23)) |> 
    mutate(death_lead24 = lead(AggDeathCount, 24)) |> 
    mutate(death_lead25 = lead(AggDeathCount, 25))
}

# drop rows with missing case data - due to rolling ave

DropMissingCases <- function(x) {
  x |> drop_na(AggCaseCountR7)
}

for (k in ranges) {
  leads[[k]] <- lapply(list[[k]], CreateLeadsD)
}

# drop rows with missing cases
leads_v2 <- list()

for (k in ranges) {
  leads_v2[[k]] <- lapply(leads[[k]], DropMissingCases)
}

rm(k)

# resulting format
# leads_v2$Range7$Age_0_17$death_lead1
# list$list$tib$col

##### dcor main calculations

leadtimes <- leads_v2$Range7$Age_0_17 |> ungroup() |> select(starts_with("death_lead")) |> colnames()
OutputListFULL <- list()
OutputListPsFULL <- list()

# for each range, agecat, and lead time, we calculate dcor for cases and deaths

for (x in ranges) {
  
  d <- leads_v2[[x]]
  
  for (y in agecats) {
    
    age <- paste("Age", y, sep = "_")
    age <- gsub("-","_", age)
    age <- gsub("\\+","_", age)
    
    for (z in leadtimes) {
      d2 <- d[[age]]
      d_cases <- d2$AggCaseCountR7 
      d_deaths <- d2[[z]]
      d_cases_deaths <- cbind(d_cases,d_deaths) |> as_tibble() |> drop_na()
      OutputListFULL[[x]][[y]][z] <- dcor(d_cases_deaths$d_cases, d_cases_deaths$d_deaths)
      OutputListPsFULL[[x]][[y]][z] <- dcorT.test(d_cases_deaths$d_cases, d_cases_deaths$d_deaths)$p.value
    }
    
  }
  
}

rm(x,y,z,age,d,d2,d_cases,d_deaths,d_cases_deaths)

##### clean up results function

FinalizeRes <- function(x,y) {
  xtemp <- x |> 
    mutate(days = row_number()) |>
    pivot_longer(!days, names_to = "agecat", values_to = "DistCor") |>
    mutate(agecat = gsub("CaseDeath_", "", agecat))
  
  ytemp <- y |> 
    mutate(days = row_number()) |>
    pivot_longer(!days, names_to = "agecat", values_to = "ttestPval") |>
    mutate(agecat = gsub("CaseDeath_", "", agecat)) 
  
  long <- xtemp |> 
    left_join(ytemp, by=c("days", "agecat")) |> 
    select(days, agecat, DistCor, ttestPval)
  
  long |> 
    group_by(agecat) |> 
    summarise(
      BestLag = days[which.max(DistCor)],
      DistCor = DistCor[which.max(DistCor)],
      Pval = ttestPval[which.max(DistCor)]
    )
}

##### compiling results

Res <- list()

for (i in seq(7,82,5)) {
  
  range <- paste("Range",i,sep="")
  rangep <- paste("Range",i,"p",sep="")
  
  Res[[range]] <- OutputListFULL[[range]] |> as_tibble()
  Res[[rangep]] <- OutputListPsFULL[[range]] |> as_tibble()
  
}

ResFin <- list()

for (i in seq(7,82,5)) {
  
  range <- paste("Range",i,sep="")
  rangep <- paste("Range",i,"p",sep="")
  
  ResFin[[range]] <- FinalizeRes(Res[[range]],Res[[rangep]])
  ResFin[[range]] <- ResFin[[range]] |> mutate(range=i)
  
}

rm(range,rangep,i)

ResBind <- do.call(rbind,ResFin) 

##### export results

# ResBind |>
#   write_csv("results/ResBind_omi.csv")
# last export = 2022-03-06


##### calculating CFRs

# using the same process as above, calculate CFR for each age group at each interval
# seq(7,82,5) (7 12 17 22 27 32 37 42 47 52 57 62 67 72 77 82)

# DataList contains raw data
# ResBind contains age-strat optimal lag times for each data range

# names(DataList)
# names(ResBind)

# calculate log10 CFRs - geometric mean
log10CFRlist <- list()

for (i in seq(7,82,5)) {
  
  range <- paste("CaseDeath",i,sep="_")
  
  for (j in agecats) {
    
    LEAD <- ResBind$BestLag[which(ResBind$agecat==j & ResBind$range==i)]
    
    cfr <- DataList[[range]] |> 
      ungroup() |> 
      filter(agecat==j) |> 
      mutate(CaseCountR7 = rollmean(CaseCount, k=7, fill = NA)) |> 
      mutate(DeathLead = lead(DeathCount,LEAD)) |>
      mutate(CFR = DeathLead/CaseCountR7) |>
      mutate(CFR=ifelse(CFR<0.0001,0.0001,CFR)) |> 
      summarise(CFR=10^(mean(log10(CFR),na.rm=T))) #geometric mean
    
    log10CFRlist[[range]][j] <- cfr
  }
  
}

rm(i,j,range,LEAD,cfr)

L10CFRbyRange <- lapply(log10CFRlist, as_tibble)
L10CFRbyRange <- do.call(rbind,L10CFRbyRange)
L10CFRbyRange$range <- seq(7,82,5)

##### export data

# L10CFRbyRange |>
#   select(range, `0-17`:`75+`) |>
#   write_csv("results/L10CFRbyRange_omi.csv")
# last export = 2022-03-06





