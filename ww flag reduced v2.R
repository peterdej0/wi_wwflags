##########################################################################################################
#
# create wastewater data increase flags with quantile threshold
#
##########################################################################################################

# Load libraries -----------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(zoo)
library(grid)
library(BINCOR)
library(scales)
library(officer)

# Working directory --------------------------------------------------------------------------------------

setwd("L:/Beoh_Fellows/EIS Fellows/P DeJonge/WIDHS Projects/Wastewater surveillance/Datasets")

# Load most current version of wastewater data -----------------------------------------------------------

setwd("L:/COVID Wastewater/Tableau Visualizations/")
wwtp_data = read.csv("wwtp_covid_data_v3.csv") %>%
  mutate(
    WWTP = wwtp_name,
    # flow_rate_l = average_flow_rate*3.7854*1e6,
    # sars_cov2_adj_load = (geoconc*flow_rate_l)/population_served,
    # sars_cov2_adj_load_log10 = log(sars_cov2_adj_load, base=10), 
    date = as.Date(sample_collect_date))

length(unique(wwtp_data$WWTP)) # 74 sewersheds
max(wwtp_data$date)

wwtp_data_v2 = wwtp_data[order(wwtp_data$WWTP, wwtp_data$date),]

wwtp_labs = wwtp_data %>%
  select(WWTP, lab_submitter) %>%
  distinct()

##########################################################################################################
#
# create wastewater data variables, some for use later, others to characterize WWTP traits
#
##########################################################################################################

wwtp_data_grouped = wwtp_data_v2[which(wwtp_data_v2$WWTP != "WWTPName"),] %>% 
  group_by(WWTP) %>%
  mutate(firstdif.ww = sars_cov2_adj_load - lag(sars_cov2_adj_load, 1),
         linkeddif.ww = sars_cov2_adj_load / lag(sars_cov2_adj_load, 1),
         
         # These create rolling and de-trended rolling variables, but aren't used in analysis
         # But this way we have them
         past3avg.ww = rollsumr(sars_cov2_adj_load, 3, fill=NA),
         past3avg.firstdif.ww = past3avg.ww - lag(past3avg.ww, 1),
         past3avg.linkeddif.ww = past3avg.ww / lag(past3avg.ww, 1),
         
         days_since_last_sample = date - lag(date,1),
         
         # Some WWTPs changes sampling frequency over time, assign value to frequency
         # to compare correlation within specific WWTPs 
         
         frequency_rollmean = rollmean(as.numeric(days_since_last_sample), 5, fill = NA, align = "right"),
         
         frequency_category = case_when(frequency_rollmean < 6 ~ "Multiple times per week", # Every day to every 5 days
                                        frequency_rollmean < 9 ~ "About once per week", # Every 6, 7, or 8 days
                                        frequency_rollmean > 8 ~ "Less than once per week") # Every 9 days or more
  ) %>%
  
  # Remove unnecessary variables
  
  select(-c(sewershed_trend:status, pdc)) %>%
  
  ungroup()

# Only consider cases after Jan 2021, when lab methods switched at state lab singleplex to multiplex (current)
postJAN_wwtp = wwtp_data_grouped[which(wwtp_data_grouped$date > as.Date("2021-01-31")),]
min(as.Date(postJAN_wwtp$date))

abstract_wwtp = postJAN_wwtp
max(as.Date(abstract_wwtp$date))

abstract_wwtp = abstract_wwtp[order(abstract_wwtp$WWTP, abstract_wwtp$date),] %>%
  group_by(WWTP) %>%
  mutate(mean_pop_served = mean(population_served, na.rm=T), #population served should be static in dataset, but just in case)
         mean_days_between_samples = mean(as.numeric(days_since_last_sample), na.rm=T),
         # mean_flowrate = mean(average_flow_rate), #avg_flow_rate itself is average of the 24h composite sample period
         # mean_flowrate_liters = mean(flow_rate_l), #this is avg_flow_rate in liters
  ) %>%
  ungroup() %>%
  mutate(population_quartile = ntile(mean_pop_served, 4),
         population_category = case_when(mean_pop_served < 10000 ~ "1: <10k",
                                         mean_pop_served < 20000 ~ "2: 10 to <20k",
                                         mean_pop_served < 50000 ~ "3: 20k to <50k",
                                         mean_pop_served >= 50000 ~ "4: 50k plus"))

# Did sampling frequency change during surveillance period
abstract_wwtp = abstract_wwtp %>%
  group_by(WWTP) %>%
  mutate(freq.diff = max(frequency_rollmean) - min(frequency_rollmean),
         changed_freq_yn = case_when(freq.diff < 2 ~ "No", freq.diff >= 2 ~ "Yes"))


##########################################################################################################
#
# create wastewater data quantiles, window of time related to actual days - not number of samples
# based on analyses, 90d window seems to work best
#
##########################################################################################################

workset = abstract_wwtp

# To calculate the rolling quantiles over past N days (rather than time points) --------------------------

# This should help account for WWTPs with different sampling frequencies
# Create list of all days between start and stop of workset for each WWTP
days = (rep(seq(min(workset$date),max(workset$date), by="days"), times = length(unique(workset$WWTP))))
wwtps = rep(unique(workset$wwtp_name), each=length(seq(min(workset$date),max(workset$date), by="days")))
wwtps_days = data.frame(days, wwtps)
names(wwtps_days) = c("date", "WWTP")
workset2 = left_join(wwtps_days, workset, by=c("date", "WWTP"))

# Two duplicate entries removed
df_dups = workset2[c("date", "WWTP")]
workset2 = workset2[!duplicated(df_dups),]

# Fill in missing ww data for each date with most recent ww sample value --------------------------------

# This does not impact the quantile value over time. 
workset2 = workset2 %>%
  mutate(sars_cov2_adj_load_carryforward = zoo::na.locf(sars_cov2_adj_load, na.rm = F),
         sars_cov2_adj_load_log10 = log(sars_cov2_adj_load_carryforward, base = 10))

##################################################################################################################################
#
# Create quantile cut-offs, this one at 80% (i.e., what is the 80th tile over the past N days)
#
##################################################################################################################################

workset3 = workset2[order(workset2$WWTP, workset2$date),] %>%
  group_by(WWTP) %>%
  mutate(
    ntile_80 = rollapply(sars_cov2_adj_load_log10, width = 30, FUN = "quantile", p = 0.80, align = "right", na.rm=T, fill=NA)) %>%
  select(WWTP, population_served, lab_submitter, date, 
         past3avg.ww:past3avg.linkeddif.ww,
         days_since_last_sample,
         sars_cov2_adj_load, sars_cov2_adj_load_log10, sars_cov2_adj_load_carryforward,
         ntile_80) %>%
  ungroup()


##################################################################################################################################
#
# Create ranked percentile for individual ww concentration levels, to show relative magnitude of ww compared to past N days
#
##################################################################################################################################

# MEthod 1 (not ideal, considers days instead of N observations, which means %tiles not exactly right) --------------------------------------------------------------

# # Define window of observations to use
# N = 30
# 
# # Remove WWTPS without sufficient number of observations
# filtered = workset3[order(workset3$WWTP, workset3$date),] 
# 
# # Create rolling rank variable based on past N observations
# # Interpreted as: "This ww value was higher than X% of values in the past N observations"
# workset4 = filtered[order(filtered$WWTP, filtered$date),] %>%
#   group_by(WWTP) %>%
#   mutate(roll_rank = TTR::runPercentRank(sars_cov2_adj_load_log10, n = N)) 
# 
# test60 = workset4 %>%
#   filter(WWTP == "Appleton WWTF")
# 
# # Check to see that this function is calculating ranked percentile okay
# test60 = test60[c(as.numeric(nrow(test60) - (N-1)):nrow(test60)),]
# test60$rank = (rank(test60$sars_cov2_adj_load_log10) - 0.5)/N

# Method 2 (consider only the N most recent dates from each WWTP and rejoin with main set) ------------------------------------------------------------

# Define window of observations to use
N = 60

# Consider only data collected in the past N days
filtered = workset3[order(workset3$WWTP, workset3$date),] %>%
  filter(date >= max(date) - (N - 1)) %>%
  # And consider only non-missing ww data
  filter(!is.na(sars_cov2_adj_load))

# Create for loop to capture the percentile for each ww treatment facility
unique_wwtps = unique(filtered$WWTP)
rolling_rank_estimates = as.data.frame(matrix(ncol=3, nrow=0))
colnames(reg_estimates) = c("WWTP", "date", "rolling_rank")

for (i in 1:length(unique_wwtps)){
  
  print(paste(unique_wwtps[i]))
  # Isolate each WWTP in turn
  test = filtered %>% filter(WWTP == unique_wwtps[[i]])
  # Calculate the percentile rank for each ww value in the past N days
  test$rolling_rank = (rank(test$sars_cov2_adj_load_log10))/nrow(test)
  # Extract data for the most recent date
  lastday = test %>% filter(date == max(date)) %>% select(WWTP, date, rolling_rank)
  # Bind to rolling rank estimates dataset
  rolling_rank_estimates = rbind(rolling_rank_estimates, lastday)
  
}

# Incorporate the most recent day's ntile with the full dataset
workset4 = left_join(workset3, rolling_rank_estimates, by = c("WWTP", "date"))

# Once quantiles have been created, remove all days without actual ww sample data -----------------------------
# This keeps the quantile levels based on actual days, but subsequent regression model based on # previous samples (not days)

workset4 = workset4 %>%
  filter(!is.na(sars_cov2_adj_load))

# Finally, create variables based on average of previous 3 ww samples ------------------------------------------

# If avg above the established quantile, (and slope passes threshold), then flag raised

K = 3

workset5 = workset4[order(workset4$WWTP, workset4$date),] %>%
  group_by(WWTP) %>%
  mutate(pastKavg.wwlog10 = rollmean(sars_cov2_adj_load_log10, K, align = "right", na.pad = T),
         pastKpercentchange.ww = (pastKavg.wwlog10 - lag(pastKavg.wwlog10,1))/lag(pastKavg.wwlog10)*100)

##########################################################################################################
#
# Calculate trends in ww data based on regression line of previous samples
# NOT day based, unlike the quantile calculation
#
##########################################################################################################

# Set up empty dataframe for regression estimate calculations in FOR LOOP
reg_estimates = as.data.frame(matrix(ncol=7, nrow=0))
colnames(reg_estimates) = c("WWTP", "date", "ww_days_elapsed", "lmreg_n" , "lmreg_slope", "lmreg_sig", "modeled_percentchange")

# Establish date parameter - STATIC
maxdate = max(abstract_wwtp$date)
maxdate

distinct_wwtps = workset4 %>%
  distinct(WWTP)

# Define window to calculate regression (define as number - 1)
# (i.e. for 5 observation window, S = 4)
S = 4

for (i in 1:nrow(distinct_wwtps)){
  
  print(paste(distinct_wwtps[i,1]))
  ww.x = workset5 %>%
    filter(WWTP==paste(distinct_wwtps[i,1]))
  
  for (k in 1:(nrow(ww.x) - S)){
    ww.x.subset = ww.x[c(k:(k+S)),] 
    lm.subset = lm(sars_cov2_adj_load_log10 ~ date, # date included works same as days_since_last_sample 
                   data = ww.x.subset) 
    summary(lm.subset)
    
    # Extract row to bind with workset
    ww.x.tobind = ww.x.subset %>%
      filter(date == max(date)) %>%
      select(WWTP, date) %>%
      mutate(ww_days_elapsed = as.numeric(max(ww.x.subset$date) - min(ww.x.subset$date)),
             lmreg_n = nrow(ww.x.subset), 
             lmreg_slope = summary(lm.subset)$coefficients[2,1],
             lmreg_sig = summary(lm.subset)$coefficients[2,4],
             modeled_percentchange = (10^(lmreg_slope*ww_days_elapsed)-1)*100
      )
    
    # Join with full set of reg estimates
    reg_estimates = rbind(reg_estimates, ww.x.tobind)
    
  }
}

# Bind regression estimates to workset 

workset6 = left_join(workset5, reg_estimates, by=c("WWTP", "date"))

# Create flag indicators (CDC trend, trend+quantile, trend+quantile+pvalue) --------------------------------------------------------

pval = 0.3

workset7 = workset6 %>%
  
  mutate(
    
    wastewatersampledate = 1,
    
    change_m2 = case_when(modeled_percentchange > 100 ~ "INCREASE - Major",
                          modeled_percentchange > 9 ~ "INCREASE - Moderate", 
                          modeled_percentchange > -10 ~ "Fluctuating", 
                          modeled_percentchange > -100  ~ "DECREASE - Moderate",
                          modeled_percentchange < -99  ~ "DECREASE - Major"),
    
    # CDC trend --------------------------------------------------------------------------
    cdc_flag = case_when(change_m2 == "INCREASE - Major" ~ "Flag"),
    
    # CDC trend plus quantile ------------------------------------------------------------
    flag_ntile80 = case_when(pastKavg.wwlog10 > ntile_80 & cdc_flag == "Flag" ~ "Flag"),
    
    # CDC trend plus quantile plus p value -----------------------------------------------
    flag_ntile80_pval = case_when(pastKavg.wwlog10 > ntile_80 & cdc_flag == "Flag" & lmreg_sig < pval ~ "Flag"))

# Above set contains flags on certain days of wastewater surveillance for each WWTP

# Interpretation, "We are seeing increasing wastewater concentrations (based on reg slope) 
# which are, on average (based on pastkavg.wwlog10), considerably higher than levels
# observed in the past 90 days (based on 80th ntile)"



