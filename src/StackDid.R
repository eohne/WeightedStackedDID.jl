install.packages("data.table")
install.packages("tidyverse")
install.packages("ggthemes")
install.packages("rio")
install.packages("geomtextpath")
install.packages("gghighlight")
install.packages("collapse")
install.packages("fixest")
install.packages("modelsummary")


library(data.table)
library(tidyverse)
library(ggthemes)
library(rio)
library(geomtextpath)
library(gghighlight)
library(collapse)
library(fixest)
library(modelsummary)

dtc = fread("data/acs1860_unins_2008_2021.csv")


create_sub_exp = function(dataset, timeID, groupID, adoptionTime, focalAdoptionTime, kappa_pre, kappa_post){
  
  # Copy dataset 
  dt_temp = copy(dataset)

  # Determine earliest and latest time in the data. 
        # Used for feasibility check later
  minTime = dt_temp[, fmin(get(timeID))]
  maxTime = dt_temp[, fmax(get(timeID))]
  
  # Include only the treated groups and the clean controls that adopt at least kappa_post periods after the focal atime.
  dt_temp = dt_temp[get(adoptionTime) == focalAdoptionTime | get(adoptionTime) > focalAdoptionTime + kappa_post | get(adoptionTime) == TRUE | is.na(get(adoptionTime))]
  
  # Limit to time periods inside the event window defined by the kappas
  dt_temp = dt_temp[get(timeID) %in% (focalAdoptionTime - kappa_pre):(focalAdoptionTime + kappa_post)]
  
  # Make treatment group dummy
  dt_temp[, treat := 0]
  dt_temp[get(adoptionTime) == focalAdoptionTime, treat := 1] 
  
  # Make a post variable
  dt_temp[, post := fifelse(get(timeID) >= focalAdoptionTime, 1, 0)]
  
  # Make event time variable
  dt_temp[, event_time := get(timeID) - focalAdoptionTime]
  
  # Create a feasible variable
  dt_temp[, feasible := fifelse(focalAdoptionTime - kappa_pre >= minTime & focalAdoptionTime + kappa_post <= maxTime, 1, 0)]
  
  # Make a sub experiment ID.
  dt_temp[, sub_exp := focalAdoptionTime]
  
  return(dt_temp)
}; 

subexp2014 = create_sub_exp(
              dataset = dtc,
              timeID = "year",
              groupID = "statefips", 
              adoptionTime = "adopt_year", 
              focalAdoptionTime = 2014,
              kappa_pre = 3,
              kappa_post = 2)
nrow(subexp2014)
# Summarize
datasummary(All(subexp2014) ~ N + Mean + SD + Min + Max,
            data = subexp2014)


# create the sub-experimental data sets
events = dtc[is.na(adopt_year) == FALSE, funique(adopt_year)]

# make a list to store the sub experiments in.
sub_experiments = list()

# Loop over the events and make a data set for each one
for (j in events) {
  sub_name = paste0("sub_",j) 
  sub_experiments[[sub_name]] = create_sub_exp(
              dataset = dtc,
              timeID = "year",
              groupID = "statefips", 
              adoptionTime = "adopt_year", 
              focalAdoptionTime = j,
              kappa_pre = 3,
              kappa_post = 2)
}

# Vertically concatenate the sub-experiments
stackfull = rbindlist(sub_experiments)

# Remove the sub-experiments that are not feasible
stacked_dtc = stackfull[feasible == 1]

# Summarize
datasummary(All(stacked_dtc) ~ N + Mean + SD + Min + Max,
            data = stacked_dtc)




compute_weights = function(dataset, treatedVar, eventTimeVar, subexpVar) {

  # Create a copy of the underlying dataset
  stack_dt_temp = copy(dataset)

  # Step 1: Compute stack - time counts for treated and control
  stack_dt_temp[, `:=` (stack_n = .N,
                     stack_treat_n = sum(get(treatedVar)),
                     stack_control_n = sum(1 - get(treatedVar))), 
             by = get(eventTimeVar)
             ]  
  # Step 2: Compute sub_exp-level counts
  stack_dt_temp[, `:=` (sub_n = .N,
                     sub_treat_n = sum(get(treatedVar)),
                     sub_control_n = sum(1 - get(treatedVar))
                     ), 
             by = list(get(subexpVar), get(eventTimeVar))
             ]
  
  # Step 3: Compute sub-experiment share of totals
  stack_dt_temp[, sub_share := sub_n / stack_n]
  
  stack_dt_temp[, `:=` (sub_treat_share = sub_treat_n / stack_treat_n,
                     sub_control_share = sub_control_n / stack_control_n
                     )
             ]
  
  # Step 4: Compute weights for treated and control groups
  stack_dt_temp[get(treatedVar) == 1, stack_weight := 1]
  stack_dt_temp[get(treatedVar) == 0, stack_weight := sub_treat_share/sub_control_share]
  
  return(stack_dt_temp)
}  



stacked_dtc2 = compute_weights(
      dataset = stacked_dtc,
      treatedVar = "treat",
      eventTimeVar = "event_time",
      subexpVar = "sub_exp")

# Summarize
stacked_dtc2[event_time==0 & treat==0, 
             .(avg_control_weight = mean(stack_weight)), 
             by = sub_exp][order(sub_exp)]
mean(stacked_dtc2$stack_weight)
mean(stacked_dtc2$stack_weight[stacked_dtc2$sub_exp==2016])


weight_stack = feols(unins ~ i(event_time, treat, ref = -1) | treat + event_time, 
                              data = stacked_dtc2)

summary(weight_stack)
