rm(list=ls())

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))

# Load packages, global functions and variables included in config script
source("config.r")

#########################################
# Define script-specific vars
#########################################

# list of variables to keep in the moving window dataframe (keep things trim within the loops)
window.vars <- c("posix_time", "DoY", "Tair_K", "PPFD", "NEE_U05_f", "NEE_U50_f",
                 "NEE_U95_f", "AMT_K", "AMNT_K")

# robustness parameters
window_buff <- floor(length_window_days_pass1/7/2) #+/- number of days from center day to form moving window

#########################################
# Begin main
#########################################

# get list of all files that exist
files.list <- lapply(site.list, FUN = function(x) list.files(data.main, pattern = x))
names(files.list) <- site.list

# Loop through each site
for(site in site.list){
  
  site.files <- files.list[[site]] #get files corresponding to this site
  
  # Create object to hold the E0 values as we loop through this site's years of data
  E0_annual <- data.frame(matrix(NA, nrow = 3, ncol = length(site.files))) #3 rows for 5, 50 and 95% u* thresholds
  rownames(E0_annual) <- c("U05", "U50", "U95")
  
  # for this site, loop through all the available data files (1 per year)
  for(file in site.files){
   
    message(paste0("Loading file: ", file)) #print to console to give some play by play
    
    ### Prepare Data for partitioning
    # read data
    dat <- readRDS(paste0(data.main, file))
    year <- dat$Year[1]
    iyear <- which(site.files %in% file)
    colnames(E0_annual)[iyear] <- toString(year)
    
    # Define the Day of Year windows based on length of moving window
    dat.WoY_windows <- seq(from = min(dat$DoY, na.rm=T) + window_buff,
                           to = max(dat$DoY, na.rm=T) - window_buff)
    
    # now filter data for nighttime only data
    dat.night <- dat %>%
      filter(daynight == "night")
    
    
    
    ### TO-DO:
    # Do weekly averaging here
    # think about gaps and how the weekly averaging would affect things
    
    dat.night %<>%
      mutate(Week = lubridate::week(posix_time)) %>%
      relocate(Week, .after = Year)
    
    weekly.df <- dat.night %>% 
      group_by(Week) %>%
      summarize_all(mean, na.rm = T)
    
    
    # Set up empty vectors for E0, Rref to be filled in the next for-loop
    E0_v_U05 <- rep(NA, ndays)
    unc_v_U05 <- rep(NA, ndays)
    E0_v_U50 <- rep(NA, ndays)
    unc_v_U50 <- rep(NA, ndays)
    E0_v_U95 <- rep(NA, ndays)
    unc_v_U95 <- rep(NA, ndays)
    
    ###########################################################
    ## Loop through moving window and solve for fit parameters
    ###########################################################
    message(paste0("Calculating daily estimates of E_0 for file: ", file))
    
    # Loop through the moving 7-day window
    for(day in dat.DoY_windows){
      
      #Filter data for this 7-day window
      DoY_window <- seq(day - window_buff, day + window_buff)
      dat.window <- dat.night %>%
        filter(DoY %in% DoY_window) %>%
        select(any_of(window.vars))
      
      if(nrow(dat.window) == 0)
        next
      
      x <- dat.window$Tair_K #grab your "x variable", i.e. input to Reco function
      y05 <- dat.window$NEE_U05_f #grab your "y variable", i.e. output to Reco function
      y50 <- dat.window$NEE_U50_f #grab your "y variable", i.e. output to Reco function
      y95 <- dat.window$NEE_U95_f #grab your "y variable", i.e. output to Reco function
      Tref_K <- unique(dat.window$AMNT_K) #can only be one value
      
      
      # Define a "percent data utility" to see what percent of the HH timestamps in this
      # window are considered "nighttime", and actually have usable data
      perc_util <- round(length(which(!is.na(y50)))*100/(length_window_days_pass1*24*2), 1)
      
      # If percent utility is <10%, let's consider the data to be too sparse to 
      # determine meaningful model parameters for this window. Skip this 7-day window.
      # Also skip if there is NEE data but no Temp data (does happen sometimes)
      if(perc_util < 10 | all(is.na(x)))
        next
      
      if(!all(is.na(y05))){
        # run first round of model and assign to object
        mod05_primer <- minpack.lm::nlsLM(y05 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
        mod05_resid <- resid(mod05_primer)

        # Trim the first round of model to exclude upper and lower 5% of residuals
        itrim <- which(mod05_resid < quantile(mod05_resid, resid_low) |
                      mod05_resid > quantile(mod05_resid, resid_high))
        y05_trim <- y05
        y05_trim[itrim] <- NA

        ## Rerun model on trimmed dataset
        mod05 <- minpack.lm::nlsLM(y05_trim ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                 T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
        
        
        #use broomExtra to capture model output parameters
        tidymod05 <- broomExtra::tidy(mod05)
        E0_est <- tidymod05$estimate[2]
        uncert <- tidymod05$std.error[2]
        
        # Now get stats on temperature range, # of records and temp sensitivity
        # to determine if estimate is valid
        temp_range <- diff(range(x[-itrim], na.rm = T))
        num_records <- length(which(!is.na(y05_trim)))
        isValid <- (temp_range > min_T_range) & (num_records > min_records) &
                     (E0_est > E0_low) & (E0_est < E0_high)
        
        
        if(isValid){
          E0_v_U05[day] <- E0_est # save the daily estimated parameter into E0 vector
          unc_v_U05[day] <- uncert
        }
        
        
      }
      
      if(!all(is.na(y50))){
        # run first round of model and assign to object
        mod50_primer <- minpack.lm::nlsLM(y50 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                        T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
        mod50_resid <- resid(mod50_primer)
        
        # Trim the first round of model to exclude upper and lower 5% of residuals
        itrim <- which(mod50_resid < quantile(mod50_resid, resid_low) |
                         mod50_resid > quantile(mod50_resid, resid_high))
        y50_trim <- y50
        y50_trim[itrim] <- NA
        
        ## Rerun model on trimmed dataset
        mod50 <- minpack.lm::nlsLM(y50_trim ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                          T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
        
        # use broomExtra to capture model output parameters
        tidymod50 <-broomExtra::tidy(mod50)
        E0_est <- tidymod50$estimate[2]
        uncert <- tidymod50$std.error[2]
        
        # Now get stats on temperature range, # of records and temp sensitivity
        # to determine if estimate is valid
        temp_range <- diff(range(x[-itrim], na.rm = T))
        num_records <- length(which(!is.na(y50_trim)))
        isValid <- (temp_range > min_T_range) & (num_records > min_records) &
          (E0_est > E0_low) & (E0_est < E0_high)
        
        
        if(isValid){
          E0_v_U50[day] <- E0_est # save the daily estimated parameter into E0 vector
          unc_v_U50[day] <- uncert
        }

      }
      
      if(!all(is.na(y95))){
        # run first round of model and assign to object
        mod95_primer <- minpack.lm::nlsLM(y95 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                                                                  T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
        mod95_resid <- resid(mod95_primer)
        
        # Trim the first round of model to exclude upper and lower 5% of residuals
        itrim <- which(mod95_resid < quantile(mod95_resid, resid_low) |
                         mod95_resid > quantile(mod95_resid, resid_high))
        y95_trim <- y95
        y95_trim[itrim] <- NA
        
        ## Rerun model on trimmed dataset
        mod95 <- minpack.lm::nlsLM(y95_trim ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                                                                T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
        
        # use broomExtra to capture model output parameters
        tidymod95 <-broomExtra::tidy(mod95)
        E0_est <- tidymod95$estimate[2]
        uncert <- tidymod95$std.error[2]
        
        # Now get stats on temperature range, # of records and temp sensitivity
        # to determine if estimate is valid
        temp_range <- diff(range(x[-itrim], na.rm = T))
        num_records <- length(which(!is.na(y95_trim)))
        isValid <- (temp_range > min_T_range) & (num_records > min_records) &
          (E0_est > E0_low) & (E0_est < E0_high)

        if(isValid){
          E0_v_U95[day] <- E0_est # save the daily estimated parameter into E0 vector
          unc_v_U95[day] <- uncert
        }
      }
      
    } #end rolling 7-day window for-loop
    
    #get the 3 least uncertain values of E0 for each u* threshold
    E0_v_U05_top3 <- head(E0_v_U05[order(unc_v_U05)], 3)
    E0_v_U50_top3 <- head(E0_v_U50[order(unc_v_U50)], 3)
    E0_v_U95_top3 <- head(E0_v_U95[order(unc_v_U95)], 3)
    
    E0_annual[[toString(year)]][1] <- mean(E0_v_U05_top3, na.rm = T)
    E0_annual[[toString(year)]][2] <- mean(E0_v_U50_top3, na.rm = T)
    E0_annual[[toString(year)]][3] <- mean(E0_v_U95_top3, na.rm = T)
    
  } #end this site's files for-loop
  
  # Save annual E0 vals
  outfile <- paste0(data.E0.case, site, "_E0_annual_df.rds")
  saveRDS(E0_annual, outfile)
  
} #end sites for-loop