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
window_buff <- floor(length_window_days_pass1/2) #+/- number of days from center day to form moving window

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

    #get number of days
    dat.DoY_windows <- seq(from = min(dat$DoY, na.rm=T) + window_buff,
                           to = max(dat$DoY, na.rm=T) - window_buff)

    # now filter data for nighttime only data
    dat.night <- dat %>%
      filter(daynight == "night")

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
          mod05 <- minpack.lm::nlsLM(y05 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                  T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)

          #use broomExtra to capture model output parameters
          tidymod05 <- broomExtra::tidy(mod05)
          E0_est <- tidymod05$estimate[2]
          uncert <- tidymod05$std.error[2]

          isValid <- T #likely to remove this

          if(isValid){
            E0_v_U05[day] <- E0_est # save the daily estimated parameter into E0 vector
            unc_v_U05[day] <- uncert
          }


      }

      if(!all(is.na(y50))){
          # run model and assign to object
          mod50 <- minpack.lm::nlsLM(y50 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                   T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
          # use broomExtra to capture model output parameters
          tidymod50 <-broomExtra::tidy(mod50)
          E0_est <- tidymod50$estimate[2]
          uncert <- tidymod50$std.error[2]

          isValid <- T #likely to remove this

          if(isValid){
            E0_v_U50[day] <- E0_est # save the daily estimated parameter into E0 vector
            unc_v_U50[day] <- uncert
          }

      }

      if(!all(is.na(y95))){
        # run model and assign to object
        mod95 <- minpack.lm::nlsLM(y95 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                                   T0 = T0_K, E0), start = list(Rref = 1, E0 = 100), model = TRUE)
        # use broomExtra to capture model output parameters
        tidymod95 <-broomExtra::tidy(mod95)
        E0_est <- tidymod95$estimate[2]
        uncert <- tidymod95$std.error[2]

        isValid <- T #likely to remove this

        if(isValid){
          E0_v_U95[day] <- E0_est # save the daily estimated parameter into E0 vector
          unc_v_U95[day] <- uncert
        }
      }

    } #end rolling 7-day window for-loop

  E0_annual[[toString(year)]][1] <- mean(E0_v_U05, na.rm = T)
  E0_annual[[toString(year)]][2] <- mean(E0_v_U50, na.rm = T)
  E0_annual[[toString(year)]][3] <- mean(E0_v_U95, na.rm = T)

  } #end this site's files for-loop

  # Save annual E0 vals
  outfile <- paste0(data.E0.case, site, "_E0_annual_df.rds")
  saveRDS(E0_annual, outfile)

} #end sites for-loop
