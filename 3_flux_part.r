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

window_buff <- floor(length_window_days_pass2/2) #+/- number of days from center day to form moving window

#########################################
# Begin main
#########################################

# get list of all files that exist
files.list <- lapply(site.list, FUN = function(x) list.files(data.main, pattern = x))
names(files.list) <- site.list

# Loop through each site
for(site in site.list){
  
  site.files <- files.list[[site]] #get files corresponding to this site
  
  # Load E0 Estimates from previous step
  E0_annual_file <- paste0(data.E0.case, site, "_E0_annual_df.rds")
  E0_annual <- readRDS(E0_annual_file)
  
  # for this site, loop through all the available data files (1 per year)
  for(file in site.files){
   
    message(paste0("Loading file: ", file)) #print to console to give some play by play
    
    ### Prepare Data for partitioning
    # read data
    dat <- readRDS(paste0(data.main, file))
    year <- dat$Year[1]
    E0_vals <- E0_annual[[toString(year)]]
    
    # Define the Day of Year windows based on length of moving window
    dat.DoY_windows <- seq(from = min(dat$DoY, na.rm=T) + window_buff,
                           to = max(dat$DoY, na.rm=T) - window_buff)
    
    # now filter data for nighttime only data
    dat.night <- dat %>%
      filter(daynight == "night")
    
    # Set up empty vectors for E0, Rref to be filled in the next for-loop
    E0_v_U05 <- rep(NA, ndays)
    Rref_v_U05 <- rep(NA, ndays)
    E0_v_U50 <- rep(NA, ndays)
    Rref_v_U50 <- rep(NA, ndays)
    E0_v_U95 <- rep(NA, ndays)
    Rref_v_U95 <- rep(NA, ndays)
    
    ###########################################################
    ## Loop through moving window and solve for fit parameters
    ###########################################################
    message(paste0("Calculating daily estimates of R_ref & E_0 for file: ", file))
    
    # Loop through the moving 7-day window
    for(day in dat.DoY_windows){
      
      #Filter data for this 7-day window
      DoY_window <- seq(day - window_buff, day + window_buff)
      dat.window <- dat.night %>%
        filter(DoY %in% DoY_window) %>%
        dplyr::select(any_of(window.vars))
      
      if(nrow(dat.window) == 0)
        next
      
      x <- dat.window$Tair_K #grab your "x variable", i.e. input to Reco function
      y05 <- dat.window$NEE_U05_f #grab your "y variable", i.e. output to Reco function
      y50 <- dat.window$NEE_U50_f #grab your "y variable", i.e. output to Reco function
      y95 <- dat.window$NEE_U95_f #grab your "y variable", i.e. output to Reco function
      Tref_K <- unique(dat.window$AMNT_K) #can only be one value
      
      
      # Define a "percent data utility" to see what percent of the HH timestamps in this
      # window are considered "nighttime", and actually have usable data
      perc_util <- round(length(which(!is.na(y50)))*100/(length_window_days_pass2*24*2), 1)
      
      # If percent utility is <10%, let's consider the data to be too sparse to 
      # determine meaningful model parameters for this window. Skip this 7-day window.
      # Also skip if there is NEE data but no Temp data (does happen sometimes)
      if(perc_util < 10 | all(is.na(x)))
        next
      
      if(!all(is.na(y05)) & !is.na(E0_vals[1])){
        # run model and assign to object
        mod05 <- minpack.lm::nlsLM(y05 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                 T0 = T0_K, E0 = E0_vals[1]), start = list(Rref = 1), model = TRUE)
        # use broomExtra to capture model output parameters
        tidymod05 <- broomExtra::tidy(mod05)
        # save the daily estimated parameters into the Rref and E0 vectors
        Rref_v_U05[day] <- tidymod05$estimate[1]
        E0_v_U05[day] <- E0_vals[1]#tidymod05$estimate[2]
      }
      
      if(!all(is.na(y50)) & !is.na(E0_vals[2])){
        # run model and assign to object
        mod50 <- minpack.lm::nlsLM(y50 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                 T0 = T0_K, E0 = E0_vals[2]), start = list(Rref = 1), model = TRUE)
        # use broomExtra to capture model output parameters
        tidymod50 <-broomExtra::tidy(mod50)
        # save the daily estimated parameters into the Rref and E0 vectors
        Rref_v_U50[day] <- tidymod50$estimate[1]
        E0_v_U50[day] <- E0_vals[2]#tidymod50$estimate[2]
        # message(paste0("Performing parameter estimation, day ", day, ", ",
        #                perc_util, "% data used. Rref = ", round(Rref_v_U50[day], 1), 
        #                ", E0 = ", round(E0_v_U50[day], 2)))
        
        # if you wanted to plot as debug, do this:
        # combine initial vectors with predicted values and create plot in ggplot
        if(plot_debug){
          Sys.sleep(1)
          graphics.off()
          g1 <- as_tibble(cbind("Tair"=x,"NEE"=y50,".fitted"=predict(mod50))) %>% # bind variables
          ggplot() + # pipe ggplot
          geom_point(aes(x=Tair, y=NEE), size=0.5) +
          geom_line(aes(x=Tair, y=.fitted), color = "purple", size=1) + # plot predicted values as a geom_line
          ggtitle(paste("DoY", day)) + 
          theme_FS()
          print(g1)
          
          }
      } #end NA check on y50
      
      if(!all(is.na(y95)) & !is.na(E0_vals[3])){
        # run model and assign to object
        mod95 <- minpack.lm::nlsLM(y95 ~ fNight_Reichstein(Temp = x, Rref, Tref = Tref_K,
                                   T0 = T0_K, E0 = E0_vals[3]), start = list(Rref = 1), model = TRUE)
        # use broomExtra to capture model output parameters
        tidymod95 <-broomExtra::tidy(mod95)
        # save the daily estimated parameters into the Rref and E0 vectors
        Rref_v_U95[day] <- tidymod95$estimate[1]
        E0_v_U95[day] <- E0_vals[3]#tidymod95$estimate[2]
      }
      
    } #end rolling 7-day window for-loop

    #When Rref is less than zero, set it to zero
    ## IMPORTANT: SKIP FOR NOW!
    ## TODO: uncomment if you wish to set all negative GPP to zero
    Rref_v_U05 = replace(Rref_v_U05, Rref_v_U05 < 0, 0) 
    Rref_v_U50 = replace(Rref_v_U50, Rref_v_U50 < 0, 0)
    Rref_v_U95 = replace(Rref_v_U95, Rref_v_U95 < 0, 0)
    
    # Apply the daily Rref and E0 values to the main dataset; calculate Reco and GPP
    message(paste0("Calculating Reco and GPP for file: ", file))
    dat <- dat %>%
      mutate(Rref_U05 = Rref_v_U05[DoY], E0_U05 = E0_v_U05[DoY], 
             Rref_U50 = Rref_v_U50[DoY], E0_U50 = E0_v_U50[DoY],
             Rref_U95 = Rref_v_U95[DoY], E0_U95 = E0_v_U95[DoY])
    
    ### 5% u* threshold
    # check if the E0 value is NULL. if so, save list of NA vals to Reco and GPP
    if(is.na(E0_vals[1])){
      dat <- dat %>%
        mutate(Reco_U05 = as.numeric(rep(NA, nrow(dat))),
               GPP_U05_f = as.numeric(rep(NA, nrow(dat))))
    } else{
      dat <- dat %>%
        mutate(Reco_U05 = fNight_Reichstein(Tair_K, Rref_U05, Tref_K, T0_K, E0_U05),
               GPP_U05_f = -1 * (NEE_U05_f - Reco_U05))
        #iNight <- which(dat$daynight == "night")
        #dat$GPP_U05_f[iNight] <- 0
        #dat$Reco_U05[iNight] <- dat$NEE_U05_f[iNight]
      
      ## TODO: uncomment if you wish to set all negative GPP to zero
      # inegGPP <- which(dat$GPP_U05_f < 0)
      # dat$GPP_U05_f[inegGPP] <- 0
      # dat$Reco_U05[inegGPP] <- dat$NEE_U05_f[inegGPP]
    }
    
    ### 50% u* threshold
    # check if the E0 value is NULL. if so, save list of NA vals to Reco and GPP
    if(is.na(E0_vals[2])){
      dat <- dat %>%
        mutate(Reco_U50 = as.numeric(rep(NA, nrow(dat))),
               GPP_U50_f = as.numeric(rep(NA, nrow(dat))))
    } else{
      dat <- dat %>%
        mutate(Reco_U50 = fNight_Reichstein(Tair_K, Rref_U50, Tref_K, T0_K, E0_U50),
               GPP_U50_f = -1 * (NEE_U50_f - Reco_U50))
      #iNight <- which(dat$daynight == "night")
      #$GPP_U50_f[iNight] <- 0
      #dat$Reco_U50[iNight] <- dat$NEE_U50_f[iNight]
      
      ## TODO: uncomment if you wish to set all negative GPP to zero
      # inegGPP <- which(dat$GPP_U50_f < 0)
      # dat$GPP_U50_f[inegGPP] <- 0
      # dat$Reco_U50[inegGPP] <- dat$NEE_U50_f[inegGPP]
    }
    
    ### 95% u* threshold
    # check if the E0 value is NULL. if so, save list of NA vals to Reco and GPP
    if(is.na(E0_vals[3])){
      dat <- dat %>%
        mutate(Reco_U95 = as.numeric(rep(NA, nrow(dat))),
               GPP_U95_f = as.numeric(rep(NA, nrow(dat))))
    } else{
      dat <- dat %>%
        mutate(Reco_U95 = fNight_Reichstein(Tair_K, Rref_U95, Tref_K, T0_K, E0_U95),
               GPP_U95_f = -1 * (NEE_U95_f - Reco_U95))
      #iNight <- which(dat$daynight == "night")
      #dat$GPP_U95_f[iNight] <- 0
      #dat$Reco_U95[iNight] <- dat$NEE_U95_f[iNight]
      
      ## TODO: uncomment if you wish to set all negative GPP to zero
      # inegGPP <- which(dat$GPP_U95_f < 0)
      # dat$GPP_U95_f[inegGPP] <- 0
      # dat$Reco_U95[inegGPP] <- dat$NEE_U95_f[inegGPP]
    }

    
    # Save data object
    outfile <- paste0(data.case, site, "_flux_part_", year, "_df.rds")
    saveRDS(dat, outfile)
    
    # # TODO: find a fix for this. Currently unable to plot if all GPP and Reco vals are NA
    # if(is.na(E0_vals[2]))
    #   next
    
    # plot data
    # grob <- dat %>% # take data
    #   ggplot() + theme_FS() + # pipe in ggplot call
    #   geom_point(aes(x=posix_time, y=NEE_U50_f), size=0.5, color = "gray50") + # plot NEE vs time
    #   geom_point(aes(x=posix_time, y=GPP_U50_f), size = 0.5, color = "red") + # add GPP vs time
    #   geom_point(aes(x=posix_time, y=Reco_U50), size = 0.5, color = "blue") + # add Respiration vs time
    #   xlab("Time of year") + scale_x_datetime(labels = date_format("%b")) +
    #   ylab("umol m-2 s-1") +
    #   ggtitle(paste("NEE (gray) & calculated GPP (red), Reco (blue) -", site, "-", year))
    # 
    # pdf_outfile <- paste0(plot.dir, site, "_flux_", year, ".pdf") #define plot filename
    # pdf(pdf_outfile, width = 8, height = 8) #open pdf
    # print(grob) #print plot to pdf
    # dev.off() #close pdf
    
  } #end this site's files for-loop
} #end sites for-loop
