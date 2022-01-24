rm(list=ls())

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))


# Load packages, global functions and variables included in config script
source("config.r")

#########################################
# Begin main
#########################################

# get list of all files that exist
files.list <- lapply(site.list, FUN = function(x) list.files(data.dir, pattern = x))
names(files.list) <- site.list

# Loop through each site
for(site in site.list){
  
  site.files <- files.list[[site]] #get files corresponding to this site
  
  # If site has PAR missing (i.e. AMF CZO sites), add it in from Mike Goulden's data
  PAR.file <- PAR.files[[site]] #NEON sites and AMF_US-NR1 don't need this so PAR.file will be NA
  PAR.df <- NA
  
  if(!is.na(PAR.file)){
    PAR.nc <- nc_open(PAR.file)
    # retrieve time variable, convert to POSIX, round to nearest 30 mins 
    df.time <- round(julian.to.POSIX(ncvar_get(PAR.nc, "TIME"), 
                    epoch = ISOdatetime(2006,1,1,0,0,0,tz = "UTC")), "30mins") 
    df.PAR <- ncvar_get(PAR.nc, "PAR_IN") #get incoming PAR (which is probably just SW rad in, but oh well)
    nc_close(PAR.nc) #close netcdf connection
    
    # convert to data.frame, then to tibble, because now we like tibbles :)
    PAR.df <- data.frame(posix_time = df.time, PAR = df.PAR) %>% as_tibble()
    # convert NaN values to NA. NaN is a matlab thing, but R plays nicer with NA
    PAR.df %<>% mutate_at(everything(.), ~replace(., is.nan(.), NA))
  }
  
  # for this site, loop through all the available data files (1 per year)
  for(file in site.files){
   
    message(paste0("Loading file: ", file)) #print to console to give some play by play
    
    ### Prepare Data for partitioning
    # read data
    dat <- read.table(paste0(data.dir, file), header = T, stringsAsFactors = F)[-1,]
    dat <- dat %>% as_tibble() #convert to tibble so we can work with this in dplyr language
    dat %<>% mutate_at(everything(.), as.numeric)  #goddamn tibble has all cols as "char" type. Change to numeric
    dat[dat == -9999] <- NA #now turn all -9999
    year <- dat$Year[1] #grab the year from the first index of the Year field
    
    # Filter out bad gap-filling flags. Also filters out NA NEE data, which is ok
    dat %<>% filter(NEE_U05_fqc < gf_flag_thresh,
                    NEE_U50_fqc < gf_flag_thresh,
                    NEE_U95_fqc < gf_flag_thresh)
    
    # extract year from first row of Year column and define start of decimal days for conversion to POSIX
    posix_epoch <- ISOdatetime(dat$Year[1], 1, 1, 0, 0, 0, tz = "US/Pacific")
    
    #get number of times
    ntime <- nrow(dat)
    
    # add columns for new time fields
    dat <- dat %>% mutate(Tair_K = Tair_f + 273.15, 
                    dDoY = DoY + (Hour/24),
                    posix_time = julian.to.POSIX(dDoY, posix_epoch))
    
    # change timezone of posix time to UTC so the join operation will assign
    # PAR to the correct timestamps
    attr(dat$posix_time, "tzone") <- "UTC" 
    
    dat <- dat %>% relocate(posix_time, .after = Hour)
    
    # If PAR.df is not NA, add the Goulden PAR data
    if(!all(is.na(PAR.df))){
      # Join in the column 
      dat <- left_join(dat, PAR.df, by = "posix_time")
      
      # Replace PPFD with new PAR
      dat %<>% mutate(PPFD = PAR)
    }
    
    #define "daynight" i.e day vs night based on PPFD
    dat %<>% mutate(daynight = ifelse(PPFD < night_PPFD_thresh, "night", "day"))
    
    # now filter data for nighttime only data
    dat.night <- dat %>%
      filter(daynight == "night")
    
    if(plot_debug){
      #dat.night %>%
      dat %>%
        #mutate(Year = year(posix_time), DoY = yday(posix_time), Hour = 
        #         hour(posix_time)) %>%
        ggplot() + 
        #geom_point(aes(x = posix_time, y = PAR))
        geom_point(aes(x = Hour, y = NEE_U50_f)) 
        #geom_point(. %>% filter(PPFD < 5), mapping = aes(x = posix_time,
        #          y = NEE_U50_f), col = "red")
      
    }
    
    # Get annual mean temp (AMT) and annual mean night temp (AMNT) to plug into the reference temperature
    AMT_K <- mean(dat$Tair_K, na.rm = T)
    AMNT_K <- mean(dat.night$Tair_K, na.rm = T)
    
    # Apply the daily Rref and E0 values to the main dataset; calculate Reco and GPP
    dat <- dat %>%
      mutate(AMT_K = rep(AMT_K, ntime), 
             AMNT_K = rep(AMNT_K, ntime))
    
    print(paste0("Annual Mean Temperature = ", round(AMT_K - 273.15, 1), "°C (", round(AMT_K, 1), " K)"))
    print(paste0("Annual Mean Nighttime Temperature = ", round(AMNT_K - 273.15, 1), "°C (", round(AMNT_K, 1), " K)"))
    
    # Save data object
    outfile <- paste0(data.main, site, "_flux_part_", year, "_df.rds")
    #saveRDS(dat, outfile)

  } #end this site's files for-loop
} #end sites for-loop