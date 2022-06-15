rm(list=ls())

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))

# Load packages, global functions and variables included in config script
source("config.r")

LR_cols <- c("Month", "PAR", "PPFD", "NEE_U50_f", "Month", "Hour", "posix_time")

#########################################
# Begin main
#########################################

# create empty lists of ggplot objects ("grobs") to add to within loops
grob.list <- vector(mode = "list", length = 12)
grob.list.midday <- vector(mode = "list", length = length(site.list))
names(grob.list.midday) <- site.list


# get list of all files that exist
files.list <- lapply(site.list, FUN = function(x) list.files(data.case, pattern = x,
                     full.names = T))
names(files.list) <- site.list


# Loop through each site
for(site in site.list){

  message(paste0(site, ":"))
  site.files <- files.list[[site]] #get files corresponding to this site
  site.df <- data.frame() %>% as_tibble()

  for(file in site.files){
    
    # for now, we are only going to make plots for CZ2, 2014
    if(file != "./data//AMF_US-CZ2_flux_part_2011_df.rds" &
       file != "./data//AMF_US-CZ2_flux_part_2013_df.rds" &
       file != "./data//AMF_US-CZ2_flux_part_2015_df.rds" &
       file != "./data//AMF_US-NR1_flux_part_2007_df.rds")
      next
    
    #if its a summary file, skip it
    if(grepl("part_s", file))
      next
    
    flux.df <- readRDS(file) %>% #load file 
      mutate(Week = lubridate::week(posix_time),
             Month = lubridate::month(posix_time)) #grab the week and month of the year so we can summarize
    
    year <- flux.df$Year[1] #get year

    #####################################################
    ### overall weekly midday light response curve
    #####################################################
    # 
    # #Calculate weekly mean NEE -midday hours only
    # flux.sum <- flux.df %>%
    #   filter(Hour > 10 & Hour < 14) %>% #filter for midday
    #   group_by(Week) %>%
    #   summarise_at(vars(NEE), list(NEE_midday = mean))
    # 
    # #Calculate weekly mean light -midday hours only
    # flux.sum <- flux.df %>%
    #   filter(Hour > 10 & Hour < 14) %>% #filter for midday
    #   group_by(Week) %>%
    #   summarise_at(vars(PPFD), list(PPFD_midday = mean)) %>%
    #   left_join(flux.sum, by = "Week")
    # 
    # main.title <- paste0("Light response curve @ ", site, " - ", year, 
    #                      "\nMidday only, weekly averaged, 50% u* threshold")
    # png_outfile <- paste0(plot.dir, "light_response_weekly_midday_", site, "-", 
    #                       year, ".png") #define plot filename
    # ylims <- c(-20, 10)
    # 
    # #plot the weekly NEE light response at midday
    # grob <- flux.sum %>%
    #   ggplot() + theme_FS() + # pipe in ggplot call
    #   geom_point(aes(x=PPFD_midday, y=NEE_midday), size=0.5) + # plot E0 vs time
    #   ylim(ylims) + #scale_y_log10() +
    #   xlab("PPFD [umol m-2 s-1]") +
    #   ylab("NEE [umol m-2 s-1]") +
    #   ggtitle(main.title)
    # png(png_outfile, width = 600, height = 500) #open pdf
    # print(grob)
    # dev.off() #close pdf
    # 
    #site.df %<>%
    #  bind_rows(flux.df)
    
    #####################################################
    ### monthly light response curves
    #####################################################
    
    
    flux.df %<>%
      dplyr::select(any_of(LR_cols)) #select columns we need for this function
    
    ylims <- c(-10, 5) #set y axis limits
    
    # Loop through months of year
    for(month in 1:12){
      
      month.df <- flux.df %>%
        filter(Month == month) #filter for this given month
      
      fit.df <- month.df %>%
        filter(!is.na(PPFD), !is.na(NEE_U50_f), PPFD > 0) #filter out bad data for the log fit
      
      # Estimate logarithmic relationship
      logEst <- lm(NEE_U50_f~log(PPFD), data = fit.df)
      xvec <- seq(0, max(month.df$PPFD, na.rm=T), by = 1)
      logPred <- predict(logEst, newdata = data.frame(PPFD=xvec))
      log.df <- data.frame(xvec, logPred)
      lab.txt <- paste0("NEE = ", round(logEst$coefficients[2], 1), "*log(PPFD) ", 
                        ifelse(logEst$coefficients[1] < 0,"- ", "+ "), 
                        round(abs(logEst$coefficients[1]), 1))
      
      grob <- month.df %>%
        ggplot() + theme_LK() + # pipe in ggplot call
        geom_point(aes(x=PPFD, y=NEE_U50_f), size=0.5, color = "gray40") + # plot E0 vs time
        geom_line(data = log.df, aes(x = xvec, y=logPred, color = "black"), size = 0.75) +
        ylim(ylims) + #scale_y_log10() +
        xlab("PPFD [umol m-2 s-1]") +
        ylab("NEE [umol m-2 s-1]") +
        scale_color_identity(breaks = "black", labels = lab.txt, guide = "legend") +
        ggtitle(month.abb[month])
      grob.list[[month]] <- grob #assign to the corresponding list of ggplot objects ("grobs")
      
    }
    
    gridRows <- 4
    gridCols <- 3
    main.title <- paste0("Light response & log fit @ ", site, " - ", year, ", 50% u* threshold")
    
    # plot monthly light response curve
    png_outfile <- paste0(plot.dir, "light_response_bymonth_", site, "-", year, ".png") #define plot filename
    png(png_outfile, width = 800, height = 1000) #open pdf
    grid.arrange(grobs = grob.list, nrow = gridRows, ncols = gridCols,
                 top = textGrob(main.title, gp = gpar(fontsize = 20)))
    dev.off() #close pdf
    
  } #end files for-loop
} #end sites for-loop
