rm(list=ls())

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))

# Load packages, global functions and variables included in config script
source("config.r")


#ggplot format function
legend_theme <- function() {
  theme(
    legend.title=element_blank(),
    legend.text=element_text(size=10,color="black"),
    legend.background=element_rect(fill = "white", color = "black"),
    legend.key=element_blank(),
    legend.key.size=unit(1,"lines"),
    #legend.spacing.y=unit(0.1, "cm"),
    legend.position=c(0.87, 0.2),
  )
}



#########################################
# Begin main
#########################################

# create empty lists of ggplot objects ("grobs") to add to within loops
grob.list.all <- vector(mode = "list", length = length(site.list))
grob.list.all.weekly <- vector(mode = "list", length = length(site.list))
grob.list.nee <- vector(mode = "list", length = length(site.list))
names(grob.list.all) <- names(grob.list.nee) <- names(grob.list.all.weekly) <- site.list

# get list of all files that exist
files.list <- lapply(site.list, FUN = function(x) list.files(data.case, pattern = x,
                     full.names = T))
names(files.list) <- site.list

# Loop through each site
for(site in site.list){
  
  site.files <- files.list[[site]] #get files corresponding to this site
  site.df <- data.frame() %>% as_tibble()
  
  for(file in site.files){
    #if its a summary file, skip it
    if(grepl("part_s", file))
      next
    flux.df <- readRDS(file)
    site.df %<>%
      bind_rows(flux.df)
  }
  
  year1 <- head(site.df$Year, 1)
  yearn <- tail(site.df$Year, 1)
  filename <- paste0(data.case, site, "_flux_part_s", year1, "-", yearn, "_df.rds")
  saveRDS(site.df, filename)
  
  
  summary.df <- site.df %>%
    mutate(Week = lubridate::week(posix_time),
           YearWeek = paste(lubridate::year(posix_time), Week, sep = "_")) %>% #get week of year
    relocate(Week, .after = Year) %>% #put Week in the right column order after Year, before DoY
    relocate(YearWeek, .after = Year) #put the unique YearWeek column after Year
    
  
  # get weekly avg posix times
  weekly.df <- summary.df %>% 
    group_by(YearWeek) %>%
    summarise_at(vars(posix_time), list(mean_posix = mean))
  
  # get weekly avg NEE
  weekly.df <- summary.df %>%
    #filter(Hour > 10 & Hour < 14) %>%
    #filter(daynight == "night") %>%
    group_by(YearWeek) %>%
    summarise_at(vars(NEE_U50_f), list(weekly_NEE = mean)) %>% 
    left_join(weekly.df, by = "YearWeek")
  
  # get weekly avg GPP
  weekly.df <- summary.df %>%
    group_by(YearWeek) %>%
    summarise_at(vars(GPP_U50_f), list(weekly_GPP = mean)) %>% 
    left_join(weekly.df, by = "YearWeek")
  
  # get weekly avg Reco
  weekly.df <- summary.df %>%
    group_by(YearWeek) %>%
    summarise_at(vars(Reco_U50), list(weekly_Reco = mean)) %>% 
    left_join(weekly.df, by = "YearWeek")


  ylims <- c(min(site.df$NEE_U50_f, na.rm = T), min(max(site.df$GPP_U50_f, na.rm = T), 50))
  
  #### plot data
  
  # NEE, GPP, Reco
  grob <- site.df %>% # take data
    ggplot() + theme_FS() + # pipe in ggplot call
    geom_point(aes(x=posix_time, y=NEE_U50_f, color = "gray50"), size=0.5) + # plot NEE vs time
    geom_point(aes(x=posix_time, y=GPP_U50_f, color = "red"), size = 0.5) + # add GPP vs time
    geom_point(aes(x=posix_time, y=Reco_U50, color = "blue"), size = 0.5) + # add Respiration vs time
    ylim(ylims) +
    xlab("") + scale_x_datetime(labels = label_date("%b %Y")) +
    ylab("[umol m-2 s-1]") +
    ggtitle(site) +
    legend_theme() +
    scale_color_identity(breaks = c("gray50", "red", "blue"), 
                         labels = c("NEE", "GPP", "Reco"), guide = "legend")
  
  grob.list.all[[site]] <- grob
  
  # NEE, weekly GPP, weekly Reco
  grob <- site.df %>% # take data
    ggplot() + theme_FS() + # pipe in ggplot call
    geom_point(data = weekly.df, aes(x=mean_posix, y=weekly_NEE, color = "gray50"),  size = 0.75) + # add NEE vs time
    geom_line(data = weekly.df, aes(x=mean_posix, y=weekly_NEE, color = "gray50"), size = 0.75) + # add NEE vs time
    geom_point(data = weekly.df, aes(x=mean_posix, y=weekly_GPP, color = "red"),  size = 0.75, alpha = 0.5) + # add GPP vs time
    geom_line(data = weekly.df, aes(x=mean_posix, y=weekly_GPP, color = "red"), size = 0.75, alpha = 0.5) + # add GPP vs time
    geom_point(data = weekly.df, aes(x=mean_posix, y=weekly_Reco, color = "blue"), size = 0.75, alpha = 0.5) + # add Reco vs time
    geom_line(data = weekly.df, aes(x=mean_posix, y=weekly_Reco, color = "blue"), size = 0.75, alpha = 0.5) + # add Reco vs time
    ylim(c(-6, 6)) +
    xlab("") + scale_x_datetime(labels = label_date("%b %Y")) +
    ylab("[umol m-2 s-1]") +
    ggtitle(site) +
    legend_theme() +
    scale_color_identity(breaks = c("gray50", "red", "blue"), 
                         labels = c("NEE", "GPP", "Reco"), guide = "legend")
  
  grob.list.all.weekly[[site]] <- grob
  
  # Just NEE
  grob <- site.df %>% # take data
    ggplot() + theme_FS() + # pipe in ggplot call
    geom_point(aes(x=posix_time, y=NEE_U50_f, color = "gray50"), size=0.5) + # plot NEE vs time
    ylim(ylims) +
    xlab("") + scale_x_datetime(labels = label_date("%b %Y")) +
    ylab("[umol m-2 s-1]") +
    ggtitle(site) +
    legend_theme() +
    scale_color_identity(breaks = c("gray50"), labels = c("NEE"), guide = "legend")
  
  grob.list.nee[[site]] <- grob
}

gridRows <- ceiling(length(site.list)/2)
gridCols <- 2
main.title <- "Half-hourly partitioned fluxes at AMF & NEON sites\n(derived using 50% u* threshold)"

#print PNG

# HH NEE, GPP, Reco
png_outfile <- paste0(plot.dir, "flux_all_sites_all_years.png") #define plot filename
png(png_outfile, width = 800, height = 1000) #open pdf
grid.arrange(grobs = grob.list.all, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)))
dev.off() #close pdf


main.title <- "Weekly partitioned fluxes at AMF & NEON sites\n(derived using 50% u* threshold)"

# HH NEE, weekly GPP, weekly Reco
png_outfile <- paste0(plot.dir, "flux_all_sites_all_years_weekly.png") #define plot filename
png(png_outfile, width = 800, height = 1000) #open pdf
grid.arrange(grobs = grob.list.all.weekly, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)))
dev.off() #close pdf

# Now plot just NEE
main.title <- "Half-hourly net ecosystem exchange (NEE) at AMF & NEON sites\n(derived using 50% u* threshold)"

png_outfile <- paste0(plot.dir, "NEE_all_sites_all_years.png") #define plot filename
png(png_outfile, width = 800, height = 1000) #open pdf
grid.arrange(grobs = grob.list.nee, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)))
dev.off() #close pdf
