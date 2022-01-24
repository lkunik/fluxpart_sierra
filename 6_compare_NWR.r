rm(list=ls())

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))

# Load packages, global functions and variables included in config script
source("config.r")

#########################################
# Define script-specific vars
#########################################

compare.site <- "AMF_US-NR1"

# define NR1 REddyProc output data directory and file to use as input for partitioning
infiles <- sapply(NWR_years, FUN = function(x) list.files(path = data.dir.compare, pattern = paste0("AMF_US-NR1-", x)))


#########################################
# Begin main
#########################################


# get list of all files that exist
site.files.lk <- list.files(data.case, pattern = compare.site ,full.names = T)

dat.lk <- data.frame() %>% as_tibble()
  
for(file in site.files.lk){
  
  flux.df <- readRDS(file)
  dat.lk %<>%
    bind_rows(flux.df)
  
}

dat.compare <- data.frame() %>% as_tibble()

for(file in infiles){
  
  dat <- read.table(paste0(data.dir.compare, file), header = T, stringsAsFactors = F)[-1,]
  dat <- dat %>% as_tibble() #convert to tibble so we can work with this in dplyr language
  dat %<>% mutate_at(everything(.), as.numeric) #goddamn tibble has all cols as "char" type. Change to numeric
  
  year <- dat$Year[1]
  # define start of decimal days for conversion to POSIX
  posix_epoch <- ISOdatetime(year, 1, 1, 0, 0, 0, tz = "UTC")
  
  #get number of days
  ndays <- length(unique(dat$DoY))
  dat.DoY_windows <- seq(from = min(dat$DoY, na.rm=T) + 3,
                         to = max(dat$DoY, na.rm=T) - 3)
  
  
  # add column to indicate nighttime vs. daytime data
  dat <- dat %>%
    mutate(dDoY = DoY + (Hour/24),
           posix_time = julian.to.POSIX(dDoY, posix_epoch))
  
  dat.compare %<>%
    bind_rows(dat)
  
}



# plot GPP compare
grob <- dat.compare %>% # take data
  ggplot() + theme_FS() + # pipe in ggplot call
  geom_point(aes(x=posix_time, y=GPP_U50_f), size=0.5, color = "red") + # plot NEE vs time
  geom_point(data = dat.lk, mapping = aes(x=posix_time, y=GPP_U50_f), size = 0.5, color = "blue") + # add Respiration vs time
  xlab("Time of year") + scale_x_datetime(labels = date_format("%b %Y")) +
  ylab("umol m-2 s-1") +
  ggtitle(paste("Manual (blue) & REddyProc (red) GPP - US-NR1 (50% u* thresh)"))

pdf_outfile <- paste0(plot.dir, "US-NR1_GPP_compare.pdf") #define plot filename
pdf(pdf_outfile, width = 12, height = 8) #open pdf
print(grob) #print plot to pdf
dev.off() #close pdf

# plot Reco compare
grob <- dat.compare %>% # take data
  ggplot() + theme_FS() + # pipe in ggplot call
  geom_point(aes(x=posix_time, y=Reco_U50), size=0.5, color = "red") + # plot NEE vs time
  geom_point(data = dat.lk, mapping = aes(x=posix_time, y=Reco_U50), size = 0.5, color = "blue") + # add Respiration vs time
  xlab("Time of year") + scale_x_datetime(labels = date_format("%b %Y")) +
  ylab("umol m-2 s-1") +
  ggtitle(paste("Manual (blue) & REddyProc (red) Reco - US-NR1 (50% u* thresh)"))

pdf_outfile <- paste0(plot.dir, "US-NR1_Reco_compare.pdf") #define plot filename
pdf(pdf_outfile, width = 12, height = 8) #open pdf
print(grob) #print plot to pdf
dev.off() #close pdf
