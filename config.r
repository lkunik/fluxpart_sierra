#######################
### Configuration file (config.r)
### Eddy covariance flux partitioning - Sierra sites
#######################


#########################################
# Load packages
#########################################

library(tidyverse) # essential packages
library(broomExtra) # better model fitting functionality
library(minpack.lm) # better fits for non-linear models
library(FluxSynthU) # for plotting theme function theme_FS()
library(scales) 
library(ncdf4) #file reading/writing
library(flux) 
library(grid)
library(gridExtra) #for plotting multiple plots in same figure

#########################################
# Define Global Functions
#########################################

# Define POSIX conversion function
julian.to.POSIX <- function(julians,epoch=ISOdatetime(0001,1,1,0,0,0,tz="UTC")) {
  # Convention is that the Julian date of the epoch is 1.0 (not 0!)
  retval <- epoch+86400*(julians-1.0)
  attributes(retval)$tzone <- attributes(epoch)$tzone
  return(retval)
}

# Define Reco nighttime function: Temperature v NEE (from Reichstein et al 2005)
fNight_Reichstein <- function(Temp, Rref, Tref, T0, E0) {
  Rref * exp(E0*(1/(Tref - T0) - 1/(Temp - T0)))
}

#ggplot format function
theme_LK <- function() {
  theme(
    plot.title = element_text(color="black", hjust = 0.5, size=14, face="bold"),
    plot.subtitle = element_text(color="black", hjust = 0.5, size=12),
    axis.line=element_line(size=0.75),
    axis.title.x=element_text(size=15,color="black"),
    axis.text.x=element_text(size=12,color="black"),
    axis.title.y=element_text(size=15,color="black"),
    axis.text.y=element_text(size=12,color="black"),
    #legend.title=element_text(size=12,color="black"),
    legend.title=element_blank(),
    #legend.title.align = 0.5,
    legend.text=element_text(size=10,color="black"),
    legend.background=element_rect(fill = "white", color = "black"),
    legend.key=element_blank(),
    legend.key.size=unit(0.8,"lines"),
    #legend.spacing.y=unit(0.1, "cm"),
    legend.position=c(0.6, 0.88),
    panel.grid.major = element_line(linetype = "dotted", color = "gray"),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    strip.placement = "outside"
  )
}


#########################################
# Define Global Filepaths
#########################################

# list of sites to partition fluxes for:
#site.list <- c("AMF_US-NR1")
site.list <- c("AMF_US-CZ2", "AMF_US-CZ3", "AMF_US-CZ4", "NEO_SJER", "NEO_SOAP", 
               "NEO_TEAK", "AMF_US-NR1")

### Root dir path
caseroot <- "./"

# location of input data 
data.dir <- "/Users/lkunik/Documents/Eddy_Covariance/Kenny_code/1_Processed data/2_Output_from_ProcessingThruREddyProc/1_REddyProc outputs/1_after-night-part/" #failed REddyProc files which have gap-filled NEE, etc
data.dir.compare <- "/Users/lkunik/Library/CloudStorage/Box-Box/FLUXNET synthesis/Kenny/KS_FluxSynth_v8.0/1_Processed data/2_Output_from_ProcessingThruREddyProc/1_REddyProc outputs/1_after-night-part/" #successful REddyProc files for select sites (using for NWR to compare)

# Location of processed data saved during Step 1 (accessed in most scripts that follow)
data.main <- "/Users/lkunik/Documents/Eddy_Covariance/Flux_part/data/"

# Location of processed data saved after flux partitioning
data.case <- paste0(caseroot, "data/") #location of output data

# location of E0 data saved during flux partitioning pass 1
data.E0.case <- paste0(caseroot, "data/E0/")

# location of plot output
plot.dir <- paste0(caseroot, "plots/") #location of output plots

# Define the additional Mike Goulden files that contain PAR (for determining nighttime data)
PAR.files <- c(paste0(caseroot, "include/SolRad_US-CZ2.nc"),
               paste0(caseroot, "include/SolRad_US-CZ3.nc"),
               paste0(caseroot, "include/SolRad_US-CZ4.nc"),
               NA, NA, NA, NA)

names(PAR.files) <- site.list



#########################################
# Define Global Variables and constants
#########################################

NWR_years <- c(2005, 2006, 2007)

# reference and base temperatures outside of all loops:
#Tref_K <- 273.15 + 5 # +5°C. Note, different from Wutzler et al 2018 (15°C) and Reichstein et al 2005 (10°C)

night_PPFD_thresh <- 5 #umol/m2/ms
gf_flag_thresh <- 2 #QC threshold for retaining gap-filled data (NEE*_fqc must be LESS THAN (not ≤) threshold value)

T0_K <- 273.15 - 46.02 #T0 is set to -46.02°C (Lloyd & Taylor, 1994)
#note - reference temp is currently set as the annual mean nighttime temp for each site, not what's below
#Tref_K <- 273.15 + 5 # +5°C. Note, different from Wutzler et al 2018 (15°C) and Reichstein et al 2005 (10°C)

ndays <- 366 # maximum number of days in year
nweeks <- 52 # number of weeks in year

#### Fit routine (determination of E0 and Rref) params
length_window_days_pass1 <- 15 #length of moving window during first pass of fit routine
length_window_days_pass2 <- 7 #length of moving window during second pass of fit routine

resid_low <- 0.05 #lower residual quantile to trim off (for Pass 1)
resid_high <- 0.95 #upper residual quantile to trim off (for Pass 1)
min_records <- 6 #num records (for Pass 1)
min_T_range <- 5 #deg C (for Pass 1)
E0_low <- 30 #Kelvin (for Pass 1 and 2)
E0_high <- 450 #Kelvin

#########################################
# Define Runtime Options
#########################################

plot_debug <- F #Boolean: should we plot the moving 7-day window fit curves as we loop?

