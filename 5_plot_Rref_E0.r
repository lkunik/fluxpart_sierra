rm(list=ls())

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))

# Load packages, global functions and variables included in config script
source("config.r")


#########################################
# Begin main
#########################################

# create empty lists of ggplot objects ("grobs") to add to within loops
grob.list.E0 <- vector(mode = "list", length = length(site.list))
grob.list.Rref <- vector(mode = "list", length = length(site.list))
grob.list.NEE <- vector(mode = "list", length = length(site.list))
grob.list.Tair <- vector(mode = "list", length = length(site.list))
names(grob.list.E0) <- site.list
names(grob.list.Rref) <- site.list
names(grob.list.NEE) <- site.list
names(grob.list.Tair) <- site.list


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
    
    flux.df <- readRDS(file)
    site.df %<>%
      bind_rows(flux.df)
    
  }
  
  # Create out of bounds (OOB) columns in site.df
  E0_min <- -500
  E0_max <- 500
  E0_OOB <- site.df$E0_U50
  E0_OOB[E0_OOB > E0_min & E0_OOB < E0_max] <- NA
  E0_OOB[E0_OOB <= E0_min] <- E0_min + 0.0001
  E0_OOB[E0_OOB >= E0_max] <- E0_max - 0.0001
  site.df %<>% mutate(E0_U50_OOB = E0_OOB)
  
  perc_E0_NA <- round(length(which(is.na(site.df$E0_U50)))*100/length(site.df$E0_U50), 2)
  message(paste0(perc_E0_NA, "% of all E0 daily values are NA"))

  
  ### Plot E0
  ylims <- c(max(min(site.df$E0_U50, na.rm = T), -500), min(max(site.df$E0_U50, na.rm = T), 500))
  # plot data
  grob <- site.df %>% # take data
    ggplot() + theme_FS() + # pipe in ggplot call
    geom_hline(aes(yintercept = 0), size = 0.5, linetype = "dashed", color = "darkslateblue", alpha = 0.7) +
    geom_point(aes(x=posix_time, y=E0_U50), size=0.5) + # plot E0 vs time
    geom_point(aes(x=posix_time, y=E0_U50_OOB), size=0.5, color = "red") +
    ylim(ylims) + #scale_y_log10() +
    xlab("") + scale_x_datetime(labels = label_date("%b %Y")) + #label should be e.g. "Jan 2013"
    ylab("") +
    ggtitle(site)
  grob.list.E0[[site]] <- grob #assign to the corresponding list of ggplot objects ("grobs")
  
  ### Plot Rref

  perc_Rref_NA <- round(length(which(is.na(site.df$Rref_U50)))*100/length(site.df$Rref_U50), 2)
  message(paste0(perc_Rref_NA, "% of all Rref daily values are NA\n"))
  
  # Create out of bounds (OOB) column in site.df
  Rref_min <- -50
  Rref_max <- 50
  Rref_OOB <- site.df$Rref_U50
  Rref_OOB[Rref_OOB > Rref_min & Rref_OOB < Rref_max] <- NA
  Rref_OOB[Rref_OOB <= Rref_min] <- Rref_min + 0.0001
  Rref_OOB[Rref_OOB >= Rref_max] <- Rref_max - 0.0001
  site.df %<>% mutate(Rref_U50_OOB = Rref_OOB)
  
  ylims <- c(max(min(site.df$Rref_U50, na.rm = T), Rref_min), min(max(site.df$Rref_U50, na.rm = T), Rref_max))
  # plot data
  grob <- site.df %>% # take data
    ggplot() + theme_FS() + # pipe in ggplot call
    geom_hline(aes(yintercept = 0), size = 0.5, linetype = "dashed", color = "darkslateblue", alpha = 0.7) +
    geom_point(aes(x=posix_time, y=Rref_U50), size=0.5) + # plot Rref vs time
    geom_point(aes(x=posix_time, y=Rref_U50_OOB), size=0.5, color = "red") +
    ylim(ylims) + #scale_y_log10() +
    xlab("") + scale_x_datetime(labels = label_date("%b %Y")) + #label should be e.g. "Jan 2013"
    ylab("") +
    ggtitle(site)
  grob.list.Rref[[site]] <- grob #assign to the corresponding list of ggplot objects ("grobs")
  
  
  ### Plot NEE
  
  NEE_min <- -20
  NEE_max <- 20
  NEE_OOB <- (site.df %>% filter(daynight == "night"))$NEE_U50_f
  NEE_OOB[NEE_OOB > NEE_min & NEE_OOB < NEE_max] <- NA
  NEE_OOB[NEE_OOB <= NEE_min] <- NEE_min + 0.0001
  NEE_OOB[NEE_OOB >= NEE_max] <- NEE_max - 0.0001
  site.df %<>% filter(daynight == "night") %>%
    mutate(NEE_U50_OOB = NEE_OOB) %>%
    filter(NEE_U50_fqc < 2)
  
  ylims <- c(NEE_min, NEE_max)
  # plot data
  grob <- site.df %>% # take data
    filter(daynight == "night", NEE_U50_fqc < 2) %>%
    ggplot() + theme_FS() + # pipe in ggplot call
    geom_hline(aes(yintercept = 0), size = 0.5, linetype = "dashed", color = "darkslateblue", alpha = 0.7) +
    geom_point(aes(x=posix_time, y=NEE_U50_f), size=0.5) + # plot NEE vs time
    geom_point(aes(x=posix_time, y=NEE_U50_OOB), size=0.5, color = "red") +
    ylim(ylims) + #scale_y_log10() +
    xlab("") + scale_x_datetime(labels = label_date("%b %Y")) + #label should be e.g. "Jan 2013"
    ylab("") +
    ggtitle(site)
  grob.list.NEE[[site]] <- grob #assign to the corresponding list of ggplot objects ("grobs")
  
  
  ### Plot Tair
  
  Tair_min <- -20
  Tair_max <- 35
  Tair_OOB <- (site.df %>% filter(daynight == "night"))$Tair_f
  Tair_OOB[Tair_OOB > Tair_min & Tair_OOB < Tair_max] <- NA
  Tair_OOB[Tair_OOB <= Tair_min] <- Tair_min + 0.0001
  Tair_OOB[Tair_OOB >= Tair_max] <- Tair_max - 0.0001
  site.df %<>% filter(daynight == "night") %>%
    mutate(Tair_OOB = Tair_OOB) %>%
    filter(Tair_fqc < 2)
  
  ylims <- c(Tair_min, Tair_max)
  # plot data
  grob <- site.df %>% # take data
    filter(daynight == "night", Tair_fqc < 2) %>%
    ggplot() + theme_FS() + # pipe in ggplot call
    geom_hline(aes(yintercept = 0), size = 0.5, linetype = "dashed", color = "darkslateblue", alpha = 0.7) +
    geom_point(aes(x=posix_time, y=Tair_f), size=0.5) + # plot Tair vs time
    geom_point(aes(x=posix_time, y=Tair_OOB), size=0.5, color = "red") +
    ylim(ylims) + #scale_y_log10() +
    xlab("") + scale_x_datetime(labels = label_date("%b %Y")) + #label should be e.g. "Jan 2013"
    ylab("") +
    ggtitle(site)
  grob.list.Tair[[site]] <- grob #assign to the corresponding list of ggplot objects ("grobs")
  
  
} #end sites for-loop

gridRows <- ceiling(length(site.list)/2)
gridCols <- 2


### make multi-panel plot of E0 v time for all sites
main.title <- "E0, All available years, 50% u* threshold"
yaxis.lab <-  "[Kelvin]" #units are K

message("Printing multi-panel E0 plot to PDF")

#print PDF
pdf_outfile <- paste0(plot.dir, "E0_all_sites_all_years.pdf") #define plot filename
pdf(pdf_outfile, width = 12, height = 16) #open pdf
grid.arrange(grobs = grob.list.E0, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf

#print PNG
message("Printing multi-panel E0 plot to PNG")

png_outfile <- paste0(plot.dir, "E0_all_sites_all_years.png") #define plot filename
png(png_outfile, width = 800, height = 1000) #open pdf
grid.arrange(grobs = grob.list.E0, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf


### make multi-panel plot of Rref v time for all sites

main.title <- "R_ref, All available years, 50% u* threshold"
yaxis.lab <-  "[umol m-2 s-1]" #units are K

message("Printing multi-panel Rref plot to PDF")

#print PDF
pdf_outfile <- paste0(plot.dir, "Rref_all_sites_all_years.pdf") #define plot filename
pdf(pdf_outfile, width = 12, height = 16) #open pdf
grid.arrange(grobs = grob.list.Rref, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf

#print PNG

message("Printing multi-panel Rref plot to PNG")

png_outfile <- paste0(plot.dir, "Rref_all_sites_all_years.png") #define plot filename
png(png_outfile, width = 800, height = 1000) #open pdf
grid.arrange(grobs = grob.list.Rref, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf


### make multi-panel plot of NEE v nighttime for all sites
main.title <- "Nighttime NEE (PPFD < 5), All available years, 50% u* threshold"
yaxis.lab <-  "[umol m-2 s-1]" #units are K

message("Printing multi-panel NEE plot to PDF")

#print PDF
pdf_outfile <- paste0(plot.dir, "NEE_night_all_sites_all_years.pdf") #define plot filename
pdf(pdf_outfile, width = 12, height = 16) #open pdf
grid.arrange(grobs = grob.list.NEE, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf

#print PNG

message("Printing multi-panel NEE plot to PNG")

png_outfile <- paste0(plot.dir, "NEE_night_all_sites_all_years.png") #define plot filename
png(png_outfile, width = 800, height = 1000) #open pdf
grid.arrange(grobs = grob.list.NEE, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf




### make multi-panel plot of Tair v nighttime for all sites
main.title <- "Nighttime Tair (PPFD < 5), All available years, 50% u* threshold"
yaxis.lab <-  "[°Celsius]" #units are deg C

message("Printing multi-panel Tair plot to PDF")

#print PDF
pdf_outfile <- paste0(plot.dir, "Tair_night_all_sites_all_years.pdf") #define plot filename
pdf(pdf_outfile, width = 12, height = 16) #open pdf
grid.arrange(grobs = grob.list.Tair, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf

#print PNG

message("Printing multi-panel Tair plot to PNG")

png_outfile <- paste0(plot.dir, "Tair_night_all_sites_all_years.png") #define plot filename
png(png_outfile, width = 800, height = 1000) #open pdf
grid.arrange(grobs = grob.list.Tair, nrow = gridRows, ncols = gridCols,
             top = textGrob(main.title, gp = gpar(fontsize = 20)),
             left = textGrob(yaxis.lab, gp = gpar(fontsize = 20), rot = 90))
dev.off() #close pdf