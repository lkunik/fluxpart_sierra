# Run all scripts in sequence

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))

#source("1_prep_data.r")
source("2_get_E0_vals.r")
source("3_flux_part.r")
source("4_combine_years.r")
source("5_plot_Rref_E0.r")
source("6_compare_NWR.r")
