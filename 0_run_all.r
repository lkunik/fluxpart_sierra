# Run all scripts in sequence

# Set working directory to the path of the calling script (won't work if running
# line by line but will work if script is sourced)
setwd(dirname(sys.frame(1)$ofile))

# Set toggle - plot the extra goodies? T/F


source("1_prep_data.r")
source("2_get_E0_vals.r")
source("3_flux_part.r")
source("4_combine_years.r")
source("5_plot_Rref_E0.r")
source("6_compare_NWR.r")

plot_extras <- T

if(plot_extras){
  source("7_light_response.r")
  source("8a_plot_extra_flux_ts.r")
  source("8b_plot_extra_flux_doy.r")
  source("8c_plot_extra_goulden_fig.r")
  source("8d_plot_fluxpart_figure.r")
  source("8e_plot_fluxpart_figure_v2.R")
  # source("8f_plot_fluxpart_compare_flags.r")
}
