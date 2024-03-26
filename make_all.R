
# -------------------------------------------------------------------------------
# code for analysing 18 benthic state indicators calculated for 14 gradients of 
# commercial bottom trawling intensity, one eutrophication, one oxygen depletion,
# and one pollution
#
# Contact: Daniel van Denderen (pdvd@aqua.dtu.dk)
# -------------------------------------------------------------------------------

# create output directory
outdir <- "Output"
if (! dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)

# obtain the data used
indic    <- read.csv("Data/Indicator_outputs.csv")
stations <- read.csv("Data/Station_information.csv")

# make sure all packages are installed and source libraries
source("Analysis/Libraries_analysis.R")

# source correlation plot
source("Analysis/Correlation_plot.R")

# source meta-analysis
source("Analysis/Meanresponse_metafor.R")

# source overview plot
source("Analysis/Overview_plot.R")

# source individual gradients
source("Analysis/Trawl_gradient_plots.R")

# source other pressure gradients
source("Analysis/Nontrawl_gradient_plots.R")
