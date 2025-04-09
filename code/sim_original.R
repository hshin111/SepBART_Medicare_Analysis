library("ModBart")
library("SoftBart") # MFM version
# Load a specific version of SoftBART
#library("SoftBart", lib.loc = "/n/home_fasse/kirene/R/ifxrstudio/RELEASE_3_14/00LOCK-SoftBART-MFM/00new/SoftBart/libs")
library("mvtnorm")
library("Rcpp")
library(argparse)

# define parser arguments ----
# args <- list()
# args$year <- 2010
parser <- ArgumentParser()
parser$add_argument("-y", "--year", default=2010, type="integer")  
parser$add_argument("-m", "--model", default="vanilla0", type="character")
parser$add_argument("-w", "--w0_quantile", default=0.25, type="double")    
args = parser$parse_args()

## Call the SepBART function
source('code/original-solution.R')

## Load the data
input_prefix = "../michelle/data/intermediate/SepBART_data/preprocessed_data/all_merged_"
data <- readRDS(paste0(input_prefix, args$year, ".rds"))
cat(paste0("fitting data for year ", args$year))

# colnames <- c("zip", "year", "pm25_ensemble", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", 
#               "medhouseholdincome", "medianhousevalue", "poverty", "education", "popdensity", 
#               "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "pm25", 
#               "ozone", "br", "ca", "cu", "ec", "fe", "k", "nh4", "ni", "no3", "oc", "pb", "si", "so4", 
#               "v", "z", "zcta", "Year", "randall_martin_PM25", "randall_martin_BC", "randall_martin_NH4", 
#               "randall_martin_NIT", "randall_martin_OM", "randall_martin_SO4", "randall_martin_SOIL", 
#               "randall_martin_SS", "total_count", "female_percentage", "dual_percentage", "mean_age", 
#               "death_rate", "percentage_race_labelAsian", "percentage_race_labelBlack", 
#               "percentage_race_labelOther", "percentage_race_labelWhite", 
#               "percentage_race_labelHispanic", "percentage_race_labelNorth_American_Native")

data <- data.frame(lapply(data, as.numeric))
data_prcd <- as.matrix(data)

Y <- data_prcd[,49] # death rate (percentage)

if (args$model == "vanilla0") {
    #### vanilla0 ####

    # original sepbart function 
    # dual, mean_age, pct_white
    # "randall_martin_PM25","ozone"

    X <- data_prcd[, c("dual_percentage", "mean_age", "percentage_race_labelWhite")] 

    #PM2.5, Ozone, EC, OC, NH4, NO3, and SO4
    W <- data_prcd[, c("pm25", "ozone")]
    Mod_ind <- NULL
} else if (args$model == "vanilla1") {
    #### vanilla1 ####

    # original sepbart function,  
    # "randall_martin_PM25","ozone", "ec", "oc", "randall_martin_NH4", "randall_martin_NIT", "randall_martin_SO4"
    # remove: poverty, pop_density, winter_tmp, winter_prp, pct_black
    
    X <- data_prcd[, c("female_percentage", "dual_percentage", "mean_age", "percentage_race_labelWhite", 
    "mean_bmi", "smoke_rate", 
    "medhouseholdincome", "medianhousevalue", "education", "pct_owner_occ", 
    "summer_tmmx", "summer_rmax")] 
    W <- data_prcd[, c("pm25", "ozone", "ec", "oc", "nh4", "no3", "so4")]
    Mod_ind <- NULL
} else if (args$model == "vanilla2") {
    #### vanilla2 ####

    # original sepbart function,  
    # "randall_martin_PM25","ozone", "ec", "oc", "randall_martin_NH4", "randall_martin_NIT", "randall_martin_SO4"
    # remove: poverty, winter_tmp, winter_prp, pct_black
    
    X <- data_prcd[, c("female_percentage", "dual_percentage", "mean_age", "percentage_race_labelWhite", 
    "mean_bmi", "smoke_rate", 
    "medhouseholdincome", "medianhousevalue", "education", "popdensity", "pct_owner_occ", 
    "summer_tmmx", "summer_rmax")] 
    W <- data_prcd[, c("pm25", "ozone", "ec", "oc", "nh4", "no3", "so4")]
    Mod_ind <- NULL
} else if (args$model == "vanilla3") {
    #### vanilla2 ####

    # original sepbart function,  
    # "randall_martin_PM25","ozone", "ec", "oc", "randall_martin_NH4", "randall_martin_NIT", "randall_martin_SO4"
    # remove: poverty, pct_black
    # remove from modifiers: summer_tmp, summer_prp, winter_tmp, winter_prp
    
    X <- data_prcd[, c("female_percentage", "dual_percentage", "mean_age", "percentage_race_labelWhite", 
    "mean_bmi", "smoke_rate", "medhouseholdincome", "medianhousevalue", 
    "education", "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx",
    "summer_rmax", "winter_rmax")] 
    W <- data_prcd[, c("pm25", "ozone", "ec", "oc", "nh4", "no3", "so4")]
    Mod_ind <- c(1:11)
}

w0_quantile = args$w0_quantile

## Fit the model
fit <- SepBART_blk(Y = Y, X = X, W = W,
                   Mod_ind = Mod_ind,
                   W1 = NULL, W0 = NULL, W0_quantile = w0_quantile,
                   nMCMC=5000,BurnIn_portion=0.2,stepsize=8)

# # Define the file path to save the 'fit' object
# if (args$model == "vanilla0") {
#     output_prefix <- paste0("data/output/original/vanilla0/fit_result_")
# } else if (args$model == "vanilla1") {
#     output_prefix <- paste0("data/output/original/vanilla1/fit_result_")
# } else if (args$model == "vanilla2") {
#     output_prefix <- paste0("data/output/original/vanilla2/fit_result_")
# } else if (args$model == "vanilla3") {
#     output_prefix <- paste0("data/output/original/vanilla3/fit_result_")
# }

output_prefix <- paste0("data/output/original/quantile",w0_quantile*100,"/fit_result_")

# Save the 'fit' object as an RDS file
saveRDS(fit, paste0(output_prefix, args$year, ".rds"))

# Print a message to confirm successful saving
cat("The 'fit' object has been saved to:", output_prefix, "\n")
