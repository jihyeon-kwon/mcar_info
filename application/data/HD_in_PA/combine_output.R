####################
# combined results #
####################




# setup -----

library(here)
source(here("fun/combine_results_functions.R"))


# combine the summary -----

load(here("application/data/hd-35-64.Rdata"))
suff <- label_sexrace

path <- here("application/output/pa_hd_mcar_ar_res_35-64_")

combine_summary(path, suff, type = "mcar",
                save_to = here("application/output/combined/pa_hd_mcar_ar_res_35-64.rdata"))



# add mode to the summary files -----

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

path <- here("application/output/pa_hd_mcar_ar_res_35-64_")
combined_info <- combine_info(path, suff)
mode <- apply(combined_info$info, c(1,3), estimate_mode) # get mode
load(here("application/output/combined/pa_hd_mcar_ar_res_35-64.rdata"))
combined_summary$info_new <- array(dim = (dim(combined_summary$info) + c(1, 0, 0)))
combined_summary$info_new[6, , ] <- mode
combined_summary$info_new[1:5, , ] <- combined_summary$info
dimnames(combined_summary$info_new)[[1]] <- c("2.5%", "50%", "97.5%", "mean", "sd", "mode")      
dimnames(combined_summary$info_new)[[2]] <- dimnames(combined_summary$info)[[2]]
dimnames(combined_summary$info_new)[[3]] <- dimnames(combined_summary$info)[[3]]
combined_summary$info <- combined_summary$info_new
save(combined_summary, file = here("application/output/combined/pa_hd_mcar_ar_res_35-64.rdata"))
