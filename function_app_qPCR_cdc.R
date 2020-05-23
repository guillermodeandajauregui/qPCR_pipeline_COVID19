################################################################################
#
#qPCR analysis pipeline script for COVID19 detection
#CDC protocol
#by: INMEGEN Computational Genomics Dept
#Guillermo de Anda J??uregui gdeanda@inmegen.edu.mx
#
################################################################################

################################################################################
#libraries
################################################################################
source("src/functions.R")
source("src/functions_sanitizing.R")
source("src/functions_adjustment.R")
source("src/plots.R")
source("src/reports.R")
################################################################################
#read data
################################################################################

#resultados <- "data/mock_cdc_2probes.eds"

################################################################################
#wrap it up in a nice, shiny-compliant function
################################################################################

qpcr_pipeline.cdc <- function(input, output){
  
#####################
#Sanity checks 
#####################
  
### check eds has results
  
has_results <- CheckResultsEDS(eds = input)
if(has_results == FALSE){
  return("no_results")
}
  
### check sample names in eds have no special characters (only AZaz09 and -)

has_CleanNames <- CheckNamesEDS(eds = input)
if(has_CleanNames == FALSE){
  return("special_characters_in_names")
}
  
my_deltaRN <- tidy_deltaRN(input) #read deltaRN from EDS file 

########
#define probes
########

cdc_probes <- c("RP", "N1", "N2")
names(cdc_probes) <- cdc_probes

### check that all samples have all probes

has_allProbes <- CheckProbesEDS(eds = input, my_probes = cdc_probes)
if(has_allProbes == FALSE){
  return("some_samples_are_missing_probes")
}

################################################################################
#Separate sample wells from QC wells 
################################################################################
#wells.ntc <- grep(pattern = "NTC", x = colnames(my_deltaRN))
#wells.ptc <- grep(pattern = "PTC", x = colnames(my_deltaRN))
################################################################################
#Plate Quality Control
################################################################################

#curve plot
plate_curves <- 
  pivot_deltaRN(my_deltaRN) %>% 
  plot_deltaRN.long()

qc_results <- plate_qc(tdrn = my_deltaRN, 
                       all_probes = cdc_probes)

################################################################################
#Sample analysis
################################################################################


test.results <- test.plate(tdrn = my_deltaRN, probes = cdc_probes)

################################################################################
#Plot preparation
################################################################################

threshold_list = lapply(X = cdc_probes, 
                        FUN = function(i){
  get_probeThreshold(tdrn_long = split_longtdrn(pivot_deltaRN(my_deltaRN)), 
                     my_probe = i)
})

plots.qc <- plot.curves(tdrn = my_deltaRN, 
                        probes = cdc_probes, 
                        threshold_list = threshold_list,
                        qc = T)

plots.samples <- plot.curves(tdrn = my_deltaRN, 
                             probes = cdc_probes, 
                             threshold_list = threshold_list,
                             qc = F)

triplets.qc <- triplets(plots.qc)

triplets.samples <- triplets(plots.samples)

################################################################################
#Write individual reports
################################################################################

make_reports(plot_list = triplets.samples, 
             result_table = qc_results$qc.values, 
             input = input,
             outdir = output, 
             qc_results = qc_results$QC,
             qc = F)

make_reports(plot_list = triplets.qc, 
             result_table = qc_results$qc.values, 
             input = input,
             outdir = output, 
             qc_results = qc_results$QC,
             qc = T)

################################################################################
#Get the plate name
################################################################################


plate <- stringr::str_remove(string = basename(input), pattern = ".eds")


################################################################################
#Sample diagnostic
################################################################################

test_diagnosis <- cdc_classification(test.results)

test_diagnosis <- 
data.frame(plate    = plate,
           ntc.pass = qc_results$ntc.pass, 
           ptc.pass = qc_results$ptc.pass, 
           ec.pass = qc_results$ec.pass,
           qc      = qc_results$QC,
           test_diagnosis 
) %>% as.tbl() %>% select(sample, everything())

################################################################################
#Write plate Booklet
################################################################################

plateBooklet(results = test_diagnosis,
             qc_results = qc_results,
             outdir = output)

################################################################################
#Write QC output
################################################################################

out_path_qc = paste(output, "/", plate, ".qc", sep="")
write_delim(x = qc_results$qc.values, path = out_path_qc)
################################################################################
#Write run output
################################################################################

out_path_results = paste(output, "/", plate, ".res", sep="")
write_delim(x = test_diagnosis, path = out_path_results)
################################################################################
#Done
################################################################################

################################################################################
#Create list output
################################################################################

results_list <- list(
  test_results = test_diagnosis,
  qc_results = qc_results, 
  plots_qc = plots.qc, 
  plots_samples = plots.samples,
  triplets.qc = triplets.qc, 
  triplets.samples = triplets.samples
)

return(results_list)
}
