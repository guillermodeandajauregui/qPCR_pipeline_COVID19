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
source("src/plots.R")
################################################################################
#read data
################################################################################

#resultados <- "data/mock_cdc_2probes.eds"

################################################################################
#wrap it up in a nice, shiny-compliant function
################################################################################

qpcr_pipeline.cdc <- function(input, output){
  

my_deltaRN <- tidy_deltaRN(input) #read deltaRN from EDS file 

########
#define probes
########

cdc_probes <- c("RP", "N1", "N2")
names(cdc_probes) <- cdc_probes
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

plots.qc <- plot.curves(tdrn = my_deltaRN, 
                        probes = cdc_probes, 
                        qc = T)

plots.samples <- plot.curves(tdrn = my_deltaRN, 
                             probes = cdc_probes, 
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
#Sample diagnostic
################################################################################

test_diagnosis <- cdc_classification(test.results)

################################################################################
#Write QC output
################################################################################

out_path_qc = paste0(output, "qc_results.txt", collapse="/")
write_delim(x = qc_results$qc.values, path = out_path_qc)
################################################################################
#Write run output
################################################################################

out_path_results = paste0(output, "test.results.txt", collapse="/")
write_delim(x = test.results, path = out_path_results)
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
