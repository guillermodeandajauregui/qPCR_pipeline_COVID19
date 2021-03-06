################################################################################
#
#qPCR analysis pipeline script for COVID19 detection
#CDC protocol
#by: INMEGEN Computational Genomics Dept
#Guillermo de Anda Jáuregui gdeanda@inmegen.edu.mx
#
################################################################################

################################################################################
#libraries
################################################################################
source("src/functions.R")

################################################################################
#read data
################################################################################

resultados <- "data/mock_cdc_3.eds"

my_deltaRN <- tidy_deltaRN(resultados) #read deltaRN from EDS file 

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
################################################################################
#Sample diagnostic
################################################################################

test_diagnosis <- cdc_classification(test.results)

################################################################################
#Write individual reports
################################################################################

make_reports(plot_list = plots.qc, 
             result_table = qc_results$qc.values, 
             outdir = "results/otra_prueba/", 
             qc = T)


################################################################################
#Write QC output
################################################################################

write_delim(x = qc_results$qc.values, path = "results/qc_results.txt")
################################################################################
#Write run output
################################################################################

write_delim(x = test.results, path = "results/test.results.txt")
################################################################################
#Done
################################################################################
