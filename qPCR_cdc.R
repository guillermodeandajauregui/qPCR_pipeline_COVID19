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

resultados <- "data/mock_cdc_2probes.eds"

my_deltaRN <- tidy_deltaRN(resultados) #read deltaRN from EDS file 

########
#define probes
########

all_probes <- c("RP", "N1", "N2")
names(all_probes) <- all_probes
################################################################################
#Separate sample wells from QC wells 
################################################################################
wells.ntc <- grep(pattern = "NTC", x = colnames(my_deltaRN))
wells.ptc <- grep(pattern = "PTC", x = colnames(my_deltaRN))
################################################################################
#Plate Quality Control
################################################################################

#curve plot
plate_curves <- 
pivot_deltaRN(my_deltaRN) %>% 
  plot_deltaRN.long()

qc_results <- plate_qc(my_deltaRN)

################################################################################
#Sample analysis
################################################################################
test.df <-
  my_deltaRN %>% 
  select(!c(wells.ntc, wells.ptc), cycles) %>% 
  pivot_deltaRN %>% 
  separate(col = sample.id, sep = "_", into = c("well", 
                                                "sample.label", 
                                                "probe"), 
           remove = F)

test.samples <- unique(test.df$sample.label)
names(test.samples) <- test.samples

test.results <- 
  lapply(test.samples, FUN = function(my_sample){
    
    sample_data <-
      test.df %>% 
      filter(sample.label == my_sample) 
    
    #we evaluate all probes
    
    lapply(X = all_probes, FUN = function(my_probe){
      #we use a trycatch to get NAs for probes not measured in well
      tryCatch(
        {
          #extract curve
          the_curve     <- extract_curve(sample_data, probe == my_probe)
          #extract threshold
          
          
          
          the_threshold <- get_threshold.rg(the_curve)
          print(the_threshold)
          #does the curve crosses the threshold?
          #threshold for RP should be crossed at time 35
          any(the_curve$value[1:40] > the_threshold)
        },
        error =function(cond){
          message("well does not have that probe")
          return(NA)
        }
      )
    }) %>% bind_rows()
  })%>% bind_rows(.id = "sample")

test.results  
################################################################################
#Plot preparation
################################################################################

################################################################################
#Sample diagnostic
################################################################################

################################################################################
#Write individual reports
################################################################################

################################################################################
#Write QC output
################################################################################

################################################################################
#Write run output
################################################################################

################################################################################
#Done
################################################################################
