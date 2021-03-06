################################################################################
#
#functions for qPCR analysis - COVID19 detection
#by: INMEGEN Computational Genomics Dept
#Guillermo de Anda J??uregui gdeanda@inmegen.edu.mx
#
################################################################################


################################################################################
#required libs
################################################################################

library(tidyverse)
library("rmarkdown")


################################################################################
#define analysis functions here
################################################################################

##### tidy_deltaRN
tidy_deltaRN <- function(eds){
  #takes an eds file  
  #returns a tidy table with each experiment as column and each cycle as a row
  #a tdrn 
  eds.file <- unz(description = eds, filename = "apldbio/sds/analysis_result.txt")
  analysis_text <- read_lines(eds.file)
  #find the lines with the Delta Rn
  delta_rn.id <- which(str_detect(string = analysis_text, pattern = "Delta Rn"))
  
  #extract those, put them as a vector  
  ### splits them in a list
  list_deltaRN <- analysis_text[delta_rn.id] %>% 
    str_replace(pattern = "Delta Rn values\t", replacement = "") %>% 
    sapply(FUN = function(i){
      str_split(string = i, pattern = "\t") 
    }) 
  
  ### makes them numeric
  list_deltaRN <- 
    lapply(list_deltaRN, FUN = function(i){
      i <- as.numeric(i)
      names(i) <- 1:length(i)
      return(i)
    })
  
  ### name them with the info that is two rows before
  
  nomen <- 
    analysis_text[delta_rn.id - 2] %>% 
    sapply(FUN = function(i){
      intermedio <- str_split(string = i, pattern = "\t")
      intermedio <- unlist(intermedio, recursive = F)
      paste(intermedio[1], intermedio[2], intermedio[3], sep = "_")
    })
  
  names(list_deltaRN) <- nomen
  
  ### make a data frame 
  df_deltaRN <- 
  list_deltaRN %>% 
    bind_rows()
  
  ### add the cycles
  
  df_deltaRN <-
    df_deltaRN %>% 
    mutate(cycles = 1:nrow(df_deltaRN))
  
  return(df_deltaRN)
}


##### pivot_deltaRN

pivot_deltaRN <- function(tdrn){
  #takes a tidy delta RN data frame 
  #pivots it longer -> long_tdrn
  #(syntactic sugar)
  tdrn %>% 
    pivot_longer(cols = -cycles, 
                 names_to = "sample.id", 
                 values_to = "value")
  
}
##### split_longtdrn
split_longtdrn <- function(longtdrn){
  #takes a long_tdrn
  #splits sample.id int well, sample.label, probe columns
  #for maximum tidyness
  longtdrn %>% 
    separate(col = sample.id, sep = "_", into = c("well", 
                                                  "sample.label", 
                                                  "probe"),
             remove = F)
}

##### get_plate
get_plate <- function(tdrn){
  #takes a tdrn 
  #and fills a table with positions
  
  midf <- 
    tdrn %>% 
    pivot_deltaRN %>% 
    split_longtdrn %>% 
    filter(cycles==1) #to have just the one
  
  #make an empty matrix to represent the plates
  mx <- matrix(data = "empty", ncol = 8, nrow = 12)
  rownames(mx) <- 1:12
  colnames(mx) <- LETTERS[1:8]
  
  #loop over used wells to fill the matrix with sample names
  for(i in unique(midf$well)){
    filler <- 
      midf %>% 
      filter(cycles==1) %>% 
      filter(well == i) %>% 
      pull(sample.label) 
    
    mx[as.numeric(i)+1] <- filler
  }
  
  #pivot matrix and present as data frame
  mx <-as.data.frame(t(mx))
  return(mx)
}

##### extract_curve

extract_curve <- function(tdrn_long, ...){
  #takes a tdrn_long, and some filters
  #checks that filters pull only single curve 
  #returns a single curve 
  curve <-
    tdrn_long %>% 
    filter(...) %>% 
    arrange(cycles)
  
  if(nrow(curve)<1){
    stop("filters did not select a curve")
  }
  
  my_cycles <- 
  curve %>% 
    pull(cycles) 
  
  unique_cycles <-
    unique(my_cycles)
  
  if(length(my_cycles) != length(unique_cycles)){
    stop("filters did not identify a unique curve")
  }
  
  return(curve)
}



##### get_threshold.rg

get_threshold.rg <- function(curve){
  #takes a curve, 
  #calculates its threshold value
  #based on 
  #Quant Studio 6 and 7 user manual
  #calculates the mean and standard deviation (SD) of the signal in the cycles 
  #3-10 and places the threshold at mean+10*SD 
  
  mean.ini <- 
  curve %>% 
    filter(cycles%in%3:10) %>% 
    pull(value) %>% 
    mean
  
  sd.ini <- 
    curve %>% 
    filter(cycles%in%3:10) %>% 
    pull(value) %>% 
    sd
  
  my_threshold <- mean.ini + (10 * sd.ini)
  return(my_threshold)

}

##### get_plateThreshold

get_plateThreshold <- function(tdrn_long){
  #sets the baseline as the average value of all wells 
  #in the initial times
  #based on 
  #Quant Studio 6 and 7 user manual
  #calculates the mean and standard deviation (SD) of the signal in the cycles 
  #3-10 and places the threshold at mean+10*SD 
  
  mean.ini <- 
    tdrn_long %>% 
    filter(cycles%in%3:10) %>% 
    pull(value) %>% 
    mean
  
  sd.ini <- 
    tdrn_long %>% 
    filter(cycles%in%3:10) %>% 
    pull(value) %>% 
    sd
  
  my_threshold <- mean.ini + (10 * sd.ini)
  return(my_threshold)
  
}

##### analyze_sample
analyze_sample <- function(tdrn_sample, probes, threshold){
  #takes a filtered tdrn for a single sample
  #and a plate threshold
  #returns a data frame with Ct for each probe
  #returns Inf if value never crosses threshold
  
  lapply(X = probes, FUN = function(my_probe){
    #we use a trycatch to get NAs for probes not measured in well
    tryCatch(
      {
        #extract curve
        the_curve     <- extract_curve(tdrn_sample, probe == my_probe)
        #extract threshold
        
        
        #the_threshold <- get_threshold.rg(the_curve)
        the_threshold = threshold
        #the_threshold <- 725
        #print(the_threshold)
        #does the curve crosses the threshold?
        #threshold for RP should be crossed at time 35
        #any(the_curve$value[1:40] > the_threshold)
        Ct <- 
        the_curve %>% 
          filter(value >= the_threshold) %>% 
          pull(cycles) %>% min
        
      },
      error =function(cond){
        message("well does not have that probe")
        return(NA)
      }
    )
  }) %>% bind_rows()
}

##### plate_qc
plate_qc <- function(tdrn, all_probes){
  #takes an tdrn file 
  #extracts quality control wells
  #analyzes quality control 
  #returns a list with qc table and qc results
  
  #extracts quality control wells
  wells.ntc <- grep(pattern = "_NTC", x = colnames(tdrn))
  wells.ptc <- grep(pattern = "_PTC", x = colnames(tdrn))
  wells.exc <- grep(pattern = "_EC", x = colnames(tdrn))
  ##and filter the tdrn
  qc.df <-
    tdrn %>% 
    select(c(wells.ntc, wells.ptc, wells.exc), cycles) %>% 
    pivot_deltaRN %>% 
    split_longtdrn
  
  #get threshold
  my_threshold <- get_plateThreshold(pivot_deltaRN(tdrn))
  
  
  #name for iteration
  qc.samples <- unique(qc.df$sample.label)
  names(qc.samples) <- qc.samples
  
  #analyze qc
  qc.results <- 
    lapply(qc.samples, FUN = function(my_sample){
      
      sample_data <-
        qc.df %>% 
        filter(sample.label == my_sample) 
      
      #we evaluate all probes
      analyze_sample(tdrn_sample = sample_data, 
                     probes = all_probes, 
                     threshold = my_threshold)
      
    })%>% bind_rows(.id = "sample")  
  
  
  
  #evaluate logic
  
  #check that all probes for NTC DO NOT cross threshold
  
  ntc.all <-
    qc.results %>%
    filter(grepl(pattern = "NTC", x = sample)) %>%
    #select(!sample) %>%
    select(-sample) %>%
    map_dfr(.f = function(i){all(i==Inf)}) %>% #all probes dont cross threshold
    unlist %>% all(. == T) #this should be all true
  
  #check that all probes for PTC DO cross threshold
  ptc.all <-
    qc.results %>%
    filter(grepl(pattern = "PTC", x = sample)) %>%
    #select(!sample)
    select(-sample)
  
  #check that all probes for EC DO NOT cross threshold; same as NTC
  
  ec.all <-
    qc.results %>%
    filter(grepl(pattern = "EC", x = sample)) %>%
    #select(!sample) %>%
    select(-sample) %>%
    map_dfr(.f = function(i){all(i==Inf)}) %>% #all probes dont cross threshold
    unlist %>% all(. == T) #this should be all true
  
  ptc.all <-
    list(RP = all(ptc.all[["RP"]]>=35), #this should be T, do not amplify
         N1 = all(ptc.all[["N1"]]<=38), #this should be T, amplify
         N2 = all(ptc.all[["N2"]]<=38)  #this should be T, amplify
         ) %>% 
    unlist %>% all(. == T) #this should be all true
  
  qc.assess <- ifelse(ntc.all & ptc.all & ec.all, "PASS", "FAIL")
  
  result_final <- list(qc.values = qc.results,
                       ntc.pass = ntc.all,
                       ptc.pass = ptc.all,
                       ec.pass  = ec.all,
                       QC = qc.assess)
  
  if(!ntc.all){print("NTC failed")}
  if(!ptc.all){print("PTC failed")}
  if(!ec.all){print("EX failed")}
  
  return(result_final)

}

##### test.plate

test.plate <- function(tdrn, probes){
  #takes a tdrn
  #returns a table with Ct for each probe, for each sample 
  
  wells.ntc <- grep(pattern = "_NTC", x = colnames(tdrn))
  wells.ptc <- grep(pattern = "_PTC", x = colnames(tdrn))
  wells.exc <- grep(pattern = "_EC", x = colnames(tdrn))
  
  #extract the samples
  test.df <-
    tdrn %>% 
    #select(!c(wells.ntc, wells.ptc), cycles) %>%
    select(-c(wells.ntc, wells.ptc, wells.exc), cycles) %>%
    pivot_deltaRN %>% 
    separate(col = sample.id, sep = "_", into = c("well", 
                                                  "sample.label", 
                                                  "probe"), 
             remove = F)
  
  test.samples <- unique(test.df$sample.label)
  names(test.samples) <- test.samples
  
  #get threshold
  my_threshold <- get_plateThreshold(pivot_deltaRN(tdrn))
  
  #analyze tests
  test.results <- 
    lapply(test.samples, FUN = function(my_sample){
      
      sample_data <-
        test.df %>% 
        filter(sample.label == my_sample) 
      
      #we evaluate all probes
      analyze_sample(tdrn_sample = sample_data, 
                     probes = probes, 
                     threshold = my_threshold)
      
    })%>% bind_rows(.id = "sample")
  
  return(test.results)
  
}


##### cdc_classification
cdc_classification <- function(SampleResults){
  #takes results 
  #assigns classification
  SampleResults %>% 
    mutate(classification = ifelse(test = RP > 35, 
                                   "invalid",
                                   ifelse(N1 <= 38 & N2 <= 38, 
                                          "positive", 
                                          ifelse(N1 > 38 & N2 > 38,
                                                 "negative", 
                                                 "inconclusive"
                                          )
                                   )
    )
    )
}

