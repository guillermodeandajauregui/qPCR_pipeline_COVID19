################################################################################
#
#functions for qPCR analysis - COVID19 detection
#by: INMEGEN Computational Genomics Dept
#Guillermo de Anda Jáuregui gdeanda@inmegen.edu.mx
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
##### split_tidyRN.long
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
##### plot_deltaRN.long 

plot_deltaRN.long <- function(tdrn_long, 
                              guide_title = "muestra", 
                              y_title = "Delta_RN"){
  #takes a long_tdrn
  #plots all curves in it, grouped by sample.id
  #(syntactic sugar)
  tdrn_long %>% 
    ggplot(mapping = aes(x = cycles, 
                         y = value, 
                         colour = as.factor(sample.id)
    )
    ) + 
    geom_line() +
    guides(color = guide_legend(title = guide_title)) +
    ylab(y_title) +
    theme_minimal()
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

##### analyze_sample
analyze_sample <- function(tdrn_sample, probes){
  #takes a filtered tdrn for a single sample
  #returns a data frame with Ct for each probe
  #returns Inf if value never crosses threshold
  
  lapply(X = probes, FUN = function(my_probe){
    #we use a trycatch to get NAs for probes not measured in well
    tryCatch(
      {
        #extract curve
        the_curve     <- extract_curve(tdrn_sample, probe == my_probe)
        #extract threshold
        
        
        the_threshold <- get_threshold.rg(the_curve)
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
  wells.ntc <- grep(pattern = "NTC", x = colnames(tdrn))
  wells.ptc <- grep(pattern = "PTC", x = colnames(tdrn))
  
  ##and filter the tdrn
  qc.df <-
    tdrn %>% 
    select(c(wells.ntc, wells.ptc), cycles) %>% 
    pivot_deltaRN %>% 
    split_tidyRN.long
  
  
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
      analyze_sample(tdrn_sample = sample_data, probes = all_probes)
      
    })%>% bind_rows(.id = "sample")  
  
  
  
  #evaluate logic
  
  #check that all probes for NTC DO NOT cross threshold
  
  ntc.all <-
    qc.results %>%
    filter(grepl(pattern = "NTC", x = sample)) %>%
    select(!sample) %>%
    map_dfr(.f = function(i){all(i==Inf)}) %>% #all probes dont cross threshold
    unlist %>% all(. == T) #this should be all true
  
  #check that all probes for PTC DO cross threshold
  ptc.all <-
    qc.results %>%
    filter(grepl(pattern = "PTC", x = sample)) %>%
    select(!sample)
  
  ptc.all <-
    list(RP = all(ptc.all[["RP"]]<=35), #this should be T
         N1 = all(ptc.all[["N1"]]<=40), #this should be T
         N2 = all(ptc.all[["N2"]]<=40)  #this should be T
         ) %>% 
    unlist %>% all(. == T) #this should be all true
  
  
  qc.assess <- ifelse(ntc.all & ptc.all, "PASS", "FAIL")
  
  result_final <- list(qc.values = qc.results,
                       ntc.pass = ntc.all,
                       ptc.pass = ptc.all, 
                       QC = qc.assess)
  
  if(!ntc.all){print("NTC failed")}
  if(!ptc.all){print("PTC failed")}
  
  return(result_final)

}

##### test.plate

test.plate <- function(tdrn, probes){
  #takes a tdrn
  #returns a table with Ct for each probe, for each sample 
  
  wells.ntc <- grep(pattern = "NTC", x = colnames(my_deltaRN))
  wells.ptc <- grep(pattern = "PTC", x = colnames(my_deltaRN))
  
  #extract the samples
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
  
  #analyze tests
  test.results <- 
    lapply(test.samples, FUN = function(my_sample){
      
      sample_data <-
        test.df %>% 
        filter(sample.label == my_sample) 
      
      #we evaluate all probes
      analyze_sample(tdrn_sample = sample_data, probes = probes)
      
    })%>% bind_rows(.id = "sample")
  
  return(test.results)
  
}

##### plot.curves
plot.curves <- function(tdrn, probes, qc = TRUE){
  #takes an tdrn file
  #makes plots for either qc wells
  #or non-qc wells
  
  #extracts quality control wells
  wells.ntc <- grep(pattern = "NTC", x = colnames(tdrn))
  wells.ptc <- grep(pattern = "PTC", x = colnames(tdrn))
  
  ##and filter the tdrn
  if(qc == TRUE){
    qc.df <-
      tdrn %>% 
      select(c(wells.ntc, wells.ptc), cycles) %>% 
      pivot_deltaRN %>% 
      split_tidyRN.long
  }else{
    qc.df <-
      tdrn %>% 
      select(!c(wells.ntc, wells.ptc), cycles) %>% 
      pivot_deltaRN %>% 
      split_tidyRN.long
  }
  
  
  #name for iteration
  qc.samples <- unique(qc.df$sample.label)
  names(qc.samples) <- qc.samples
  
  #analyze qc
  qc.results <- 
    lapply(qc.samples, FUN = function(my_sample){
      
      sample_data <-
        qc.df %>% 
        filter(sample.label == my_sample) 
      
      lapply(X = probes, FUN = function(my_probe){
        
        the_curve     <- extract_curve(sample_data, probe == my_probe)
        #extract threshold
        
        
        the_threshold <- get_threshold.rg(the_curve)
        
        p <- 
          plot_deltaRN.long(tdrn_long = the_curve) + 
          geom_hline(yintercept = the_threshold, linetype = 2)
        
      })
      
      
    })
}

##### cdc_classification
cdc_classification <- function(SampleResults){
  #takes results 
  #assigns classification
  SampleResults %>% 
    mutate(classification = ifelse(test = RP > 35, 
                                   "invalid",
                                   ifelse(N1 <= 40 & N2 <= 40, 
                                          "positive", 
                                          ifelse(N1 > 40 & N2 > 40,
                                                 "negative", 
                                                 "inconclusive"
                                          )
                                   )
    )
    )
}


######make_reports

make_reports <- function(plot_list, 
                         result_table,
                         outdir, 
                         #qc.result,
                         qc = F){
  #makes reports from a list of plots and some result table
  lapply(seq_along(plot_list), function(i){
    
    the_sample_is <- names(plot_list)[i]
    
    my_r <- 
      result_table %>% 
      filter(sample == the_sample_is)
    
    my_name <- names(plot_list)[i]
    mea_plote <- plot_list[[i]]
    
    outpath <- paste0(outdir, "/", Sys.Date(), "_", my_name, ".pdf")
    if(qc==F){
      render("template.Rmd",output_file = outpath)
    }else{
      render("template_qc.Rmd",output_file = outpath)
    }
  })
}
