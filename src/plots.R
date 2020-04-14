library("cowplot")



##### plot.curves
plot.curves <- function(tdrn, probes, qc = TRUE){
  #takes an tdrn file
  #makes plots for either qc wells
  #or non-qc wells
  
  #extracts quality control wells
  wells.ntc <- grep(pattern = "_NTC", x = colnames(tdrn))
  wells.ptc <- grep(pattern = "_PTC", x = colnames(tdrn))
  wells.exc <- grep(pattern = "_EC", x = colnames(tdrn))
  
  ##and filter the tdrn
  if(qc == TRUE){
    qc.df <-
      tdrn %>% 
      select(c(wells.ntc, wells.ptc, wells.exc), cycles) %>% 
      pivot_deltaRN %>% 
      split_longtdrn
  }else{
    qc.df <-
      tdrn %>% 
      #select(!c(wells.ntc, wells.ptc), cycles) %>%
      select(-c(wells.ntc, wells.ptc, wells.exc), cycles) %>%
      pivot_deltaRN %>% 
      split_longtdrn
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
        
        color <- ifelse(unique(the_curve$probe) == "RP", "#e41a1c", ifelse(unique(the_curve$probe) == "N1", "#377eb8", "#4daf4a"))
        
        the_threshold <- get_threshold.rg(the_curve)
        
        p <- 
          plot_deltaRN.long(tdrn_long = the_curve) + 
          geom_line(colour = color) +
          geom_hline(yintercept = the_threshold, linetype = 2)
        
      })
      
      
    })
}

