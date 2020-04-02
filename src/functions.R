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
    filter(...)
  
  my_cycles <- 
  curve %>% 
    pull(cycles) 
  
  unique_cycles <-
    unique(my_cycles)
  
  if(length(my_cycles) != length(unique_cycles)){
    stop("filters do not identify a unique curve")
  }
  
  return(curve)
}



##### get_threshold.rg

get_threshold.rg <- function(curve){
  #takes a curve, 
  #calculates its threshold value
  #based on 
  #https://www.researchgate.net/post/How_can_I_set_the_threshold_in_a_Real-Time_PCR_result
  #which says: 
  #I know that the SDS (for PerkinElmer/Applied Biosystems) 
  #calculates the mean and standard deviation (SD) of the signal in the cycles 
  #3-10 and places the threshold at mean+10*SD 
  #(often not quite sensible, though!).
  
  
}