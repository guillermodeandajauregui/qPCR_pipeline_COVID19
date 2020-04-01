################################################################################
#
#functions for qPCR analysis - COVID19 detection
#by: INMEGEN Computational Genomics Dept
#Guillermo de Anda Jáuregui gdeanda@inmegen.edu.mx
#
################################################################################

################################################################################
#define analysis functions here
################################################################################

tidy_deltaRN <- function(analysis_text){
  #takes an "analysis_result.txt" from an eds 
  #returns a tidy table with each experiment as column and each cycle as a row
  
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

tidy_deltaRN(resultados)$cycles
