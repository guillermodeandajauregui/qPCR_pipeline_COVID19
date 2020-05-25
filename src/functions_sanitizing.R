################################################################################
#
#functions for qPCR analysis - COVID19 detection
#Sanitizing inputs
#by: INMEGEN Computational Genomics Dept
#Guillermo de Anda Jauregui gdeanda@inmegen.edu.mx
#
################################################################################

#for ease of convention, these functions evaluate 
#TRUE  if they pass the test
#FALSE  if they fail the test

###CheckResultsEDS 
CheckResultsEDS <- function(eds){
  #takes a path to an eds 
  #if it can read analysis_result, returns TRUE, else FALSE
  eds.file <- try(unz(description = eds, filename = "apldbio/sds/analysis_result.txt"), silent = T)
  #if(class(eds.file)=="try-error"){
  #  message("the error is at the unz level")
  #  return(FALSE)
  #}
  
  analysis_text <- try(suppressWarnings(read_lines(eds.file)), silent = T)
  #print(analysis_text)
  if(class(analysis_text)=="try-error"){
    message("EDS does not contain analysis_result.txt file")
    return(FALSE)
  }
  
  return(TRUE)
}

###CheckNamesEDS

CheckNamesEDS <- function(eds){
  #takes a path to an eds 
  #extracts and checks that sample names are legible 
  eds.file <- unz(description = eds, filename = "apldbio/sds/analysis_result.txt")
  analysis_text <- read_lines(eds.file)
  #find the lines with the Delta Rn
  delta_rn.id <- which(str_detect(string = analysis_text, pattern = "Delta Rn"))
  
  ### name them with the info that is two rows before
  
  nomen <- 
    analysis_text[delta_rn.id - 2] %>% 
    sapply(FUN = function(i){
      intermedio <- str_split(string = i, pattern = "\t")
      intermedio <- unlist(intermedio, recursive = F)
      intermedio[2]
    })
  
  
  #test that the sample id has no special characters
  test_for_special = grepl(pattern = "[^a-zA-Z0-9-]", x = nomen)
  my_r = all(test_for_special == FALSE) #if there is one sample with bad names, it evaluates FALSE
  #pass is TRUE
  #fail is FALSE
  return(my_r)
  
}

### Check correct probes 

CheckProbesEDS <- function(eds, my_probes){
  #takes a path to an eds 
  #extracts and checks that sample names are legible 
  eds.file <- unz(description = eds, filename = "apldbio/sds/analysis_result.txt")
  analysis_text <- read_lines(eds.file)
  #find the lines with the Delta Rn
  delta_rn.id <- which(str_detect(string = analysis_text, pattern = "Delta Rn"))
  
  ### name them with the info that is two rows before
  
  nomen <- 
    analysis_text[delta_rn.id - 2] %>% 
    lapply(FUN = function(i){
      intermedio <- str_split(string = i, pattern = "\t")
      intermedio <- unlist(intermedio, recursive = F)
      df = data.frame(sample = intermedio[2], probe = intermedio[3])
    }) 
  
  nomen = suppressWarnings(bind_rows(nomen))
  
  los_probes = 
  nomen %>% 
    group_by(sample) %>% 
    group_map(~ unique(.x$probe))
  
  #check they are the same probes with jaccard
  j_probes <-
  lapply(los_probes, function(i){
    
    J = length(intersect(i, my_probes))/length(union(i, my_probes))
    
    return(J==1)
  })
  
  #if all samples have all probes, TRUE
  resultado = suppressWarnings(all(j_probes))
  return(resultado)
  
}

