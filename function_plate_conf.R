################################################################################
#
#Wrap up function to get plate configuration
#CDC protocol
#by: INMEGEN Computational Genomics Dept
#Laura Gómez-Romero lgomez@inmegen.gob.mx
#
################################################################################

################################################################################
#libraries
################################################################################
source("src/functions.R")


################################################################################
#wrap it up in a nice, shiny-compliant function
################################################################################

get_plate_conf <- function(input){
  
  
  my_deltaRN <- tidy_deltaRN(input) #read deltaRN from EDS file 
  
  ########
  #Run function to get plate configuration
  ########
  
  plate_conf <- get_plate(my_deltaRN)
  
  return(plate_conf)
  
}
  