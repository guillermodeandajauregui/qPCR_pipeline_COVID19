################################################################################
#
#qPCR analysis pipeline script for COVID19 detection
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

resultados <- "data/DEMO_COVID2020_THERMO.eds"

################################################################################
#analysis steps
################################################################################

#####1 extract delta RN curves #################################################

my_deltaRN <- tidy_deltaRN(resultados)

##########Plot them for that nice visual inspection feel #######################

pivot_deltaRN(my_deltaRN) %>% 
  plot_deltaRN.long()

  
### TO DO: WHY IS THE Y LABEL DIFFERENT TO THE ONE IN JP's PICTURE? 
### SEEMS LIKE THE VALUES FROM THE FILE ARE MULTIPLIED TIMES 1000?
### CHECK EDS SPECIFICATIONS

#####2 automatically do sanity check of the curves (NOT READY YET) #############

#for i in curves
#check that they are curve?

### TO DO: ASK ANALYTICS TEAM, WHAT ARE PROPER SANITY CHECKS



#####3 Find the threshold (NOT READY YET) ######################################

### TO DO: FIND IN THE EDS SPECIFICATION, WHERE CAN WE FIND THE THRESHOLD

my_threshold <- 10

##### 4 split results by probe (NOT READY YET)  ################################

##### this part to make the data tidy
##### each row is an observation, each column a feature

my_deltaRN.long <- 
my_deltaRN %>% 
  pivot_deltaRN %>% 
  separate(col = sample.id, sep = "_", into = c("well", 
                                                "sample.label", 
                                                "probe"), 
           remove = F)


#### this part to extract the unique values, to use them for iterations

my_wells <- unique(my_deltaRN.long$well)
names(my_wells) <- my_wells

all_probes <- unique(my_deltaRN.long$probe)
names(all_probes) <- all_probes

#####5 Evaluate curves (NOT READY YET) #########################################

### here we iterate over each well
### returning whether the probe passed the threshold or not in that well

my_results <- 
lapply(X = my_wells, FUN = function(my_well){
  
  #get only curves for that well
  well_data <- 
    my_deltaRN.long %>% 
    filter(well == my_well)
  
  
  #we evaluate all probes
  lapply(X = all_probes, FUN = function(my_probe){
    
    #we use a trycatch to get NAs for probes not measured in well
    tryCatch(
      {
        #extract curve
        the_curve     <- extract_curve(well_data, probe == my_probe)
        #extract threshold
        the_threshold <- get_threshold.rg(the_curve)
        print(the_threshold)
        #does the curve crosses the threshold?
        any(the_curve$value > the_threshold)
      },
      error =function(cond){
        message("well does not have that probe")
        return(NA)
      }
    )
  }) %>% bind_rows() #we make them a row
}) %>% bind_rows(.id = "well") #and join them in a data frame

#####5 some clean #########################################



################################################################################
#dragons
################################################################################

#####5 Evaluate curves (NOT READY YET) #########################################

my_results.prep <- 
 
#lapply(X = my_samples, FUN = function(i){ #check with analytics for proper unique id
  lapply(X = my_well_label, FUN = function(i){
    
  temp.df <- 
  my_deltaRN.long %>% 
  #  filter(sample.label == i)
     filter(well_sample.label == i)
  
  #for each probe, check for amplification
  
  
  lapply(X = my_probes, FUN = function(j){
    
    temp.value <- 
    temp.df %>% 
      filter(probe == j) %>% 
      pull(value) %>% 
      max
    
    ##here we make a mock test, whether max curve value > threshold
    ##check what is the actual test that is done visually by analysts
    ##One thing that was discussed, they may want to have an "undetermined" 
    ##to trigger a reanalysis
    
    temp.value > my_threshold
    }) %>% 
    bind_rows()
  
}) %>% bind_rows(.id = "well_sample.label")

##### Finally, we evaluate the conditions based on the rules set in the PNO

#### THIS IS A MOCKUP, WHERE HAVING S == TRUE means POSITIVE for COVID
#### change this logic tests for the actual PNO conditions

my_results <- 
  my_results.prep %>% 
  mutate(diagnosis = ifelse(S, "positive", "negative"))  

################################################################################
#write out 
################################################################################

out_path = "results/prueba.txt"
write.table(my_results, out_path)
