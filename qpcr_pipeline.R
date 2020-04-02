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
  plot_deltaRN.long(y_title = "probando")

  
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
  pivot_longer(cols = -cycles, 
               names_to = "sample.id", 
               values_to = "value") %>% 
  separate(col = sample.id, sep = "_", into = c("well", 
                                                "sample.label", 
                                                "probe"), 
           remove = F)


#### this part to extract the unique values, to use them for iterations

my_probes <- unique(my_deltaRN.long$probe)
names(my_probes) <- my_probes

my_samples <- unique(my_deltaRN.long$sample.label)
names(my_samples) <- my_samples

###QUESTION: ARE SAMPLE.LABELS UNIQUE? 
#OR IS A SAMPLE ADDED TO MORE THAN ONE WELL?
#FOR NOW: Make them unique adding well_sample.label

my_deltaRN.long <- 
  my_deltaRN.long %>% 
  mutate(well_sample.label = paste0(well, "_", sample.label))

my_well_label <- unique(my_deltaRN.long$well_sample.label)
names(my_well_label) <- my_well_label

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
