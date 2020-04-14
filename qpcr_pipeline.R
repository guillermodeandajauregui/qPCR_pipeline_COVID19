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

resultados <- "data/Procolo_COVID-19_Prueba1_4Abr20.eds"

################################################################################
#analysis steps
################################################################################

#####1 extract delta RN curves #################################################

my_deltaRN <- tidy_deltaRN(resultados)

##########Plot them for that nice visual inspection feel #######################

pivot_deltaRN(my_deltaRN) %>% 
  plot_deltaRN.long()

ggsave(filename = "results/testCurves_20200404.pdf", 
       units = "mm", 
       width = 210, 
       height = 148)
  
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

my_deltaRN.long %>% 
  group_by(sample.label) %>% 
  filter(sample.label == "Control COVID 1", probe == "N1")

#### this part to extract the unique values, to use them for iterations

my_wells <- unique(my_deltaRN.long$well)
names(my_wells) <- my_wells

all_probes <- unique(my_deltaRN.long$probe)
names(all_probes) <- all_probes



#we are evaluating same sample label in many wells 
#so we have to filter that as well
#check with analytics if this will be the expected behavior 

my_deltaRN.long <- 
my_deltaRN.long %>% 
  mutate(sample.label = paste0("well_", well, " ", sample.label))

my_samples <- unique(my_deltaRN.long$sample.label)
names(my_samples) <- my_samples


#will keep only the first well for each probe in Controls 
my_deltaRN.long.fk <- 
my_deltaRN.long %>%
  group_by(sample.label, probe) %>% 
  filter(well == min(as.numeric(well))) %>% 
  arrange(well, cycles)

#####5 Evaluate curves (NOT READY YET) #########################################

### here we iterate over each well
### returning whether the probe passed the threshold or not in that well


my_results <- 
#lapply(X = my_wells, FUN = function(my_well){
  lapply(X = my_samples, FUN = function(my_sample){  
  #get only curves for that well
  #well_data <-
  
  sample_data <-
    my_deltaRN.long.fk %>% 
    filter(sample.label == my_sample) #%>% 
    #filter(well == min(as.numeric(well)))
  

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
        #make a plot to show 
        p <- 
        plot_deltaRN.long(the_curve) +
          geom_hline(yintercept = the_threshold) 
        plot(p)
          
        #does the curve crosses the threshold?
        #threshold for RP should be crossed at time 35
        any(the_curve$value[1:40] > the_threshold)
      },
      error =function(cond){
        message("well does not have that probe")
        return(NA)
      }
    )
  }) %>% bind_rows() #we make them a row
}) %>% bind_rows(.id = "well") #and join them in a data frame


##### plotting section


the_plots <- 
lapply(X = my_samples, FUN = function(my_sample){
  
  #get only curves for that well
  sample_data <-
    my_deltaRN.long.fk %>% 
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
        #make a plot to show 
        p <- 
          plot_deltaRN.long(the_curve) +
          geom_hline(yintercept = the_threshold) 
        
      },
      error =function(cond){
        message("well does not have that probe")
        return(NA)
      }
    )
  }) 
  
})




#####reannotate results  #######################################################

##### we will have an id table for the samples; let's mock one up

annotationFile <-
my_deltaRN.long %>%
  select(well) %>%
  unique() %>%
  mutate(well = as.numeric(well)) %>% 
  mutate(well.id = paste0(LETTERS[(well%/%12 + 1)], (well%%12 + 1))) %>% 
  mutate(well = as.character(well)) %>% 
  mutate(id.external = paste0("sample_", nrow(my_deltaRN.long)),
         id.internal = paste0("inmegen_", 
                              str_pad(1:nrow(my_deltaRN.long), 
                                      2, 
                                      side = "left", 
                                      pad = 0)
                              )
         )

my_results <- 
  left_join(annotationFile, my_results)

##### Finally, we evaluate the conditions based on the rules set in the PNO

#### THIS IS A MOCKUP, WHERE HAVING S == TRUE means POSITIVE for COVID
#### change this logic tests for the actual PNO conditions
my_results <-
my_results %>% 
  purrr::map_df(.f = function(i){
    tidyr::replace_na(i, replace = "undetermined")}
    ) %>% 
    mutate(diagnosis = ifelse(N1=="TRUE", "positive", 
                            ifelse(N2=="TRUE", "positive", "negative")
                            )
         )  


################################################################################
#write out - reports
################################################################################

#we rename the the_plots object so that we can iterate over it easily 


#names(the_plots) <- paste0("well_", names(the_plots))
library("rmarkdown")

lapply(seq_along(the_plots), FUN = function(i){
  
  #the_well_is <- str_split(names(the_plots)[i], pattern = "_", simplify = T)[2] %>% 
  #  as.numeric()
  
  the_sample_is <- names(the_plots)[i]
  
  my_r <- 
  my_results %>% 
    filter(well == the_sample_is) 

    my_name <- names(the_plots)[i]
    mea_plote <- the_plots[[i]]
    outpath <- paste0("results/reports/", Sys.Date(), "_", my_name, ".pdf")
    render("template.Rmd",output_file = outpath)    
})

################################################################################
#write out 
################################################################################

out_path = "results/prueba.txt"
write_delim(my_results, out_path)
