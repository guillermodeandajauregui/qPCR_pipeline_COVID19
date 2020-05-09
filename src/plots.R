################################################################################
#
# Plot and report functions for qPCR analysis - COVID19 detection
#by: INMEGEN Computational Genomics Dept
#Guillermo de Anda J??uregui gdeanda@inmegen.edu.mx
# and Hugo Tovar hatovar@inmegen.gob.mx
#
################################################################################


################################################################################
#required libs
################################################################################
library("cowplot")
library("ggpubr")

################################################################################
#define analysis functions here
################################################################################

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
    theme_minimal()+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(1, 1, 1, 1, "mm"),
          plot.background = element_rect(
    							fill = NULL,
    							colour = "darkgrey",
    							size = 0.3)
          )
}


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
          geom_hline(yintercept = the_threshold, linetype = 2, colour = "maroon") 
      })
      
      
    })
}


triplets <- function(curve.list){
	triplet <- lapply(seq_along(curve.list), function(i){
	plot_grid(curve.list[[i]]$N1, curve.list[[i]]$N2, curve.list[[i]]$RP, 
		labels = c("N1", "N2", "RP"),
		ncol = 1, nrow = 3,
		label_size = 11,
		hjust = 1)
	})
	triplets <- lapply(seq_along(triplet), function(i){
	annotate_figure(triplet[[i]], bottom = text_grob("Cycles", size = 10, hjust = 1),
                left = text_grob("Delta RN", size = 10, rot = 90))
	})
	names(triplets) <- names(curve.list)
	return(triplets)
}

######make_reports

make_reports <- function(plot_list, 
                         result_table,
                         input,
                         outdir, 
                         qc_results,
                         qc = F){
  plate <- stringr::str_remove(string = basename(input), pattern = ".eds")
  qcplate <- ifelse(qc_results == "PASS", "true", "")
  #makes reports from a list of plots and some result table
  if(qc==F){
  lapply(seq_along(plot_list), function(i){
    
    the_sample_is <- names(plot_list)[i]
    
    my_r <- 
      result_table %>% 
      filter(sample == the_sample_is)
    
    my_name <- names(plot_list)[i]
    mea_plote <- plot_list[i]
    
    outpath <- paste0(outdir, "/", Sys.Date(), "_", my_name, ".pdf")
    outpath_inf <- paste0(outdir, "/", Sys.Date(), "_", my_name, "_results.pdf")
    # render("template_inf.Rmd", output_file = outpath_inf)
    render("template.Rmd",output_file = outpath)})
  }else{
  	my_r <- as.matrix(result_table)
  	my_r[my_r == "Inf"] <- "45+"
    ntc <- grep(pattern = "NTC", x = names(plot_list))
	ptc <- grep(pattern = "PTC", x = names(plot_list))
 	exc <- grep(pattern = "EC", x = names(plot_list))
    outpath <- paste0(outdir, "/", Sys.Date(), "_", plate, ".pdf")
    render("template_qc.Rmd",output_file = outpath)
  }
}