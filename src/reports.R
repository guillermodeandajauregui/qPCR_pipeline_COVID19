################################################################################
#
# R code for generat reports - COVID19 eds app/COVID binnacle
# by: INMEGEN Computational Genomics Dept
# Hugo Tovar hatovar@inmegen.gob.mx and
# 
################################################################################




################################################################################
#required libs
################################################################################
library("rmarkdown")


######################################################################################
######################################################################################
#                                    s     r    c                                    #
######################################################################################
######################################################################################

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
    my_r <- as.matrix(result_table)[,c("sample", "N1", "N2", "RP")]
    my_r[grep("Inf", my_r)] <- "45+"
    ntc <- grep(pattern = "NTC", x = names(plot_list))
  ptc <- grep(pattern = "PTC", x = names(plot_list))
  exc <- grep(pattern = "EC", x = names(plot_list))
    outpath <- paste0(outdir, "/", Sys.Date(), "_", plate, ".pdf")
    render("template_qc.Rmd",output_file = outpath)
  }
}

## Extracts the plaque processing date from the sample INMEGEN ID
getDecade <- function(idinmegen){
  sapply(sapply(idinmegen, 
            function(x){strsplit(x, split = "-")}), 
    function(x){x[2]})
}

getMonth <- function(idinmegen){
  tail <- sapply(sapply(idinmegen, 
            function(x){strsplit(x, split = "-")}), 
    function(x){x[4]})
  monthID <- substr(tail, start = 3, stop = 3)
  month <- sapply(idinmegen, function(x){which(c("E", "F", "Z", "A", "Y", "J", "L", "G", "S", "O", "N", "D") %in% monthID)})
  return(month)
}

getDay <- function(idinmegen){
  tail <- sapply(sapply(idinmegen, 
            function(x){strsplit(x, split = "-")}), 
    function(x){x[4]})
  day <- substr(tail, start = 1, stop = 2)
  return(day)
}

getDate <- function(idinmegen){
  myDate <- as.Date(paste(paste("20", getDecade(idinmegen), sep = ""), 
                        getMonth(idinmegen), 
                        getDay(idinmegen), sep ="-"))
  return(myDate)
}

## Extract the client ID from the sample INMEGEN ID
getClientID <- function(idinmegen){
  sapply(sapply(idinmegen, 
            function(x){strsplit(x, split = "-")}), 
    function(x){x[3]})
}

# ######################################################################################
# ######################################################################################
# #                                   Inform Lab                                       #
# ######################################################################################
# ######################################################################################


# Print the plate report for the lab
plateBooklet <- function(results, outdir, qc_results){
    qcdf <- data.frame(controType = c("Positivo", "Negativo", "Extracción"), 
                    externalName = c("nCoVPC", "NTC", "HSC"), 
                    used = c("Falla sustancial del reactivo, incluida la integridad del primer y la sonda.", "Contaminación de reactivos y/o ambiente.", "Falla en el procedimiento de lisis y extracción, contaminación potencial durante la extracción."))
    qc_table <- data.frame(qc_results$qc.values[c(2,1,3),c(3,4,2)])
    my_qc <- cbind(qcdf,qc_table)
    rownames(my_qc) <- c("PTC", "NTC", "EC")
    my_qc["PTC", c("N1", "N2", "RP")] = ifelse(qc_results$ptc.pass, paste("\\color{teal}{", my_qc["PTC", c("N1", "N2", "RP")], "}", sep = "") , paste("\\cellcolor{red!50}{", my_qc["PTC", c("N1", "N2", "RP")], "}", sep = "") )
    my_qc["NTC", c("N1", "N2", "RP")] = ifelse(qc_results$ntc.pass, paste("\\color{teal}{", "+45", "}", sep = ""), paste("\\cellcolor{red!50}{", my_qc["NTC", c("N1", "N2", "RP")], "}", sep = "") )
    my_qc["EC", c("N1", "N2", "RP")] = ifelse(is.na(my_qc$N1[3]), my_qc["EC", c("N1", "N2", "RP")], ifelse(qc_results$ec.pass, paste("\\color{teal}{", "+45", "}", sep = ""), paste("\\cellcolor{red!50}{", my_qc["EC", c("N1", "N2", "RP")], "}", sep = "")))

    my_r <- results[,c("plate", "sample", "N1", "N2", "RP", "classification")]
      my_r <- cbind(nsample = sub('(^[0-9]$)','0\\1', 1:length(my_r[,1])), my_r)
      my_r[,"N1"] = ifelse(my_r[,"N1"] == "Inf", "Indeterminado", paste("\\cellcolor{red!50}{", my_r[,"N1"], "}", sep = ""))
      my_r[,"N2"] = ifelse(my_r[,"N2"] == "Inf", "Indeterminado", paste("\\cellcolor{red!50}{", my_r[,"N2"], "}", sep = ""))
      my_r[,"classification"] = ifelse(my_r[,"classification"] == "positive", "\\cellcolor{red!50}{positivo}", 
                                      ifelse(my_r[,"classification"] == "negative", "\\cellcolor{yellow!50}{negativo}", 
                                            ifelse(my_r[,"classification"] == "inconclusive", "\\cellcolor{cyan!50}{inconclusive}", 
                                                  ifelse(my_r[,"classification"] == "invalid", "\\cellcolor{cyan!50}{invalid}", 
                                                        ifelse(my_r[,"classification"] == "edge_positive", "\\cellcolor{red!25}{positivo marginal}", 
                                                          my_r[,"classification"])))))
    processingdate <- unique(getDate(my_r$sample))
    plate <- unique(results$plate)
    outpath <- paste0(outdir, "/", Sys.Date(), "_", plate, "_", "plateBooklet.pdf") 
    render("plateBooklet.Rmd",output_file = outpath)
}




# ######################################################################################
# ######################################################################################
# #                                 Inform Client                                      #
# ######################################################################################
# ######################################################################################

## Print the results report for all clients
resultsReport <- function(results, outdir){
    the_samples <- which(results$qc == "PASS" 
        & (results$classification == "positive" 
          | results$classification == "negative"))
    my_r <- results[the_samples, c("sample", "classification")]
            my_r <- cbind(nsample = sub('(^[0-9]$)','0\\1', 1:length(my_r[,1])), my_r)
            my_r$classification = ifelse(my_r$classification == "positive", paste("\\cellcolor{red!50}{", "positive", "}", sep = ""),"negative")
    client <- getClientID(my_r$sample)
    n_samples <- length(rownames(my_r))
    allsamples <- length(rownames(results))
    processingdate <- unique(getDate(my_r$sample))
    client_file <- paste0(gsub("[^[:alnum:]]", "_", unique(client[order(client)])), collapse = "-")# Claimant name institution safe for file name
    outpath <- paste0(outdir, "/", Sys.Date(), "_", client_file, "_", "report.pdf") 
    render("resultsReport.Rmd",output_file = outpath)
  }

## Print the results report by clients

makeReports <- function(table_diagnosis, # data with columns like simulated data
                         outdir # The report output directory will be defined either by the user or by the DDT
                         ){
    results <- cbind(table_diagnosis, id_client = getClientID(table_diagnosis$sample))
    lapply(split(results,results$id_client), resultsReport, outdir = outdir)
  }



