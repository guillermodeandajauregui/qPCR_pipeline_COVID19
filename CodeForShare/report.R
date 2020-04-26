################################################################################
#
# R code for generat reports - COVID19 eds app/COVID binnacle
# by: INMEGEN Computational Genomics Dept
# Hugo Tovar hatovar@inmegen.gob.mx and
# Guillermo de Anda J??uregui gdeanda@inmegen.edu.mx
#
################################################################################




################################################################################
#required libs
################################################################################
library("rmarkdown")


################################################################################
#define analysis functions here
################################################################################

report <- read.csv("test_diagnosis.csv")




## simulate data

## "batch" object contains the information that will be provided by the eds analysis system.
batch <- rbind(report, report, report, report, report, report, report, report, report, report, report, report)
batch$qc <- c(rep("PASS", 12), rep("FAIL", 12), rep("PASS", 12), rep("FAIL", 12))
batch$ptc.pass <- c(rep("TRUE", 12), rep("FALSE", 12), rep("TRUE", 12), rep("FALSE", 12))
batch$sample <- paste("INM", sub('(^[0-9]$)','0\\1', 1:48), sep = "")
batch <- cbind(plate = c(rep("plate1", 12), rep("plate2", 12), rep("plate3", 12), rep("plate4", 12)), batch)


## "metadata" object contains the information that the eds analysis system does not have and will have to be added by the COVID binnacle
metadata <- data.frame(client = c(rep("Instituto Uno", 24), rep("Instituto Dos", 24)),
                      id_client = c(rep("inst01", 24), rep("inst02", 24)), 
                      id_samples_client = c(paste("isnt1_", sub('(^[0-9]$)','0\\1', 1:24), sep = ""), paste("isnt2_", sub('(^[0-9]$)','0\\1', 1:24), sep = "")), 
                      patient_name = c("Louis Armstrong", "Duke Ellington", "Miles Davis", "Charlie Parker", "John Coltrane", "Dizzy Gillespie", "Billie Holiday", "Thelonious Monk", "Charles Mingus", "Count Basie", "Lester Young", "Coleman Hawkins", "Ella Fitzgerald", "Sonny Rollins", "Sidney Bechet", "Art Blakey", "Ornette Coleman", "Benny Goodman", "Jelly Roll Morton", "Bill Evans", "Art Tatum", "Clifford Brown", "Stan Getz", "Sarah Vaughan", "Bud Powell", "Fletcher Henderson", "Django Reinhardt", "Herbie Hancock", "Wayne Shorter", "Horace Silver", "Dave Brubeck", "Rahsaan Roland Kirk", "Cecil Taylor", "King Oliver", "Sun Ra", "Gil Evans", "Lionel Hampton", "Art Pepper", "Eric Dolphy", "Oscar Peterson", "Charlie Christian", "Ben Webster", "Fats Waller", "Earl Hines", "Woody Herman", "Wes Montgomery", "J.J. Johnson", "John McLaughlin"))



result_table <- cbind(metadata, batch)


######make_reports

# Preguntar como serán selecciona los reportes tanto por muestra como por institución 
# ¿por intervalo de fecha? ¿se generarán diariamente para todas las instituciones trabajadas?
# En esta versión se trabaja en un reporte para una institución definida por el usuario.
#

make_reports <- function(result_table, # data with columns like simulated data
                         #input, # The origin of the data will be defined by the DDT
                         id_client, # Internal ID of the requesting institution
                         outdir, # The report output directory will be defined either by the user or by the DDT
                         by_sample = FALSE){
    # Select conclusive data for specific petitioner institution
    the_samples <- which(result_table$id_client %in% id_client 
    & result_table$qc == "PASS" 
    & (result_table$classification == "positive" 
      | result_table$classification == "negative"))
  # make reports by petitioner institution with those samples that both your plate has passed 
  # the controls and its diagnosis is conclusive
  if(by_sample == FALSE){
    # select specific columns 
    my_r <- as.matrix(result_table[the_samples, c(3, 4, 11:14)])
    my_r[my_r == "Inf"] <- "45+"
    client <- unique(result_table$client[the_samples])
    n_samples <- length(rownames(my_r))
    rownames(my_r) <- sub('(^[0-9]$)','0\\1', 1:n_samples)
    outpath <- paste0(outdir, "/", Sys.Date(), "_", the_client_is, ".pdf") 
    render("report.Rmd",output_file = outpath)
  # make reports by patient with those samples that both your plate has passed the controls and its
  # diagnosis is conclusive
  }else{
    lapply(seq_along(result_table$sample[the_samples]), function(i){
      client <- result_table$client[the_samples][i]
      idsamplesclient <- result_table$id_samples_client[the_samples][i]
      patient_name <- result_table$patient_name[the_samples][i]
      diagnosis <- ifelse(result_table$classification[the_samples][i] == "positive", "SÍ", "NO")
      r_results <- data.frame(probe = c("N1", "N2", "RP"), Ct = c(ifelse(result_table$N1[the_samples][i] == Inf, "45+", result_table$N1[the_samples][i]), ifelse(result_table$N2[the_samples][i] == Inf, "45+", result_table$N2[the_samples][i]), ifelse(result_table$RP[the_samples][i] == Inf, "45+", result_table$RP[the_samples][i])))
      outpath <- paste0(outdir, "/", Sys.Date(), "_", idsamplesclient, ".pdf")
      render("sample.Rmd",output_file = outpath)
    })
  }
}