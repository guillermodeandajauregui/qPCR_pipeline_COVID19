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

#report <- read.csv("test_diagnosis.csv")
report1 <- read.table("reportstest.results.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
report2 <- read.table("reportstest.results2.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)


## simulate data

batch <- rbind(report1, report2)

## "batch" object contains the information that will be provided by the eds analysis system.
# batch <- rbind(report, report, report, report, report, report, report, report, report, report, report, report)
# batch$qc <- c(rep("PASS", 12), rep("FAIL", 12), rep("PASS", 12), rep("FAIL", 12))
# batch$ptc.pass <- c(rep("TRUE", 12), rep("FALSE", 12), rep("TRUE", 12), rep("FALSE", 12))
# batch$sample <- paste("INM", sub('(^[0-9]$)','0\\1', 1:48), sep = "")
# batch <- cbind(plate = c(rep("plate1", 12), rep("plate2", 12), rep("plate3", 12), rep("plate4", 12)), batch)

names <- c("Louis Armstrong", "Duke Ellington", "Miles Davis", "Charlie Parker", "John Coltrane", "Dizzy Gillespie", "Billie Holiday", "Thelonious Monk", "Charles Mingus", "Count Basie", "Lester Young", "Coleman Hawkins", "Ella Fitzgerald", "Sonny Rollins", "Sidney Bechet", "Art Blakey", "Ornette Coleman", "Benny Goodman", "Jelly Roll Morton", "Bill Evans", "Art Tatum", "Clifford Brown", "Stan Getz", "Sarah Vaughan", "Bud Powell", "Fletcher Henderson", "Django Reinhardt", "Herbie Hancock", "Wayne Shorter", "Horace Silver", "Dave Brubeck", "Rahsaan Roland Kirk", "Cecil Taylor", "King Oliver", "Sun Ra", "Gil Evans", "Lionel Hampton", "Art Pepper", "Eric Dolphy", "Oscar Peterson", "Charlie Christian", "Ben Webster", "Fats Waller", "Earl Hines", "Woody Herman", "Wes Montgomery", "J.J. Johnson", "John McLaughlin")

## "metadata" object contains the information that the eds analysis system does not have and will have to be added by the COVID binnacle
metadata <- data.frame(client = c(rep("Instituto Uno", 3), rep("Instituto Dos", 4), rep("Instituto Tres", 18), rep("Instituto Cuatro", 21)),
                      id_client = c(rep("A", 3), rep("D", 4), rep("E", 18), rep("F", 21)), 
                      id_samples_client = c(paste("isnt1_", sub('(^[0-9]$)','0\\1', 1:3), sep = ""), paste("isnt2_", sub('(^[0-9]$)','0\\1', 1:4), sep = ""), paste("isnt3_", sub('(^[0-9]$)','0\\1', 1:18), sep = ""), paste("isnt4_", sub('(^[0-9]$)','0\\1', 1:21), sep = "")), 
                      patient_name = names[-c(47,48)],
                      id_inmegen = c("CV-20-A-06Y000022", "CV-20-A-06Y000023", "CV-20-A-06Y000024", "CV-20-D-06Y000088", "CV-20-D-06Y000089", "CV-20-D-06Y000090", "CV-20-D-06Y000091", "CV-20-E-06Y000073", "CV-20-E-06Y000074", "CV-20-E-06Y000075", "CV-20-E-06Y000076", "CV-20-E-06Y000077", "CV-20-E-06Y000078", "CV-20-E-06Y000079", "CV-20-E-06Y000080", "CV-20-E-06Y000082", "CV-20-E-06Y000083", "CV-20-E-06Y000084", "CV-20-E-06Y000085", "CV-20-E-06Y000086", "CV-20-E-06Y000087", "CV-20-E-06Y000088", "CV-20-E-06Y000089", "CV-20-E-06Y000096", "Rep CV-20-E-05Y000041", "CV-20-F-06Y000001", "CV-20-F-06Y000002", "CV-20-F-06Y000003", "CV-20-F-06Y000004", "CV-20-F-06Y000005", "CV-20-F-06Y000006", "CV-20-F-06Y000007", "CV-20-F-06Y000008", "CV-20-F-06Y000009", "CV-20-F-06Y000010", "CV-20-F-06Y000011", "CV-20-F-06Y000012", "CV-20-F-06Y000013", "CV-20-F-06Y000014", "CV-20-F-06Y000015", "CV-20-F-06Y000016", "CV-20-F-06Y000017", "CV-20-F-06Y000018", "CV-20-F-06Y000019", "CV-20-F-06Y000020", "CV-20-F-06Y000021"),
                      method = sample(c("hisopado nasofaríngeo", "saliva", "lavado bronquial"), size = 46, replace = TRUE))

which(batch$sample %in% metadata$id_inmegen)
  
result_table <- cbind(metadata, batch[match(metadata$id_inmegen, batch$sample),])


######make_reports

# Preguntar como serán selecciona los reportes tanto por muestra como por institución 
# ¿por intervalo de fecha? ¿se generarán diariamente para todas las instituciones trabajadas?
# En esta versión se trabaja en un reporte para una institución definida por el usuario.
#

## Extrae la fecha de procesamiento de placa a partir del ID de muestra del INMEGEN
id2date <- function(idinmegen){
codeMonth <- c("E", "F", "Z", "A", "Y", "J", "L", "G", "S", "O", "N", "D")
myDate <- as.Date(paste(paste("20", substr(idinmegen, start = 4, stop = 5), sep = ""), which(codeMonth %in% substr(idinmegen, start = 11, stop = 11)), substr(idinmegen, start = 9, stop = 10), sep ="-"))
return(myDate)
}

reportSamples <- function(result_table, outdir){
    Sys.setenv(LANG = "es")
    the_samples <- which(result_table$qc == "PASS" 
        & (result_table$classification == "positive" 
          | result_table$classification == "negative"))
    my_r <- result_table[the_samples, c("id_inmegen", "id_samples_client", "id_client", "patient_name", "method", "classification")]
            my_r <- cbind(nsample = sub('(^[0-9]$)','0\\1', 1:length(my_r[,1])), my_r)
            my_r$classification = ifelse(my_r$classification == "positive", paste("\\cellcolor{red!50}{", "positive", "}", sep = ""),"negative")
            my_r$id_samples_client <- gsub("_", "-", my_r$id_samples_client)
    client <- strsplit(my_r$id_inmegen[1], split = "-")[[1]][3]
    client_file <- gsub("[^[:alnum:]]", "_", client) # Name of the claimant institution safe for file name
    n_samples <- length(rownames(my_r))
    allsamples <- length(rownames(result_table))
    processingdate <- id2date(my_r$id_inmegen[1])
    outpath <- paste0(outdir, "/", Sys.Date(), "_", client_file, "_", "report.pdf") 
    render("report.Rmd",output_file = outpath)
  }

make_reports <- function(result_table, # data with columns like simulated data
                         #input, # The origin of the data will be defined by the DDT
                         #id_client, # Internal ID of the requesting institution
                         outdir, # The report output directory will be defined either by the user or by the DDT
                         report = TRUE){
      if(report == TRUE){
        lapply(split(result_table,result_table$id_client), reportSamples, result_table = result_table, outdir = outdir)
      }else{
  # make reports by patient with those samples that both your plate has passed the controls and its
  # diagnosis is conclusive
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




df <-structure(list(Server =structure(1:2, .Label =c("Server1","Server2"), class = "factor"), CPU =c(79.17, 93), UsedMemPercent =c(16.66,18.95)), .Names =c("Server", "CPU", "UsedMemPercent"), row.names =c(NA,-2L), class = "data.frame")

df[, 2] =ifelse(df[, 2]>80,paste("\\color{red}{",round(df[, 2], 2), "}"),round(df[, 2], 2))



\cellcolor{red!50}
# What you need

kable(df, "latex", escape = F)


dicMonth <- data.frame(codeMonth = c("E", "F", "Z", "A", "Y", "J", "L", "G", "S", "O", "N", "D"), month =c("enero", "febrero", "marzo", "abril", "mayo", "junio", "julio", "agosto", "septiembre", "octubre", "noviembre", "diciembre"), numMonth = 1:12)

