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



######################################################################################
######################################################################################
#                                 M E T A D A T A                                    #
######################################################################################
######################################################################################

## "metadata" object contains the information that the eds analysis system does not have and will have to be added by the COVID binnacle

names <- c("Louis Armstrong", "Duke Ellington", "Miles Davis", "Charlie Parker", "John Coltrane", "Dizzy Gillespie", "Billie Holiday", "Thelonious Monk", "Charles Mingus", "Count Basie", "Lester Young", "Coleman Hawkins", "Ella Fitzgerald", "Sonny Rollins", "Sidney Bechet", "Art Blakey", "Ornette Coleman", "Benny Goodman", "Jelly Roll Morton", "Bill Evans", "Art Tatum", "Clifford Brown", "Stan Getz", "Sarah Vaughan", "Bud Powell", "Fletcher Henderson", "Django Reinhardt", "Herbie Hancock", "Wayne Shorter", "Horace Silver", "Dave Brubeck", "Rahsaan Roland Kirk", "Cecil Taylor", "King Oliver", "Sun Ra", "Gil Evans", "Lionel Hampton", "Art Pepper", "Eric Dolphy", "Oscar Peterson", "Charlie Christian", "Ben Webster", "Fats Waller", "Earl Hines", "Woody Herman", "Wes Montgomery", "J.J. Johnson", "John McLaughlin", "Artie Shaw", "Lee Morgan", "David Murray", "Chick Corea", "Max Roach", "Roy Eldridge", "Modern Jazz Quartet", "Anthony Braxton", "Bix Beiderbecke", "Dexter Gordon", "Keith Jarrett", "Lee Konitz", "Stan Kenton", "Joe Henderson", "Gerry Mulligan", "Benny Carter", "Teddy Wilson", "Freddie Hubbard", "Cannonball Adderley", "McCoy Tyner", "Chet Baker", "Lennie Tristano", "Jimmy Smith", "Mary Lou Williams", "George Russell", "Fats Navarro", "Bennie Moten", "Jimmie Lunceford", "Wynton Marsalis", "Albert Ayler", "Charlie Haden", "Erroll Garner", "Meade Lux Lewis", "Pat Metheny", "Jack Teagarden", "Johnny Hodges", "Chick Webb", "James P. Johnson", "Jimmy Giuffre", "Jaco Pastorius", "Hank Mobley", "Elvin Jones", "Evan Parker", "Paul Chambers", "Ron Carter", "Carla Bley", "Bennie Golson", "Jackie McLean", "James Carter", "Donald Byrd", "Johnny Dodds", "Glenn Mille")
id_inmegen <- read.table(file = "data/idINMEGEN.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]
id_client <- substr(id_inmegen, 7,7)
Replace <- data.frame(from = c("F", "A", "E", "D"), to = c("Instituto 4", "Instituto 1", "Instituto 2", "Instituto 3"))
client <- DataCombine::FindReplace(data.frame(a=substr(id_inmegen, 7,7)), Var = "a",  replaceData = Replace )[[1]]
id_samples_client <- id_client
      id_samples_client[which(id_client == "F")] <- paste("isnt4_", sub('(^[0-9]$)','0\\1', 1:length(id_samples_client[which(id_client == "F")])), sep = "")
      id_samples_client[which(id_client == "A")] <- paste("isnt1_", sub('(^[0-9]$)','0\\1', 1:length(id_samples_client[which(id_client == "A")])), sep = "")
      id_samples_client[which(id_client == "E")] <- paste("isnt2_", sub('(^[0-9]$)','0\\1', 1:length(id_samples_client[which(id_client == "E")])), sep = "")
      id_samples_client[which(id_client == "D")] <- paste("isnt3_", sub('(^[0-9]$)','0\\1', 1:length(id_samples_client[which(id_client == "D")])), sep = "")
method = sample(c("hisopado nasofaríngeo", "saliva", "lavado bronquial"), size = length(id_inmegen), replace = TRUE)
notes <- sample(c("", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "", "", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", "", "", "", "Repetir extración y RT-PCR. Solo uno de los dos blancos para CoV-2 amplificó.", "", "", ""), size = length(id_inmegen), replace = TRUE)
plate <- rep("a", 85)
plate[1:24] <- "P-01-06Y-2CQP-1" 
plate[25:48] <- "P-01-09Y-2CQP-1v2" 
plate[49:64] <- "P-01-30A-2CQP-1" 
plate[65:85] <- "P-04-06Y-2CQP-1"

metadata <- data.frame(plate = plate,
                      id_inmegen = id_inmegen,
                      client = client,
                      id_client = id_client, 
                      id_samples_client = id_samples_client, 
                      patient_name = names[1:length(id_inmegen)],
                      method = method,
                      notes= notes)


input <- "data/P-01-06Y-2CQP-1.res"
qc_input <- "data/P-01-06Y-2CQP-1.qc"

results <- read.table(input, header = TRUE, stringsAsFactors = FALSE)

# result_table contains informatión generated by ARPA plus "metadata" by plate
result_table <- cbind(metadata[metadata$id_inmegen %in% results$sample,], results[match(metadata[metadata$id_inmegen %in% results$sample,]$id_inmegen, results$sample),])


######################################################################################
######################################################################################
#                                    s     r    c                                    #
######################################################################################
######################################################################################


## Extrae la fecha de procesamiento de placa a partir del ID INMEGEN de muestra
getDate <- function(idinmegen){
  codeMonth <- c("E", "F", "Z", "A", "Y", "J", "L", "G", "S", "O", "N", "D")
  myDate <- as.Date(paste(paste("20", substr(idinmegen, start = 4, stop = 5), sep = ""), 
                        sapply(idinmegen, function(x){which(codeMonth %in% substr(x, start = 11, stop = 11))}), 
                        substr(idinmegen, start = 9, stop = 10), sep ="-"))
  return(myDate)
}

## Extrae el ID del cliente del ID INMEGEN de la muestra
getClientID <- function(idinmegen){
  sapply(sapply(idinmegen, 
            function(x){strsplit(x, split = "-")}), 
    function(x){x[3]})
}

## Imprime el informe de resultados para el cliente
resultsReport <- function(result_table, outdir){
    the_samples <- which(result_table$qc == "PASS" 
        & (result_table$classification == "positive" 
          | result_table$classification == "negative"))
    my_r <- result_table[the_samples, c("id_inmegen", "id_samples_client", "client", "patient_name", "method", "classification")]
            my_r <- cbind(nsample = sub('(^[0-9]$)','0\\1', 1:length(my_r[,1])), my_r)
            my_r$classification = ifelse(my_r$classification == "positive", paste("\\cellcolor{red!50}{", "positive", "}", sep = ""),"negative")
            my_r$id_samples_client <- gsub("[^[:alnum:]]", "-", my_r$id_samples_client) # Cleint ID sample safe for kable/latex
    client <- unique(result_table$client)
    n_samples <- length(rownames(my_r))
    allsamples <- length(rownames(result_table))
    processingdate <- unique(getDate(my_r$id_inmegen))
    client_file <- paste0(gsub("[^[:alnum:]]", "_", client[order(client)]), collapse = "-")# Claimant name institution safe for file name
    outpath <- paste0(outdir, "/", Sys.Date(), "_", client_file, "_", "report.pdf") 
    render("resultsReport.Rmd",output_file = outpath)
  }

# Imprime el informe de placa para el lab
plateBooklet <- function(result_table, outdir, qc_input){
    #plate <- stringr::str_remove(string = basename(input), pattern = ".res")
    plate <- unique(result_table$plate)

    qcdf <- data.frame(controType = c("Positivo", "Negativo", "Extracción"), 
                    externalName = c("nCoVPC", "NTC", "HSC"), 
                    used = c("Falla sustancial del reactivo, incluida la integridad del primer y la sonda.", "Contaminación de reactivos y/o ambiente.", "Falla en el procedimiento de lisis y extracción, contaminación potencial durante la extracción."))
    qc_results <- read.table(qc_input, header = TRUE, sep = " ", row.names = 1)[c("PTC", "NTC", "EC"),c("N1","N2","RP")]
    my_qc <- cbind(qcdf,qc_results)
    rownames(my_qc) <- c("PTC", "NTC", "EC")
    my_qc["PTC", c("N1", "N2", "RP")] = ifelse(my_qc["PTC", "N1"] <=38 & my_qc["PTC", "N2"] <=38 & my_qc["PTC", "RP"] >=35, paste("\\color{green}{", my_qc["PTC", c("N1", "N2", "RP")], "}", sep = "") , paste("\\cellcolor{red!50}{", my_qc["PTC", c("N1", "N2", "RP")], "}", sep = "") )
    my_qc["NTC", c("N1", "N2", "RP")] = ifelse(my_qc["NTC", c("N1", "N2", "RP")] == Inf, paste("\\color{green}{", "+45", "}", sep = ""), paste("\\cellcolor{red!50}{", my_qc["NTC", c("N1", "N2", "RP")], "}", sep = "") )
    my_qc["EC", c("N1", "N2", "RP")] = ifelse(is.na(my_qc$N1[3]), my_qc["EC", c("N1", "N2", "RP")], ifelse(my_qc["EC", c("N1", "N2", "RP")] == Inf, paste("\\color{green}{", "+45", "}", sep = ""), paste("\\cellcolor{red!50}{", my_qc["EC", c("N1", "N2", "RP")], "}", sep = "")))

    my_r <- result_table[,c("plate", "id_samples_client", "id_inmegen", "N1", "N2", "RP", "classification", "notes")]
      my_r <- cbind(nsample = sub('(^[0-9]$)','0\\1', 1:length(my_r[,1])), my_r)
      my_r$id_samples_client <- gsub("[^[:alnum:]]", "-", my_r$id_samples_client) # Cleint ID sample safe for kable/latex
      my_r[,"N1"] = ifelse(my_r[,"N1"] == "Inf", "Indeterminado", paste("\\cellcolor{red!50}{", my_r[,"N1"], "}", sep = ""))
      my_r[,"N2"] = ifelse(my_r[,"N2"] == "Inf", "Indeterminado", paste("\\cellcolor{red!50}{", my_r[,"N2"], "}", sep = ""))
      my_r[,"classification"] = ifelse(my_r[,"classification"] == "positive", "\\cellcolor{red!50}{positivo}", 
                                      ifelse(my_r[,"classification"] == "negative", "\\cellcolor{yellow!50}{negativo}", 
                                            ifelse(my_r[,"classification"] == "inconclusive", "\\cellcolor{cyan!50}{inconclusive}", 
                                                  ifelse(my_r[,"classification"] == "invalid", "\\cellcolor{cyan!50}{invalid}", 
                                                        ifelse(my_r[,"classification"] == "edge_positive", "\\cellcolor{red!25}{positivo marginal}", 
                                                          my_r[,"classification"])))))
    processingdate <- unique(getDate(my_r$id_inmegen))
    outpath <- paste0(outdir, "/", Sys.Date(), "_", plate, "_", "plateBooklet.pdf") 
    render("plateBooklet.Rmd",output_file = outpath)
}




######################################################################################
######################################################################################
#                                    a     p    p                                    #
######################################################################################
######################################################################################



make_reports <- function(result_table, # data with columns like simulated data
                         outdir, # The report output directory will be defined either by the user or by the DDT
                         report = FALSE){
      if(report == FALSE){
        lapply(split(result_table,result_table$id_client), resultsReport, outdir = outdir)
      }else{
  # make reports by patient with those samples that both your plate has passed the controls and its
  # diagnosis is conclusive
    lapply(split(result_table,result_table$plate), plateBooklet, outdir = outdir)
  }
}


