#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinyFiles)
library(DT)
library(shinyauthr)
library(shinyjs)
library(sodium)



source("function_app_qPCR_cdc.R")
source("function_plate_conf.R")
source("src/functions.R")
source("src/reports.R")

#input <- c("data/Procolo_COVID-19_Prueba1_4Abr20.eds")
#output <- c("results/")

# dataframe that holds usernames, passwords and other user data
user_base <- data_frame(
  user = c("lgomez", "efernandez", "memorales", "aivalencia", "dgutierrez", "cfresno"),
  password = c("admin1234", "covid-inmegen-ez", "covid-inmegen-ms", "covid-inmegen-aa", "covid-inmegen-dz", "admin-system"), 
  password_hash = sapply(c("admin1234","covid-inmegen-ez", "covid-inmegen-ms", "covid-inmegen-aa", "covid-inmegen-dz", "admin-system"), sodium::password_store), 
  permissions = c("admin", "standard", "standard", "standard", "standard","standard"),
  name = c("Admin", "One", "Two", "Three", "Four", "Five")
)

saveRDS(user_base, "user_base.rds")


sidebar <- sidebarPanel(uiOutput("sidebarpanel"))


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  ######## ERROR MESSSAGE COLOR
  tags$head(
    tags$style(HTML("
                    .shiny-output-error-validation {
                    color: #ff0000;
                    font-weight: bold;
                    }
                    "))
  ),
  
  theme = shinytheme("spacelab"),
  
  #titlePanel( div(column(width = 6, tags$img(src = "images/inmegen.jpg"))),
  #                column(width = 6, h2("Aplicaci??n para an??lisis de RT-PCR")), 
  #            windowTitle="MyPage"),
  
  titlePanel( div(column(width = 6, h1("Aplicativo para analizar datos de RT-PCR")), 
                  column(width = 4, tags$img(src = "images/inmegen.jpg"))),
              windowTitle="rt-PCR-analysis"),
  
  hr(),
  
  ###### SIDE BAR - CONTROLADOR DE ACCIONES
  sidebar,
  ###### DISE??O DEL PANEL PARA IMPRESION DE RESULTADOS
  mainPanel(
    
    shinyjs::useShinyjs(),
    # add logout button UI 
    div(class = "pull-right", logoutUI(id = "logout")),
    # add login panel UI function
    loginUI(id = "login"),
    
    hr(),
    hr(),
    
    tabsetPanel(
      id = "navbar",
      tabPanel(title = "Tabla Resumen",
               value = "table",
               
               ###### PROTOCOLO
               fluidRow(
                 h4("Protocolo seleccionado"),
                 textOutput("protocolo")
               ),
               hr(),
               hr(),
               
               ###### ARCHIVO A PROCESAR    
               fluidRow(
                 h4("Directorio a procesar seleccionado"),
                 textOutput("input_dir")
               ),
               hr(),
               hr(),
               
               ###### DIRECTORIO DE RESULTADOS
               fluidRow(
                 h4("Directorio de resultados seleccionado"),
                 textOutput("output_dir")
               ),
               hr(),
               hr(),
               
               ###### TABLA DE RESULTADOS
               fluidRow( 
                 h4("Tabla de resultados"),
                 #textOutput("run_ready")
                 dataTableOutput(outputId = 'run_ready')
               )
      ),
      tabPanel(title = "Curvas QC",
               value = "curves", 
               
               ########## TABLA DE QC
               fluidRow( 
                 h4("Tabla de QC"),
                 br(),
                 br(),
                 #textOutput("run_ready")
                 dataTableOutput(outputId = 'qc_table')
                ),
               
               h3("Selecciona una placa para visualizar sus curvas de calidad"),
               selectInput(inputId = "plate_qc", label = "Placa", choices = NULL),
               br(),
               br(),
               
               h3(textOutput("caption1")),
               plotOutput("plot1"),
               br(),
               br(),
               br(),
               h2(textOutput("caption2")),
               plotOutput("plot2"),
               br(),
               br(),
               br(),
               h2(textOutput("caption3")),
               plotOutput("plot3")
               #dataTableOutput(outputId = 'summary_table_gene')
      ),
      tabPanel(title = "Curvas muestra",
               value = "samples", 
               
               h3("Selecciona una muestra para visualizar sus curvas de amplificacion"),
               selectInput(inputId = "sample", label = "", choices = NULL),
               br(),
               br(),
               
               h3(textOutput("caption_sample")),
               plotOutput("plot_sample")
      ),
      tabPanel(title = "Reporte de la placa",
               value = "plate",
               
               h3("Selecciona una placa para visualizar su configuracion"),
               selectInput(inputId = "conf", label = "", choices = NULL),
               br(),
               br(),
               
               dataTableOutput(outputId = 'plate_conf')
               #dataTableOutput(outputId = 'prueba'),
               #plotlyOutput('coveragePlot', height = "400px")
      )
    )
    
  )
)


###### SERVIDOR
server <- function(input, output, session) {
  
  credentials <- callModule(shinyauthr::login, "login", 
                            data = user_base,
                            user_col = user,
                            pwd_col = password_hash,
                            sodium_hashed = TRUE,
                            log_out = reactive(logout_init()))
  
  logout_init <- callModule(shinyauthr::logout, "logout", reactive(credentials()$user_auth))
  
  output$sidebarpanel <- renderUI({
    if(credentials()$user_auth) { 
      fluidPage( 
        sidebarPanel(
      ######## SELECT ONE GENE TO PLOT
      selectInput("protocol", "Selecciona un protocolo",
                  choices = as.list(c("CDC", "Berlin")), selected = 1),
      hr(),
      
      ######## h2("Selecciona el archivo a procesar"),
      #fileInput("rtpcr", "Sube el archivo a procesar"),
      h5('Selecciona el directorio a procesar'),
      shinyDirButton('input_dir', 'Directorio a procesar', 'Selecciona el directorio a procesar'),
      hr(),
      
      ######## h2("Selecciona el directorio para los resultados"),
      h5('Selecciona el directorio para almacenar los resultados'),
      shinyDirButton('directory', 'Directorio de resultados', 'Selecciona el directorio para almacenar los resultados'),
      #shinyDirChoose(input, 'out_dir', roots = c(home = '~')),
      #fileInput("dir_out", "Selecciona el directorio para almacenar los resultados"),
      hr(),
      
      ######## BUTTON TO GENERATE SUMMARY TABLE
      h5('Presiona para analizar los datos de la corrida'),
      actionButton("analizar", "Analizar corrida"),
      width = 12
    ))
    }
  })
  ####### LEER EL PROTOCOLO A UTILIZAR
  protocol <- reactive({
    protocol <- input$protocol
    
    if (is.null( protocol))
      return(NULL)
    
    return(protocol)
  })

  ####### IMPRIMIR EL PROTOCOLO A UTILIZAR
  output$protocolo <- renderText({
    protocol()
  })
  
  ###### LEER ESTRUCTURA DE DIRECTORIOS LOCAL
  volumes <- getVolumes()
  
  ###### DESPLIEGUE PARA LA ELECCION DEL DIRECTORIO A PROCESAR

  shinyDirChoose(input, 'input_dir', roots=volumes, session=session)
  
  input_dir <- reactive({
    return(print(parseDirPath(volumes, input$input_dir)))
  })
  
  ####### IMPRIMIR EL DIRECTORIO DE SALIDA
  output$input_dir <- renderText({
    input_dir()
  })
  
  ###### DESPLIEGUE PARA LA ELECCION DEL DIRECTORIO DE SALIDA
  
  shinyDirChoose(input, 'directory', roots=volumes, session=session)
  
  output_dir <- reactive({
    return(print(parseDirPath(volumes, input$directory)))
  })
  
  ####### IMPRIMIR EL DIRECTORIO DE SALIDA
  output$output_dir <- renderText({
    output_dir()
  })
  
  ####### CORRER EL PROCESO DE CLASIFICACION AL DARLE CLICK AL BOTON ANALIZAR
  table_out_list <- eventReactive(input$analizar, {
    
    rtpcr <- input_dir()
    
    if (is.null( rtpcr))
      return("NO EXISTE ARCHIVO DE ENTRADA")
    
    output <- output_dir()
    
    ######## VALIDATE THAT  INPUT DIRECTORY WAS SELECTED
    ######## OTHERWISE PRINT TEXT DESCRIBIING THE ERROR
    validate(
      need(rtpcr != "", "NO INPUT DIRECTORY WAS SELECTED")
    )
    
    
    ######## VALIDATE THAT OUTPUT DIRECTORY WAS SELECTED
    ######## OTHERWISE PRINT TEXT DESCRIBIING THE ERROR
    validate(
      need(output != "", "NO OUTPUT DIRECTORY WAS SELECTED")
    )
    
    ######## GET PATH FOR ALL EDS FILES IN DIRECTORY 
    eds.files <- list.files(rtpcr, pattern = ".eds")
    files <- paste(rtpcr, eds.files, sep="/")

    ######## GET PLATES NAMES DERIVED FROM EDS FILES IN DIRECTORY 
    plates <- sub(".eds", "", eds.files)
    
    ######## OBTAIN RESULT FOR ALL EDS FILES IN DIRECTORY    
    withProgress(message = 'corriendo analisis', value = 0.3, {
      all_results_list <- lapply(files, qpcr_pipeline.cdc, output=paste(output, "/", sep=""))
      names(all_results_list) <- plates
      
    })
    
    all_results_status <- sapply(all_results_list, is.list)
    all_results_valid <- all_results_list[all_results_status]
    
    return(all_results_valid)
    
  })
  
  ####### GENERAR TABLA DE RESULTADOS
  results_samples <- reactive({
    table_out_list()
    
    all_results_list <- table_out_list()
    ################################################################################
    #Merge results
    ################################################################################
    
    test_results <- lapply(all_results_list, function(x){
      df <- x$test_results
      df$ec.pass[df$ec.pass == "not_run"] <- NA
      return(df)
    })
    test_results_merge <- bind_rows(test_results)
    
    withProgress(message = 'imprimiendo reportes', value = 0.3, {
      makeReports(test_results_merge, paste(output_dir(), "/", sep=""))
    })
    
    return(test_results_merge)
  })
  
  ####### DESPLEGAR TABLA DE RESULTADOS
  output$run_ready <- renderDataTable({
    results_samples()
    
    datatable(results_samples()) %>% 
        formatStyle( 'classification', 
                     target = 'row',
                     backgroundColor = styleEqual(c("positive", "negative"),
                                                  c('aquamarine', 'pink')) ) %>% 
        formatStyle( 'qc', 
                   target = 'cell',
                   backgroundColor = styleEqual(c("FAIL"),
                                                c('yellow')) )
  })
  
  ####### DESPLEGAR TABLA DE RESULTADOS
  output$qc_table <- renderDataTable({
    table_out_list()
    
    all_results_list <- table_out_list()
    test_results_merge <- results_samples()
    ################################################################################
    #Merge QC
    ################################################################################
    
    qc_results <- lapply(all_results_list, function(x){x$qc_results$qc.values})
    qc_results_merge <- bind_rows(qc_results)
    
    ################################################################################
    #Add plate name 
    ################################################################################
    
    info.qc <- test_results_merge %>% 
      select(c(plate, ntc.pass, ptc.pass, ec.pass)) %>% 
      unique() %>% 
      mutate(rep = rowSums(.[2:4], na.rm = TRUE))
    
    qc_results_merge$plate <- rep(info.qc$plate, info.qc$rep)
    
    
    ################################################################################
    #Add QC classification
    ################################################################################
    
    pass.qc <- info.qc %>%
      select(c(ntc.pass, ptc.pass, ec.pass))  %>%
      as.data.frame()  %>%
      t() %>% 
      unlist() %>% 
      as.vector() %>% 
      na.omit() %>% 
      as.vector()
    
    qc_results_merge$QC <- pass.qc
    
    datatable(qc_results_merge) %>% 
      formatStyle( 'QC', 
                   target = 'row',
                   backgroundColor = styleEqual(c("TRUE", "FAIL"),
                                                c('white', 'yellow')) ) 
 })
    
  ################################################################################
  #Get plate names for drop-down menu
  ################################################################################
  observe({
    updateSelectInput(session = session, inputId = "plate_qc", choices = names(table_out_list()))
  })
  
  
  ####### IMPRIMIR CURVAS AL DARLE CLICK AL BOTON 
  
  ################################################################################
  # TABLE AND PLLOT PTC
  ################################################################################
  
  output$caption1 <- renderText({
    table_out_list()
    plate <- input$plate_qc
    names(table_out_list()[[plate]]$triplets.qc)[1]
  })
  
  output$plot1 <- renderPlot({
    table_out_list()
    plate <- input$plate_qc
    plots <- table_out_list()[[plate]]$triplets.qc
    (plots[[1]] + ggtitle(labels(plots)[1]))
  })
  
  ################################################################################
  # TABLE AND PLLOT NTC
  ################################################################################
  
  output$caption2 <- renderText({
    table_out_list()
    plate <- input$plate_qc
    names(table_out_list()[[plate]]$triplets.qc)[2]
  })
  
  output$plot2 <- renderPlot({
    table_out_list()
    plate <- input$plate_qc
    plots <- table_out_list()[[plate]]$triplets.qc
    (plots[[2]] + ggtitle(labels(plots)[2]))
  })
  
  ################################################################################
  # TABLE AND PLLOT EC
  ################################################################################
  
  output$caption3 <- renderText({
    table_out_list()
    plate <- input$plate_qc
    text <- "ESTA PLACA NO CONTIENE EL CONTROL EC. REVISAR LA PLACA 4 CONSECUTIVA."
    if (length(names(table_out_list()[[plate]]$triplets.qc)) == 3){
      text <- names(table_out_list()[[plate]]$triplets.qc)[3]
    }
    text
  })
  
  output$plot3 <- renderPlot({
    table_out_list()
    plate <- input$plate_qc
    plots <- table_out_list()[[plate]]$triplets.qc
    if (length(plots) == 3){
      (plots[[3]] + ggtitle(labels(plots)[3]))
    }
  })
  
  
  ################################################################################
  #Get sample names for drop-down menu
  ################################################################################
  observe({
    updateSelectInput(session = session, inputId = "sample", choices = results_samples()$sample)
  })
  
  output$caption_sample <- renderText({
    results_samples()
    merge <- results_samples()
    sample <- input$sample
    sample
  })
  
  output$plot_sample <- renderPlot({
    results_samples()
    merge <- results_samples()
    sample <- input$sample
    plate <- merge$plate[merge$sample == sample]
    plots <- table_out_list()[[plate]]$triplets.samples[[sample]]
    (plots)
  })
  
  
  observe({
    updateSelectInput(session = session, inputId = "conf", choices = names(table_out_list()))
  })
  
  ####### OBTENER LA CONFIGURACION DE LA PLACA AL DARLE CLICK AL BOTON REPORTE PLACA
  plate_out <- reactive({
    table_out_list()
    
    input_dir <- input_dir()
    
    if (is.null( input_dir))
      return("NO EXISTE EL DIRECOTRIO DE ENTRADA")
    
    eds_file <- paste(paste(input_dir, input$conf, sep="/"), "eds", sep=".")
    
    plate_conf <- get_plate_conf(eds_file)
    return(plate_conf)
  })
   
  ####### IMPRIMIR TABLA DE CONFUGURACION 
  output$plate_conf <- renderDataTable({
    plate_out()
    datatable(plate_out())
  })

}

# Run the application 
shinyApp(ui = ui, server = server)

