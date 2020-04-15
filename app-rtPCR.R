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


source("function_app_qPCR_cdc.R")
source("function_plate_conf.R")
source("src/functions.R")

#input <- c("data/Procolo_COVID-19_Prueba1_4Abr20.eds")
#output <- c("results/")

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  
  #titlePanel( div(column(width = 6, tags$img(src = "images/inmegen.jpg"))),
  #                column(width = 6, h2("Aplicaci??n para an??lisis de RT-PCR")), 
  #            windowTitle="MyPage"),
  
  titlePanel( div(column(width = 6, h1("Aplicativo para analizar datos de RT-PCR")), 
                  column(width = 4, tags$img(src = "images/inmegen.jpg"))),
              windowTitle="rt-PCR-analysis"),
  
  hr(),
  
  ###### SIDE BAR - CONTROLADOR DE ACCIONES
  sidebarPanel(
    ######## SELECT ONE GENE TO PLOT
    selectInput("protocol", "Selecciona un protocolo",
                choices = as.list(c("CDC", "Berlin")), selected = 1),
    hr(),
    
    ######## h2("Selecciona el archivo a procesar"),
    #fileInput("rtpcr", "Sube el archivo a procesar"),
    h5('Selecciona el archivo a procesar'),
    shinyFilesButton('file', 'Archivo a procesar', 'Selecciona el archivo a procesar', FALSE),
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
    hr(),
    
    ######## BUTTON TO PRINT PLOTS TO WEB
    #h5('Presiona para imprimir las curvas control en la interfaz web'),
    #actionButton("plot", "Imprimir curvas control en Web"),
    #hr(),
    
    ######## BUTTON TO GENERATE PLATE CONFIGURATION TABLE
    h5('Presiona para generar el reporte de la placa'),
    actionButton("plate", "Reporte placa")

  ),
  
  ###### DISE??O DEL PANEL PARA IMPRESION DE RESULTADOS
  mainPanel(
    
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
                 h4("Archivo a procesar seleccionado"),
                 textOutput("input_file")
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
      tabPanel(title = "Curvas control",
               value = "curves", 
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
      tabPanel(title = "Reporte de la placa",
               value = "plate", 
               dataTableOutput(outputId = 'plate_conf')
               #dataTableOutput(outputId = 'prueba'),
               #plotlyOutput('coveragePlot', height = "400px")
      )
    )
    
  )
)


###### SERVIDOR
server <- function(input, output, session) {
  
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
  
  ###### DESPLIEGUE PARA LA ELECCION DEL ARCHIVO A PROCESAR
  shinyFileChoose(input,'file', roots=volumes, session=session)
  
  input_file <- reactive({
    inFile <- parseFilePaths(volumes, input$file)
    inFile.path <- as.character(inFile$datapath)
  })
  
  ####### IMPRIMIR LA RUTA DEL ARCHIVO A PROCESAR
  output$input_file <- renderText({
    input_file()
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
  table_out <- eventReactive(input$analizar, {
    
    rtpcr <- input_file()
    
    if (is.null( rtpcr))
      return("NO EXISTE ARCHIVO DE ENTRADA")
    
    output <- output_dir()
    
    withProgress(message = 'corriendo analisis', value = 0.3, {
      all_results <- qpcr_pipeline.cdc(input=rtpcr, output=paste(output, "/", sep=""))
    })
    
    return(all_results)
    
  })
  
  ####### DESPLEGAR TABLA DE RESULTADOS
  output$run_ready <- renderDataTable({
    table_out()
    
    qc <- table_out()$qc_results
    
    if (qc$QC != "PASS"){
      datatable(table_out()$test_results) %>% 
        formatStyle('N1', 
                    target='row',
                    backgroundColor = "yellow" )
    }else{
      datatable(table_out()$test_results) %>% 
        formatStyle( 'classification', 
                     target = 'row',
                     backgroundColor = styleEqual(c("positive", "negative"),
                                                  c('aquamarine', 'pink')) )
    }
  })
  
  ####### IMPRIMIR CURVAS AL DARLE CLICK AL BOTON 
  
  output$caption1 <- renderText({
    table_out()
    names(table_out()$triplets.qc)[1]
  })
  
  output$plot1 <- renderPlot({
    table_out()
    plots <- table_out()$triplets.qc
    (plots[[1]] + ggtitle(labels(plots)[1]))
  })
  
  output$caption2 <- renderText({
    table_out()
    names(table_out()$triplets.qc)[2]
  })
  
  output$plot2 <- renderPlot({
    table_out()
    plots <- table_out()$triplets.qc
    (plots[[2]] + ggtitle(labels(plots)[2]))
  })
  
  output$caption3 <- renderText({
    table_out()
    names(table_out()$triplets.qc)[3]
  })
  
  output$plot3 <- renderPlot({
    table_out()
    plots <- table_out()$triplets.qc
    (plots[[3]] + ggtitle(labels(plots)[3]))
  })
  
  
  
  ####### OBTENER LA CONFIGURACION DE LA PLACA AL DARLE CLICK AL BOTON REPORTE PLACA
  plate_out <- eventReactive(input$plate, {
    
    rtpcr <- input_file()
    
    if (is.null( rtpcr))
      return("NO EXISTE ARCHIVO DE ENTRADA")
    
    plate_conf <- get_plate_conf(rtpcr)
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

