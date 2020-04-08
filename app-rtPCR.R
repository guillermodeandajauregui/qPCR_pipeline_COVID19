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


source("function_qpcrPipeline.R")

input <- c("data/Procolo_COVID-19_Prueba1_4Abr20.eds")
output <- c("results/")

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  
  #titlePanel( div(column(width = 6, tags$img(src = "images/inmegen.jpg"))),
  #                column(width = 6, h2("Aplicación para análisis de RT-PCR")), 
  #            windowTitle="MyPage"),
  
  titlePanel( div(column(width = 6, h1("Aplicación para análisis de RT-PCR")), 
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
    actionButton("analizar", "Analizar corrida")

  ),
  
  ###### DISEÑO DEL PANEL PARA IMPRESION DE RESULTADOS
  mainPanel(
    
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
      qpcr_pipeline.cdc(rtpcr, paste(output, "/", sep=""))
      table_summary <- read.table(paste(output, "result_table.txt", sep="/"), header=T)
    })
    
    return(table_summary)
    
  })
  
  ####### DESPLEGAR TABLA DE RESULTADOS
  output$run_ready <- renderDataTable({
    table_out()
    datatable(table_out())
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

