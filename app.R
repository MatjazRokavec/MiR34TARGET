library(shiny)
library(shinydashboard)
library(shinyjs)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
#library(msigdbr)
#library(clusterProfiler)
library(DT)


datasetmiR34a <- fread("datasets/miR34a_GEO_FC_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34a_KO <- fread("datasets/miR34a_mouse_KO_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34a_SILAC <- fread("datasets/miR34a_SILAC_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34a_PULLDOWN <- fread("datasets/miR34a_PULLDOWN_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34_correlation <- fread("datasets/miR34_correlation_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34a_correlation <- fread("datasets/miR34a_correlation_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34b_correlation <- fread("datasets/miR34b_correlation_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34c_correlation <- fread("datasets/miR34c_correlation_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)

datasetmiR34a$Name <- factor(datasetmiR34a$Name, levels = datasetmiR34a$Name)
datasetmiR34a_KO$Name <- factor(datasetmiR34a_KO$Name, levels = datasetmiR34a_KO$Name)
datasetmiR34a_SILAC$Name <- factor(datasetmiR34a_SILAC$Name, levels = datasetmiR34a_SILAC$Name)
datasetmiR34a_PULLDOWN$Name <- factor(datasetmiR34a_PULLDOWN$Name, levels = datasetmiR34a_PULLDOWN$Name)
datasetmiR34_correlation$order <- factor(datasetmiR34_correlation$order, levels = datasetmiR34_correlation$order)
datasetmiR34a_correlation$order <- factor(datasetmiR34a_correlation$order, levels = datasetmiR34a_correlation$order)
datasetmiR34b_correlation$order <- factor(datasetmiR34b_correlation$order, levels = datasetmiR34b_correlation$order)
datasetmiR34c_correlation$order <- factor(datasetmiR34c_correlation$order, levels = datasetmiR34c_correlation$order)

datasetmiR34avalidate <- fread("datasets/miR34a_validated_targets_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34avalid <- fread("datasets/miR34a_validated_targets.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)

# load summarydata
summarymiR34a <- fread("datasets/miR34a datasets/miR34a_datasets_description.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
summarymiR34aKO <- fread("datasets/miR34a datasets/miR34a_datasets_KO_description.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
summarymiR34a <- fread("datasets/miR34a datasets/miR34a_datasets_description.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)

#single study miR34a OE datasets
datasetGSE7678 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE7678.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE7754 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE7754.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE16674 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE16674.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE21832 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE21832.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE34242 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE34242.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE38341 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE38341.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE49845 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE49845.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE68740 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE68740.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE71423 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE71423.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE89156 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE89156.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE98601 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE98601.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE117506 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE117506.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE185946 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE185946.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE99401 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE99401.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE58892 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE58892.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE58800 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE58800.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE10455 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE10455.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE41322 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE41322.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE57493 <- fread("datasets/miR34a datasets/GEO, miR34a OE/GSE57493.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)


#single study miR34a KO mouse datasets
datasetGSE123628c <- fread("datasets/miR34a datasets/miR34a ko mice/GSE123628, colon cells_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE123628t <- fread("datasets/miR34a datasets/miR34a ko mice/GSE123628, T-cells_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE69484 <- fread("datasets/miR34a datasets/miR34a ko mice/GSE69484_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE92296 <- fread("datasets/miR34a datasets/miR34a ko mice/GSE92296_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE99452 <- fread("datasets/miR34a datasets/miR34a ko mice/GSE99452, Gulfem, 6xAOM_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE133775 <- fread("datasets/miR34a datasets/miR34a ko mice/GSE133775_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#datasetMAC <- fread("datasets/miR34a datasets/miR34a ko mice/Gulfem, macrophages_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetGSE84138 <- fread("datasets/miR34a datasets/miR34a ko mice/GSE84138_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)

#single study TCGA correlation datasets
#BLCA <- fread("datasets/miR34a datasets/correlation/correl, TCGA BLCA, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#BRCA <- fread("datasets/miR34a datasets/correlation/correl, TCGA BRCA, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#CESC <- fread("datasets/miR34a datasets/correlation/correl, TCGA CESC, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#CRC <- fread("datasets/miR34a datasets/correlation/correl, TCGA CRC, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#ESCA <- fread("datasets/miR34a datasets/correlation/correl, TCGA ESCA, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#LIHC <- fread("datasets/miR34a datasets/correlation/correl, TCGA LIHC, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#PAAD <- fread("datasets/miR34a datasets/correlation/correl, TCGA PAAD, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
#STAD <- fread("datasets/miR34a datasets/correlation/correl, TCGA STAD, log2_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)

#load miRNA target prediction data
datasetmiR34atarget <- fread("datasets/miRNA target prediction/miR34a_targets_prediction_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34btarget <- fread("datasets/miRNA target prediction/miR34b_targets_prediction_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34ctarget <- fread("datasets/miRNA target prediction/miR34c_targets_prediction_trans.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
datasetmiR34atarget$Name <- factor(datasetmiR34atarget$Name, levels = datasetmiR34atarget$Name)
datasetmiR34btarget$Name <- factor(datasetmiR34btarget$Name, levels = datasetmiR34btarget$Name)
datasetmiR34ctarget$Name <- factor(datasetmiR34ctarget$Name, levels = datasetmiR34ctarget$Name)

#load pre-calculated scores required for the screener
precalcGEOOE <- fread("datasets/precalc/miR34a_GEO_FC_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
precalcGEOKO <- fread("datasets/precalc/miR34a_mouse_KO_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
precalcPREDICT <- fread("datasets/precalc/miR34_targets_prediction_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
precalcCORRELa <- fread("datasets/precalc/miR34a_correlation_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
precalcCORRELb <- fread("datasets/precalc/miR34b_correlation_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
precalcCORRELc <- fread("datasets/precalc/miR34c_correlation_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
precalcSILAC <- fread("datasets/precalc/miR34a_SILAC_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)
precalcPULLDOWN <- fread("datasets/precalc/miR34a_PULLDOWN_calc.txt", sep="\t", dec=".", check.names=FALSE, stringsAsFactors=FALSE)


#if you click Enter, actionbutton is clicked; https://stackoverflow.com/questions/55731607/trigger-actionbutton-with-keyboard-in-shiny-using-uioutput
jscode <- '$(document).keyup(function(e) {
    if (e.key == "Enter") {
    $("#go").click();
}});'

font_size <- 11
header <- dashboardHeader(disable=TRUE)
sidebar <- dashboardSidebar(disable=TRUE)


body <- dashboardBody(
  
  tags$head(
    tags$style(
      HTML(".shiny-notification {
             position:fixed;
             top: calc(6%);
             left: calc(6%);
             }
             "
      )
    )
  ),
  
  
  tags$head(tags$script(HTML(jscode))), #if you click Enter, actionbutton is clicked
  
  tags$head(tags$style(HTML(
  #background color white instead of default gray
    '.content-wrapper {
        background-color: #fff;
      }
    '))),
  
  tabsetPanel(id="Tabset",
    tabPanel(title="query genes", value="searchgenes",
        useShinyjs(),    ## IMPORTANT: so shiny knows to use the shinyjs library    
        br(),
        fluidRow(
                  column(2, textInput(inputId = "gene1",
                             label = NULL,
                             placeholder = "Enter a gene symbol",
                             width="auto")
                         ),
                  
                  #column(3,   textInput.typeahead(             #package ShinySky; suggests gene names; doesn't work so well
                  #                               id = "gene1",
                  #                               label = NULL,
                  #                               placeholder = "Enter a human gene symbol",
                  #                               local = precalcGEOOE$Name, 
                  #                               width="auto")
                  #       ),
                  
                  column(1, actionButton("go", "search", class = "btn-primary")),
                  column(6, textOutput("instructions")),
                  column(2, htmlOutput("validated")),
                  column(1, htmlOutput("depmap"))
        ),
        br(),
        fluidRow(
                  column(width=4, htmlOutput("FCplot_title"), #GEO FC mRNA
                  plotOutput("FCPlot", click = "FCplot_click", height = "auto")
                        ),
                  
                  column(width=5,
                           fluidRow(
                    column(width=5, htmlOutput("predictPlot_title"), #prediction 
                         plotOutput("predictPlot", width = "auto", height = "auto"),
                          ),
                    column(width=7, htmlOutput("FCplot_SILAC_title"), #SILAC and pulldown
                           plotOutput("FCPlot_SILAC", width = "auto", height = "auto"),
                           br(),    
                           htmlOutput("FCplot_PULLDOWN_title"),
                           plotOutput("FCPlot_PULLDOWN", width = "auto", height = "auto"),
                      )),
                    fluidRow(
                      br(),
                    column(width=12, 
                           htmlOutput("miR34a_KO_mice_title"), #mice KO FC mRNA
                           plotOutput("KOPlot", click = "FCplotKO_click", width = "auto", height = "auto")
                    )),
                    
                    ),
                    
                  
                  column(width=3, htmlOutput("correlation_title"), #correlation    
                         #plotOutput("correlationPlot", click = "Correlation_click", width = "auto", height = "auto")
                         plotOutput("correlationPlot", width = "auto", height = "auto")
                  ),
                  
                  
        ),
        br(), br(),
        fluidRow(style="height:1px;",
              column(10, tableOutput("summaryFC")),
              column(2, plotOutput("RawDataFC"))
        ),

        br(), br(),
        fluidRow(
              column(10, tableOutput("summaryFCKO")),
              column(2, plotOutput("RawDataFCKO")),
        ),
        
        #fluidRow(
        #  column(3, plotOutput("CorrelationPlot34a")),
        #  column(3, plotOutput("CorrelationPlot34b")),
        #  column(3, plotOutput("CorrelationPlot34c"))
        #        ),
    ),
    
    tabPanel("screen for miR-34 targets", fluid = TRUE, #SCREENER
             fluidPage(

               tags$style(".checkbox {margin-top:0px; margin-bottom:-15px;}"), #top/bottom margin of checkboxes set to 0
               
               sidebarLayout(    
                 sidebarPanel(width=5,
                              wellPanel(  #GEO FC miR34 induction
                                #checkboxInput("GEOOEcheck", "Change in mRNA expression after miR-34 induction", FALSE),
                              
                                div(style="display: flex; flex-wrap: wrap; align-items: flex-start;",  #put items side by side; https://stackoverflow.com/questions/20637248/shiny-4-small-textinput-boxes-side-by-side; https://css-tricks.com/snippets/css/a-guide-to-flexbox/
                                    div(style = "flex: 1;", checkboxInput("GEOOEcheck", NULL, FALSE)), #checkbox
                                    div(style = "flex: 20;", HTML('<b
                                                                  style="position: relative;top: 0px;">Change in mRNA expression by ectopic miR-34</b>')) #align text to checkbox; https://community.rstudio.com/t/vertical-align-checkboxinput-and-its-label/71686
                                ),
                                
                                div(style="display: flex; flex-wrap: wrap;", #put two sliders side by side
                                    div(style = "flex: 5;", p("Repression more than X fold:"),
                                        sliderInput(inputId = "GEOOEsliderFC",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 4,
                                                    value = 1.5,
                                                    step = 0.25)),
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 5;", p("In at least X studies:"),
                                        sliderInput(inputId = "GEOOEsliderNO",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 27,
                                                    value = 10)),
                                ),
                                
                                br(),     #line separator
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 1;", checkboxInput("GEOKOcheck", "", FALSE)),#GEO FC miR34 mouse KO
                                    div(style = "flex: 20;", HTML('<b
                                                      style="position: relative;top: 0px;">Change in mRNA expression in <i> miR-34 </i> knockout mice</b>'))
                                ),
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 5;", p("Induction more than X fold:"),
                                        sliderInput(inputId = "GEOKOsliderFC",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 4,
                                                    value = 1.5,
                                                    step = 0.25)),
                                    
                                    div(style = "flex: 1;", p("")),
                                    
                                    
                                    div(style = "flex: 5;", p("In at least X studies:"),
                                        sliderInput(inputId = "GEOKOsliderNO",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 7,
                                                    value = 3)),
                                ),
                                
                                br(), 
                                
                                tags$b("Predicted miR-34 target in at least X algorithms"), #Target Prediction, bold text
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 1;", checkboxInput("predictioncheckA", "miR-34a", FALSE)), #prediction 34a
                                    div(style = "flex: 5;",
                                        sliderInput(inputId = "predictAslider",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 12,
                                                    step = 1,
                                                    value = 6))
                                ),
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 1;", checkboxInput("predictioncheckB", "miR-34b", FALSE)), #prediction 34b
                                    div(style = "flex: 5;",
                                        sliderInput(inputId = "predictBslider",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 12,
                                                    step = 1,
                                                    value = 6))
                                ),
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 1;", checkboxInput("predictioncheckC", "miR-34c", FALSE)), #prediction 34c
                                    div(style = "flex: 5;",
                                        sliderInput(inputId = "predictCslider",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 12,
                                                    step = 1,
                                                    value = 6))
                                ),
                                
                                br(),     #line separator
                                
                                tags$b("Correlation with miR-34 expression (number of cancer types)"), #Correlation, bold text
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 3;", checkboxInput("correlAcheck", "miR-34a", FALSE)), #prediction miR34a
                                    div(style = "flex: 1;", p("")),                                    
                                    div(style = "flex: 10;", p("Correlation coeficient:"),
                                        sliderInput(inputId = "correlAsliderR",
                                                    label = NULL,
                                                    min = -1,
                                                    max = 0,
                                                    value = -0.1,
                                                    step = 0.1)),
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 8;", p("Number of cancer types:"),
                                        sliderInput(inputId = "correlAsliderNO",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 32,
                                                    value = 10,
                                                    step=1))
                                ),
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 3;", checkboxInput("correlBcheck", "miR-34b", FALSE)), #prediction miR34b
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 10;", p("Correlation coeficient:"),
                                        sliderInput(inputId = "correlBsliderR",
                                                    label = NULL,
                                                    min = -1,
                                                    max = 0,
                                                    value = -0.1,
                                                    step = 0.1)),
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 8;", p("Number of cancer types:"),
                                        sliderInput(inputId = "correlBsliderNO",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 32,
                                                    value = 10,
                                                    step=1))
                                ),
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 3;", checkboxInput("correlCcheck", "miR-34c", FALSE)), #prediction miR34c
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 10;", p("Correlation coeficient:"),
                                        sliderInput(inputId = "correlCsliderR",
                                                    label = NULL,
                                                    min = -1,
                                                    max = 0,
                                                    value = -0.1,
                                                    step = 0.1)),
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 8;", p("Number of cancer types:"),
                                        sliderInput(inputId = "correlCsliderNO",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 32,
                                                    value = 10,
                                                    step=1))
                                ),
                                
                                br(),     #line separator
                                
                                
                                div(style="display: flex; flex-wrap: wrap;", #SILAC miR-34a OE
                                    div(style = "flex: 1;", checkboxInput("SILACcheck", "", FALSE)),
                                    div(style = "flex: 20;", HTML('<b
                                                      style="position: relative;top: 5px;">Change in protein expression after miR-34a induction</b>'))
                                ),
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 5;", p("Repression more than X fold:"),
                                        sliderInput(inputId = "SILACsliderFC",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 4,
                                                    value = 1.5,
                                                    step = 0.25)),
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 5;", p("In at least X studies:"),
                                        sliderInput(inputId = "SILACsliderNO",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 3,
                                                    value = 1,
                                                    step =1)),
                                ),
                                
                                br(),     #line separator
                                
                                
                                
                                div(style="display: flex; flex-wrap: wrap;",  #PULLDOWN miR-34a
                                    div(style = "flex: 1;", checkboxInput("PULLDOWNcheck", "", FALSE)),
                                    div(style = "flex: 20;", HTML('<b
                                                      style="position: relative;top: 5px;">Enrichment of mRNA in miR-34a pulldown studies</b>'))
                                ),               
                                
                                div(style="display: flex; flex-wrap: wrap;",
                                    div(style = "flex: 5;", p("Enrichment more than X fold:"),
                                        sliderInput(inputId = "PULLDOWNFC",
                                                    label = NULL,
                                                    min = 2,
                                                    max = 10,
                                                    value = 2,
                                                    step = 2)),
                                    div(style = "flex: 1;", p("")),
                                    div(style = "flex: 5;", p("In at least X studies:"),
                                        sliderInput(inputId = "PULLDOWNsliderNO",
                                                    label = NULL,
                                                    min = 1,
                                                    max = 3,
                                                    value = 1)),
                                ),
                                
                                br(),     #line separator
                                
                                actionButton("reset", "reset", class = "btn-primary")
                                
                              )),
                 
                 mainPanel(width=7,
                   fluidRow(
                     column(width=8, 
                            br(),
                            htmlOutput("INTERSECTNOgenes"),
                            br(),
                            DT::dataTableOutput("INTERSECTresults")),
                     #column(br(),br(),width=2, downloadButton("downloadgenes", "download")), #button for download
                     #column(width=6,
                     #      br(),br(),
                     #      plotOutput("GSEA_H"),
                     #      br(),
                     #      plotOutput("GSEA_KEGG"),
                     #      br(),
                     #     plotOutput("GSEA_GOBP"))

                   )
                 )
               )                   
             )
    ),
    
    
    tabPanel("download datasets", fluid = TRUE,
           br(),
           downloadLink("DLmiR34OE", "Datasets with miR34 induction in cell lines (mRNA expression)"),
           br(),
           downloadLink("DLmiR34silac", "Datasets with miR34 induction in cell lines (protein expression)"),
           br(),
           downloadLink("DLmiR34KO", "Datasets with miR34 knockout mouse models (mRNA expression)"),
           br(),
           downloadLink("DLmiR34pulldown", "Datasets with miR34 pulldown in cell lines (mRNA expression)"),
           br(),
           downloadLink("DLmiR34prediction", "Predicted miR34 target genes (based on miRNA prediction algorithms)"),
           br(),
           downloadLink("DLmiR34correl", "Correlation betweenn miR34 and mRNAs (based on TCGA dataset)"),
           br(),
           downloadLink("DLmiR34validate", "Published miR34 target genes (updated on January, 2022)")
           
    ),
    
    tabPanel("about", fluid = TRUE,
             br(),
             #downloadLink("manual", "Tutorial, how to use this website (PDF)"),
             #br(), br(),
             textOutput(outputId = "about"),
             br(), br(),
             uiOutput(outputId = "contact1"),
             textOutput(outputId = "contact2"),   
             uiOutput(outputId = "contact3"),   
             textOutput(outputId = "contact4"),   
             textOutput(outputId = "contact5"),   
             textOutput(outputId = "contact6"),   
             
    )
  )
)

ui <- dashboardPage(
  header,
  sidebar,
  body
)

server = function(input, output, session) {

  #search server
  
  #instructions
  output$instructions <- renderText("The platform will report: the fold change in mRNA (A) and protein (C) expression of the queried gene after miR34a/b/c induction. (B) Prediction as miR34a/b/c target gene. (D) mRNA enrichment in miR34a pulldown studies. (E) the fold change in mRNA expression in miR34a/b/c knockout mice. (F) the correlation between the mRNA expression of the queried gene and the expression of miR34a/b/c in TCGA datasets. Click on individual bars in panels A and E to obtain more information about indicated studies.")  
  
  
  observeEvent(input$go, {
    
    gene <- toupper(input$gene1)
    
    if (!(gene %in% precalcGEOOE$Name)) showNotification("Enter a valid gene symbol", duration = 2, type = "error")
    
    req(gene %in% precalcGEOOE$Name)
    
    shinyjs::hide(id = "summaryFCKO")
    shinyjs::hide(id = "RawDataFCKO")
    #shinyjs::hide(id = "CorrelationPlot34a")   
    #shinyjs::hide(id = "CorrelationPlot34b")
    #shinyjs::hide(id = "CorrelationPlot34c")
    shinyjs::hide(id = "summaryFC")
    shinyjs::hide(id = "RawDataFC")

    miR34aFC <- datasetmiR34a[[gene]]
    miR34aKOFC <- datasetmiR34a_KO[[gene]]
    miR34atargetFC <- datasetmiR34atarget[[gene]]
    miR34avalidate <-datasetmiR34avalidate[[gene]]
    
    Name <- datasetmiR34a$Name
    ID <- datasetmiR34a$ID
    V1 <- datasetmiR34a$order
    miRNA <- datasetmiR34a$miRNA
    miR34adata <- data.frame(V1, Name, ID, miR34aFC, miRNA)
    datanoNA <- miR34adata[complete.cases(miR34adata[ , 4]),]
    datanoNA$ID <- as.character(datanoNA$ID)
    
    NameKO <- datasetmiR34a_KO$Name
    IDKO <- datasetmiR34a_KO$ID
    V1KO <- datasetmiR34a_KO$order
    miR34adataKO <- data.frame(V1KO, NameKO, IDKO, miR34aKOFC)
    dataKOnoNA <- miR34adataKO[complete.cases(miR34adataKO[ , 4]),]
    dataKOnoNA$IDKO <- as.character(dataKOnoNA$IDKO)    

    miR34acorrel <- datasetmiR34_correlation[[gene]]
    
    Name_SILAC <- datasetmiR34a_SILAC$Name
    miR34aFC_SILAC <- datasetmiR34a_SILAC[[gene]]
    miR34adata_SILAC <- data.frame(Name_SILAC, miR34aFC_SILAC)
    datanoNA_SILAC <- miR34adata_SILAC[complete.cases(miR34adata_SILAC[ , 2]),]
    
    Name_PULLDOWN <- datasetmiR34a_PULLDOWN$Name
    miR34aFC_PULLDOWN <- datasetmiR34a_PULLDOWN[[gene]]
    miR34adata_PULLDOWN <- data.frame(Name_PULLDOWN, miR34aFC_PULLDOWN)
    datanoNA_PULLDOWN <- miR34adata_PULLDOWN[complete.cases(miR34adata_PULLDOWN[ , 2]),]
    
    summaryfirstrow <- summarymiR34a$ID
    summaryfirstrowKO <- summarymiR34aKO$ID

    #chart titles, first row
    output$FCplot_title <- renderText({paste("<b>(A) Change in", "<i>", gene,"</i>mRNA expression after the induction of", "<font color=\"darkblue\">miR-34a","<font color=\"darkgreen\">, miR-34b", "<font color=\"black\">, or", "<font color=\"firebrick\">miR-34c</b>")})
    output$miR34a_KO_mice_title <- renderText({paste("<b>(E) Change in", "<i>", gene,"</i>mRNA expression in miR-34 knockout mice</b>")})
    output$predictPlot_title <- renderText({paste("<b>(B) prediction of", "<i>", gene, "</i>as a miR-34 target")})
    output$correlation_title <- renderText({paste("<b>(F) Correlation between", "<i>", gene, "</i>mRNA and miR-34a/b/c expression</b>")})
    
    #barplot of FC GEO mRNA  
    output$FCPlot <- renderPlot({
      ggplot(datanoNA, aes(x=miR34aFC, y=Name, fill=miRNA)) + 
        geom_bar(stat = "identity", color="black", width=0.75) +
        scale_fill_manual(values=c("darkblue", "darkgreen", "firebrick"))+
        scale_x_continuous(trans='log2') +
        theme(text=element_text(size=font_size)) +
        theme(axis.text.y=element_text(size=font_size))+
        geom_vline(xintercept = 1, color = "black") +
        theme(panel.background = element_rect(fill = "white", color="black")) +
        labs(x=paste (gene, "mRNA expression (fold change)")) +
        labs(y="Study\n") +
        theme(legend.position="none")
    }, height=400) #530
    
    #barplot of mice KO GEO mRNA  
    output$KOPlot <- renderPlot({
      ggplot(dataKOnoNA, aes(x=miR34aKOFC, y=NameKO)) + 
        geom_bar(stat = "identity", color="black", fill="darkblue", width=0.75) +
        scale_x_continuous(trans='log2') +
        theme(text=element_text(size=font_size)) +
        theme(axis.text.y=element_text(size=font_size)) +
        geom_vline(xintercept = 1, color = "black") +
        theme(panel.background = element_rect(fill = "white", color="black")) +
        labs(x=paste (gene, "mRNA expression (fold change)")) +
        labs(y="Study\n")
    }, height=180)
    
    #Correlationplot mRNA/miR-34 correlation  
    
    #DOT plot
    'output$correlationPlot <- renderPlot({
      ggplot(datasetmiR34_correlation, aes(x=Name, y=miR34acorrel)) + 
        geom_dotplot(aes(fill=miRNA), binwidth=(max(miR34acorrel)-min(miR34acorrel))/60, binaxis="y", stackdir="center")+ #binwidth should not be a constant number, otherwise it will change when expression values (y-axis) change
        coord_flip()+
        scale_fill_manual(values=c("darkblue", "darkgreen", "firebrick"))+
        theme(text=element_text(size=font_size)) +
        theme(axis.text.y=element_text(size=font_size))+
        theme(axis.text.x=element_text(size=font_size))+
        geom_hline(yintercept = 0, color = "black") +
        theme_minimal()+
        labs(x="Cancer type (TCGA)\n") +
        labs(y="Pearson correlation (r)")+
        theme(legend.position="none")
    }, height=170)'
    
    #heatmap correlation
    output$correlationPlot <- renderPlot({
      correlName <- datasetmiR34a_correlation$Name
      miR34a <- datasetmiR34a_correlation[[gene]]
      miR34b <- datasetmiR34b_correlation[[gene]]
      miR34c <- datasetmiR34c_correlation[[gene]]
      correlHeatMap <- data.frame(correlName, miR34a, miR34b, miR34c)
      correlHeatMap$correlName <- factor(correlHeatMap$correlName, levels = correlHeatMap$correlName[order(correlHeatMap$correlName, decreasing = TRUE)])
      melted_correlHeatMap <- reshape2::melt(correlHeatMap)          

      ggplot(melted_correlHeatMap, aes(x = variable, y = correlName, fill = value)) +
        geom_tile(aes(fill = value), colour = "gray90", lwd=1) +
        theme(text=element_text(size=font_size)) +
        theme(axis.text.y=element_text(size=font_size))+
        theme(axis.text.x=element_text(size=font_size))+
        theme(panel.background = element_rect(fill = "white", color="black")) +
        labs(x=NULL) +
        labs(y="Cancer type (TCGA)\n")+
        #theme(legend.position = "none")+
        labs(fill="Pearson\ncorrelation\ncoeficient")+
        scale_x_discrete(position = "top")+
        scale_fill_gradient2(low = "blue",
                             mid = "white",
                             high = "red")
    },height=400)
    
    #SILAC chart title
    output$FCplot_SILAC_title <- renderText({
      if (nrow(datanoNA_SILAC)<1) {paste("<b>(C) Change in", gene,"protein expression after the induction of miR-34a:", gene,"protein not detected in any study")}
      else {paste("<b>(C) Change in", gene,"protein expression after the induction of miR-34a")}})
    
    #barplot of FC SILAC  
    output$FCPlot_SILAC <- renderPlot({
      if (nrow(datanoNA_SILAC)<1) {}
      else {
        ggplot(datanoNA_SILAC, aes(x=miR34aFC_SILAC, y=Name_SILAC)) + 
          geom_bar(stat = "identity", color="black", fill="darkblue", width=0.75) +
          scale_x_continuous(trans='log2') +
          theme(text=element_text(size=font_size)) +
          theme(axis.text.y=element_text(size=font_size))+
          geom_vline(xintercept = 1, color = "black") +
          theme(panel.background = element_rect(fill = "white", color="black")) +
          labs(x=paste (gene, "protein expression (fold change)")) +
          labs(y="Study\n")}
    }, height=75) 
    
    #PULLDOWN chart title
    output$FCplot_PULLDOWN_title <- renderText({
      if (nrow(datanoNA_PULLDOWN)<1) {paste("<b>(D) Enrichment of", "<i>", gene, "</i>mRNA in miR-34a pulldown studies:", gene,"mRNA was not detected in any study")}
      else {paste("<b>(D) Enrichment of", "<i>", gene,"</i>mRNA in miR-34a pulldown studies")}})
    
    #barplot of Enrichment miR-34a PULLDOWN  
    output$FCPlot_PULLDOWN <- renderPlot({
      if (nrow(datanoNA_PULLDOWN)<1) {}
      else {
        ggplot(datanoNA_PULLDOWN, aes(x=miR34aFC_PULLDOWN, y=Name_PULLDOWN)) + 
          geom_bar(stat = "identity", color="black", fill="darkblue", width=0.75) +
          scale_x_continuous(trans='log2') +
          theme(text=element_text(size=font_size)) +
          theme(axis.text.y=element_text(size=font_size))+
          theme(axis.title.x=element_text(size=font_size))+
          geom_vline(xintercept = 1, color = "black") +
          theme(panel.background = element_rect(fill = "white", color="black")) +
          labs(x=paste (gene, "mRNA enrichment (fold change miR34a/control pulldown)")) +
          labs(y="Study\n")}
    }, height=75)     
    
    #plot of miRNA prediction heatmap   
    output$predictPlot <- renderPlot({
      predictName <- datasetmiR34atarget$Name
      miR34a <- datasetmiR34atarget[[gene]]
      miR34b <- datasetmiR34btarget[[gene]]
      miR34c <- datasetmiR34ctarget[[gene]]
      predictHeatMap <- data.frame(predictName, miR34a, miR34b, miR34c)
      melted_predictHeatMap <- reshape2::melt(predictHeatMap)          
      melted_predictHeatMap$Y1 <- cut(melted_predictHeatMap$value,breaks = c(-Inf,0:3,Inf),right = FALSE)
      
      ggplot(melted_predictHeatMap, aes(x = variable, y = predictName, fill = value)) +
        geom_tile(aes(fill = Y1), colour = "gray90", lwd=1) +
        scale_fill_manual(breaks=c("[0,1)", "[1,2)", "[2,3)", "[3, Inf)"),
                          values = c("white", "darkblue", "darkgreen", "firebrick"))+  
        theme(text=element_text(size=font_size)) +
        theme(axis.text.y=element_text(size=font_size))+
        theme(axis.text.x=element_text(size=font_size))+
        theme(panel.background = element_rect(fill = "white", color="black")) +
        labs(x=NULL) +
        labs(y="mRNA target prediction algorithm\n")+
        theme(legend.position = "none")+
        scale_x_discrete(position = "top") 
      
      
    },height=170)
    
    #Validated targets; uses trasposed miR-34a vlidated targets table (non-validated need "none")
    output$validated <- renderText({
      if (miR34avalidate == "none") 
        {paste("<b>", gene, "is not a published miR-34a target")}
      else {paste("<div style='background-color:green'><font color=white><background=green><b>", gene, "is a published miR-34a target: ", miR34avalidate)}
    })
    
    #depmap link
    output$depmap <- renderUI(a(href=paste("https://depmap.org/portal/gene/", gene, sep=""), "DepMap portal", target="_blank"))
    
    #show selected and hide unselected summary and raw data panels
    observeEvent(input$FCplot_click, {
      shinyjs::hide(id = "summaryFCKO")
      shinyjs::hide(id = "RawDataFCKO")
      #shinyjs::hide(id = "CorrelationPlot34a")   
      #shinyjs::hide(id = "CorrelationPlot34b")
      #shinyjs::hide(id = "CorrelationPlot34c")
      shinyjs::show(id = "summaryFC")
      shinyjs::show(id = "RawDataFC")
    })
    
    
    observeEvent(input$FCplotKO_click, {
      shinyjs::hide(id = "summaryFC")
      shinyjs::hide(id = "RawDataFC")
      #shinyjs::hide(id = "CorrelationPlot34a")   
      #shinyjs::hide(id = "CorrelationPlot34b")
      #shinyjs::hide(id = "CorrelationPlot34c")
      shinyjs::show(id = "summaryFCKO")   
      shinyjs::show(id = "RawDataFCKO")
    })
    
    observeEvent(input$Correlation_click, {
      shinyjs::hide(id = "summaryFC")
      shinyjs::hide(id = "RawDataFC")
      shinyjs::hide(id = "summaryFCKO")
      shinyjs::hide(id = "RawDataFCKO")
      #shinyjs::show(id = "CorrelationPlot34a")   
      #shinyjs::show(id = "CorrelationPlot34b")
      #shinyjs::show(id = "CorrelationPlot34c")
    })
    
    #barplot/dotplot of mRNA expression in miR34a OE single studies
    output$RawDataFC <- renderPlot({
      if (is.null(input$FCplot_click$y)) return("")
      else {
        #define dataset after click on the plot
        dataset <- get(eval(datanoNA$ID[round(input$FCplot_click$y)]))
        genedata <- dataset[[gene]]
        Namedataset <- dataset$Name
        V1dataset <- dataset$order
        
        #calculate mean and SD
        data <- data.frame(V1dataset, Namedataset, genedata)
        data.summary <- data %>%
          group_by(Namedataset) %>%
          summarise(
            sd = sd(genedata, na.rm = TRUE),
            genedata = mean(genedata),
            V1dataset = mean(V1dataset)
          )
        
        #make plot
        myorder <- data.summary$Namedataset[order(data.summary$V1dataset)] #order (sort) groups according to V1 column
        ggplot(data, aes(Namedataset, genedata)) +
          geom_bar(stat="identity", data = data.summary, fill = "gray75", color = "black") +
          scale_x_discrete(limits=myorder)+
          theme(text=element_text(size=font_size)) +
          theme(axis.text.y=element_text(size=font_size))+
          theme(axis.text.x=element_text(size=font_size))+
          theme(axis.text.x=element_text(angle=45, hjust=1))+
          theme(panel.background = element_rect(fill = "white", color="black")) +
          labs(x=element_blank()) +
          labs(y=paste (gene, "mRNA expression\n"))+
          geom_dotplot(binwidth=max(genedata)/30, binaxis="y", stackdir="center", fill ="black")+ #binwidth should not be a constant number, otherwise it will change when expression values (y-axis) change
          geom_errorbar( aes(ymin = genedata-sd, ymax = genedata+sd), 
                         data = data.summary, width = 0.2) 
        
      }
    }, height=250)

    #barplot/dotplot of mRNA expression in miR34a KO mice single studies
    output$RawDataFCKO <- renderPlot({
      if (is.null(input$FCplotKO_click$y)) return("")
      else {
        #define dataset after click on the plot
        datasetKO <- get(eval(dataKOnoNA$ID[round(input$FCplotKO_click$y)]))
        genedataKO <- datasetKO[[gene]]
        NamedatasetKO <- datasetKO$Name
        V1datasetKO <- datasetKO$order
        
        #calculate mean and SD
        dataKO <- data.frame(V1datasetKO, NamedatasetKO, genedataKO)
        data.summaryKO <- dataKO %>%
          group_by(NamedatasetKO) %>%
          summarise(
            sdKO = sd(genedataKO, na.rm = TRUE),
            genedataKO = mean(genedataKO),
            V1datasetKO = mean(V1datasetKO)
          )
        
        #make plot
        myorder <- data.summaryKO$NamedatasetKO[order(data.summaryKO$V1datasetKO)] #order (sort) groups according to V1 column
        ggplot(dataKO, aes(NamedatasetKO, genedataKO)) +
          geom_bar(stat="identity", data = data.summaryKO, fill = "gray75", color = "black") +
          scale_x_discrete(limits=myorder)+
          theme(text=element_text(size=font_size)) +
          theme(axis.text.y=element_text(size=font_size))+
          theme(axis.text.x=element_text(size=font_size))+
          theme(axis.text.x=element_text(angle=45, hjust=1))+
          theme(panel.background = element_rect(fill = "white", color="black")) +
          labs(x=element_blank()) +
          labs(y=paste (gene, "mRNA expression\n"))+
          geom_dotplot(binwidth=max(genedataKO)/30, binaxis="y", stackdir="center", fill ="black")+ #binwidth should not be a constant number, otherwise it will change when expression values (y-axis) change
          geom_errorbar( aes(ymin = genedataKO-sdKO, ymax = genedataKO+sdKO), 
                         data = data.summaryKO, width = 0.2)
      }
    }, height=250)
    
    
    
    output$summaryFC <- renderTable({ #GEO miR34a induction studies summary table
      if (is.null(input$FCplot_click$y)) return()
      else {
        datasetGSE <- datanoNA$ID[round(input$FCplot_click$y)]
        summaryGSE <- summarymiR34a[[datasetGSE]]
        summarytable <- data.frame(summaryfirstrow, summaryGSE)
        names(summarytable) <- NULL
        summarytable
      }
    }) 
    
    output$summaryFCKO <- renderTable({ #GEO miR34a KO mice studies summary table
      if (is.null(input$FCplotKO_click$y)) return()
      else {
        datasetGSEKO <- dataKOnoNA$ID[round(input$FCplotKO_click$y)]
        summaryGSEKO <- summarymiR34aKO[[datasetGSEKO]]
        summarytableKO <- data.frame(summaryfirstrowKO, summaryGSEKO)
        names(summarytableKO) <- NULL
        summarytableKO
      }
    }) 
    
    
    'output$CorrelationPlot34a <- renderPlot({ #Correlation plot miR34a
      if (is.null(input$Correlation_click$y)) return()
      else {
        datasetCORRELATION <- get(eval(datasetmiR34a_correlation$ID[round(input$Correlation_click$y*3)]))
        ggscatter(datasetCORRELATION, x = "mir34a", y = gene, 
                  color = "darkblue",
                  add = "reg.line", conf.int = TRUE, 
                  add.params = list(color = "black",
                                    fill = "lightgray"),
                  title = "correlation with miR-34a",
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "miR-34a expression (log2 RSEM)", ylab = paste(gene, "expression (log2 RSEM)"))
      }
    }) 
    
    output$CorrelationPlot34b <- renderPlot({ #Correlation plot miR34b
      if (is.null(input$Correlation_click$y)) return()
      else {
        datasetCORRELATION <- get(eval(datasetmiR34a_correlation$ID[round(input$Correlation_click$y*3)]))
        ggscatter(datasetCORRELATION, x = "mir34b", y = gene, 
                  color = "darkgreen",
                  add = "reg.line", conf.int = TRUE, 
                  add.params = list(color = "black",
                                    fill = "lightgray"),
                  title = "correlation with miR-34b",
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "miR-34b expression (log2 RSEM)", ylab = paste(gene, "expression (log2 RSEM)"))
      }
    }) 
    
    output$CorrelationPlot34c <- renderPlot({ #Correlation plot miR34c
      if (is.null(input$Correlation_click$y)) return()
      else {
        datasetCORRELATION <- get(eval(datasetmiR34a_correlation$ID[round(input$Correlation_click$y*3)]))
        ggscatter(datasetCORRELATION, x = "mir34c", y = gene, 
                  color = "firebrick",
                  add = "reg.line", conf.int = TRUE, 
                  add.params = list(color = "black",
                                    fill = "lightgray"),
                  title = "correlation with miR-34c",
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "miR-34c expression (log2 RSEM)", ylab = paste(gene, "expression (log2 RSEM)"))
      }
    })' 
    
    
# Download FC data
'   output$downloadData <- downloadHandler(
      filename = function() {
        paste0("Regulation of ", gene, " by miR34a.txt", sep = "")
      },
      content = function(file) {
        vroom::vroom_write(do.call(rbind, Map(data.frame, Study=miR34adata$Name, 
                                              "Fold Change"=miR34aFC)), file)
      })'
  })
  
 
#Screen server
  
  INTERSECTION <- reactive({
    if (input$GEOOEcheck==TRUE) GEOOElessgenes <- subset(precalcGEOOE, get(eval(paste0(paste("less", input$GEOOEsliderFC,sep=""))))>=input$GEOOEsliderNO)$Name
    else GEOOElessgenes <- precalcGEOOE$Name
    if (input$GEOKOcheck==TRUE) GEOKOmoregenes <- subset(precalcGEOKO, get(eval(paste0(paste("more",input$GEOKOsliderFC,sep=""))))>=input$GEOKOsliderNO) $Name   
    else GEOKOmoregenes <- precalcGEOOE$Name
    if (input$predictioncheckA==TRUE) PREDICTAmoregenes <- subset(precalcPREDICT, miR34a>=input$predictAslider)$Name
    else PREDICTAmoregenes <- precalcGEOOE$Name
    if (input$predictioncheckB==TRUE) PREDICTBmoregenes <- subset(precalcPREDICT, miR34b>=input$predictBslider)$Name
    else PREDICTBmoregenes <- precalcGEOOE$Name
    if (input$predictioncheckC==TRUE) PREDICTCmoregenes <- subset(precalcPREDICT, miR34c>=input$predictCslider)$Name
    else PREDICTCmoregenes <- precalcGEOOE$Name
    if (input$correlAcheck==TRUE) correlAlessgenes <- subset(precalcCORRELa, get(eval(paste0(paste("less",input$correlAsliderR,sep=""))))>=input$correlAsliderNO)$Name
    else correlAlessgenes <- precalcGEOOE$Name
    if (input$correlBcheck==TRUE) correlBlessgenes <- subset(precalcCORRELb, get(eval(paste0(paste("less",input$correlBsliderR,sep=""))))>=input$correlBsliderNO)$Name
    else correlBlessgenes <- precalcGEOOE$Name
    if (input$correlCcheck==TRUE) correlClessgenes <- subset(precalcCORRELc, get(eval(paste0(paste("less",input$correlCsliderR,sep=""))))>=input$correlCsliderNO)$Name
    else correlClessgenes <- precalcGEOOE$Name
    if (input$SILACcheck==TRUE) SILAClessgenes = subset(precalcSILAC, get(eval(paste0(paste("less",input$SILACsliderFC,sep=""))))>=input$SILACsliderNO)$Name
    else SILAClessgenes <- precalcGEOOE$Name
    if (input$PULLDOWNcheck==TRUE) PULLDOWNmoregenes = subset(precalcPULLDOWN, get(eval(paste0(paste("more",input$PULLDOWNFC,sep=""))))>=input$PULLDOWNsliderNO)$Name
    else PULLDOWNmoregenes <- precalcGEOOE$Name
    INTERSECTION <- as.data.frame(Reduce(intersect, list(GEOOElessgenes, GEOKOmoregenes, PREDICTAmoregenes, PREDICTBmoregenes, PREDICTCmoregenes, correlAlessgenes, correlBlessgenes, correlClessgenes, SILAClessgenes, PULLDOWNmoregenes)))
  })
  
  colnames(datasetmiR34avalid) <- c('Reduce(intersect, list(GEOOElessgenes, GEOKOmoregenes, PREDICTAmoregenes, PREDICTBmoregenes, PREDICTCmoregenes, correlAlessgenes, correlBlessgenes, correlClessgenes, SILAClessgenes, PULLDOWNmoregenes))', 'ref')

  output$INTERSECTNOgenes <- renderText({
    if (nrow(INTERSECTION()) == 24396) paste("<b>Select criteria and cutoff values for the identification of potential miR-34 targets")
    else paste("<b>", nrow(INTERSECTION()), "mRNAs fulfill the selected criteria. Click on a mRNA to get more information.")
  })
  
  output$INTERSECTresults <- DT::renderDataTable({ #uses non-trasposed miR-34 avlidated target table; non-validated ref is blank
    if (nrow(INTERSECTION()) == 24396) mergedtable = NULL
    else {
      mergedtable <- merge(INTERSECTION(), datasetmiR34avalid, by="Reduce(intersect, list(GEOOElessgenes, GEOKOmoregenes, PREDICTAmoregenes, PREDICTBmoregenes, PREDICTCmoregenes, correlAlessgenes, correlBlessgenes, correlClessgenes, SILAClessgenes, PULLDOWNmoregenes))")
      colnames(mergedtable) <- c('mRNA', 'published miR-34 target')
    }
    datatable( data = mergedtable #https://rstudio.github.io/DT/extensions.html
               , extensions = 'Buttons',
               selection = 'single'
               , options = list(
                 paging = FALSE #no pagination; http://www.baoruidata.com/examples/018-datatable-options/
                 , dom = "Blfrtip" #https://datatables.net/reference/option/dom
                 , buttons = list("copy", "csv") # only two buttons, copy, and download csv; https://rstudio.github.io/DT/003-tabletools-buttons.html
              ))
  })
  
  #store gene symbol from the selected row as selRow value; https://community.rstudio.com/t/getting-the-user-selected-entry-of-a-data-table-in-shiny-to-actually-make-a-scatterplot/52419
  observe({
    req(input$INTERSECTresults_rows_selected)
    selRow <- INTERSECTION()[input$INTERSECTresults_rows_selected,] #print(selRow[[1]]) - shows selRow variable value in console
    print(selRow)
    gene1 <- selRow
    updateTabsetPanel(session, "Tabset", selected = "searchgenes") #https://stackoverflow.com/questions/43552906/how-to-switch-between-navbar-tabs-with-a-button-r-shiny
    updateTextInput(session, "gene1", value = selRow) #https://shiny.rstudio.com/reference/shiny/0.14/updateTextInput.html
    click("go")
    proxy = dataTableProxy('INTERSECTresults') #unselect row in the INTERSECTresult table (otherwise the screener does not work well), https://yihui.shinyapps.io/DT-proxy/
    proxy %>% selectRows(NULL)
  }) 
  
  #when the reset button in the screener is pressed
  observeEvent(input$reset, { 
    updateCheckboxInput(session =  session, inputId = "GEOOEcheck", value = FALSE) #https://stackoverflow.com/questions/56820844/how-can-i-change-r-shiny-checkboxinput-value-to-false-true-programmatically
    updateCheckboxInput(session =  session, inputId = "GEOKOcheck", value = FALSE)    
    updateCheckboxInput(session =  session, inputId = "predictioncheckA", value = FALSE)
    updateCheckboxInput(session =  session, inputId = "predictioncheckB", value = FALSE)
    updateCheckboxInput(session =  session, inputId = "predictioncheckC", value = FALSE)
    updateCheckboxInput(session =  session, inputId = "correlAcheck", value = FALSE)
    updateCheckboxInput(session =  session, inputId = "correlBcheck", value = FALSE)
    updateCheckboxInput(session =  session, inputId = "correlCcheck", value = FALSE)
    updateCheckboxInput(session =  session, inputId = "SILACcheck", value = FALSE)
    updateCheckboxInput(session =  session, inputId = "PULLDOWNcheck", value = FALSE)
    
    updateSliderInput(session, "GEOOEsliderFC", value = 0.8) #https://shiny.rstudio.com/reference/shiny/1.4.0/updateSliderInput.html
    updateSliderInput(session, "GEOOEsliderNO", value = 10)
    updateSliderInput(session, "GEOKOsliderFC", value = 1.25)
    updateSliderInput(session, "GEOKOsliderNO", value = 3)
    updateSliderInput(session, "predictAslider", value = 6)
    updateSliderInput(session, "predictBslider", value = 6)
    updateSliderInput(session, "predictCslider", value = 6)
    updateSliderInput(session, "correlAsliderR", value = -0.1)
    updateSliderInput(session, "correlBsliderR", value = -0.1)
    updateSliderInput(session, "correlCsliderR", value = -0.1)
    updateSliderInput(session, "correlAsliderNO", value = 4)
    updateSliderInput(session, "correlBsliderNO", value = 4)
    updateSliderInput(session, "correlCsliderNO", value = 4)
    updateSliderInput(session, "SILACsliderFC", value = 0.8)
    updateSliderInput(session, "SILACsliderNO", value = 1)
    updateSliderInput(session, "PULLDOWNFC", value = 2)
    updateSliderInput(session, "PULLDOWNsliderNO", value = 1)
  })
  
  
  'msigdbr_H <- msigdbr(species = "human", category = "H") #Hallmark
  msigdbr_H1 = msigdbr_H %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame() ## fixing format to work with enrichr; https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
  msigdbr_KEGG <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG") #KEGG pathways
  msigdbr_KEGG1 = msigdbr_KEGG %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
  msigdbr_GOBP <- msigdbr(species = "human", category = "C2", subcategory = "CP:WIKIPATHWAYS") #GO biological process
  msigdbr_GOBP1 = msigdbr_GOBP %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
  
  
  output$GSEA_H <- renderPlot({
    if (nrow(INTERSECTION()) == 24396) NULL
    else {
    enrich_H <- enricher(gene = INTERSECTION()[[1]], TERM2GENE = msigdbr_H1, pAdjustMethod="none", pvalueCutoff = 1, qvalueCutoff = 1)  #https://www.rdocumentation.org/packages/clusterProfiler
    enrich_Hpvalue <- enrich_H$pvalue
    enrich_Hpvallog <- -log10(enrich_Hpvalue)
    enrich_H@result$pvaluelog <- enrich_Hpvallog
    barplot(enrich_H, x="pvaluelog", title = "Hallmark") #https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/dotplot,compareClusterResult-method
    }
  })
  
  output$GSEA_KEGG <- renderPlot({
    if (nrow(INTERSECTION()) == 24396) NULL
    else {
    enrich_KEGG <- enricher(gene = INTERSECTION()[[1]], TERM2GENE = msigdbr_KEGG1, pAdjustMethod="none", pvalueCutoff = 1, qvalueCutoff = 1)
    enrich_KEGGpvalue <- enrich_KEGG$pvalue
    enrich_KEGGpvallog <- -log10(enrich_KEGGpvalue)
    enrich_KEGG@result$pvaluelog <- enrich_KEGGpvallog
    barplot(enrich_KEGG, x="pvaluelog", title = "KEGG") 
    }
  })
  
  output$GSEA_GOBP <- renderPlot({
    if (nrow(INTERSECTION()) == 24396) NULL
    else {
      enrich_GOBP <- enricher(gene = INTERSECTION()[[1]], TERM2GENE = msigdbr_GOBP1, pAdjustMethod="none", pvalueCutoff = 1, qvalueCutoff = 1)
      enrich_GOBPpvalue <- enrich_GOBP$pvalue
      enrich_GOBPpvallog <- -log10(enrich_GOBPpvalue)
      enrich_GOBP@result$pvaluelog <- enrich_GOBPpvallog
      barplot(enrich_GOBP, x="pvaluelog", title = "WIKI Pathways") 
    }
  })'
  
  output$DLmiR34OE <- downloadHandler(
    filename <- "miR34_induction_FC.xlsx",
    content <- function(file) {file.copy("www/miR34_induction_FC.xlsx", file)},
    contentType = "application/zip")
  
  output$DLmiR34KO <- downloadHandler(
    filename <- "miR34_mouse_KO.xlsx",
    content <- function(file) {file.copy("www/miR34_mouse_KO.xlsx", file)},
    contentType = "application/zip")
  
  output$DLmiR34pulldown <- downloadHandler(
    filename <- "miR34a_PULLDOWN.xlsx",
    content <- function(file) {file.copy("www/miR34a_PULLDOWN.xlsx", file)},
    contentType = "application/zip")
  
  output$DLmiR34silac <- downloadHandler(
    filename <- "miR34a_SILAC.xlsx",
    content <- function(file) {file.copy("www/miR34a_SILAC.xlsx", file)},
    contentType = "application/zip")
  
  output$DLmiR34prediction <- downloadHandler(
    filename <- "miR34_target_prediction.xlsx",
    content <- function(file) {file.copy("www/miR34_target_prediction.xlsx", file)},
    contentType = "application/zip")
  
  output$DLmiR34correl <- downloadHandler(
    filename <- "miR34_correlation.xlsx",
    content <- function(file) {file.copy("www/miR34_correlation.xlsx", file)},
    contentType = "application/zip")
  
  output$DLmiR34validate <- downloadHandler(
    filename <- "published miR34a target genes (updated January 2022).xlsx",
    content <- function(file) {file.copy("www/published miR34a target genes (updated January 2022).xlsx", file)},
    contentType = "application/zip")
  
  
  output$about <- renderText("To reference information from this database, please cite the following paper:")
  
  output$contact1 <- renderText("<b>Contact:")
  output$contact2 <- renderText("Matjaz Rokavec")
  output$contact3 <- renderUI(a(href="https://www.pathologie.med.uni-muenchen.de/020wissenschaft/009ag_hermeking/engl_ag_hermeking/index.html", "Heiko Hermeking lab", target="_blank"))  
  output$contact4 <- renderText("Institute of Pathology")
  output$contact5 <- renderText("Ludwig Maximilians University, Munich Germany")
  output$contact6 <- renderText("Email: matjaz.rokavec@med.uni-muenchen.de")

  
  
  #output$manual <- downloadHandler(
  #  filename <- "manual.pdf",
  #  content <- function(file) {file.copy("www/manual.pdf", file)},
  #  contentType = "application/zip")
  
}

shinyApp(ui, server)