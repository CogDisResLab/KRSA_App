library(shiny)
library(shinydashboard)
#library(leaflet)
library(data.table)
library(DT)
#library(plyr)
library(ggplot2)
library(dplyr)
#library(gplots)
#library(igraph)
#library(stringr)
library(shinycssloaders)
library(tidyverse)
library(broom)

source("scripts/reading_data.R")
source("scripts/dropdown_btn.R")
source("scripts/krsa.R")
source("scripts/heatmap.R")
source("scripts/network.R")


###############################################################################
ui <- shinyUI(fluidPage(
    
    includeCSS("styles.css"),
    
    #Application title
    headerPanel("Kinome Random Sampling Analyzer"),
    
    #Layout type
    sidebarLayout(
        
        #Sidebar Panel
        sidebarPanel(
            
            fileInput('file0', label = h3("Upload raw data file"),
                      accept=c('text/csv', 'text/comma-separated-values','text/plain', '.csv')),
            actionButton("upload1", "upload"),
            #Options for file type including header, seperator, and quotes.
            h4("File options:"),
            checkboxInput('header1', 'Header', TRUE),
            radioButtons('sep1', 'Separator', inline=F,
                         c(Tab='\t',
                           Comma=',',
                           Semicolon=';'),
                         '\t'),
            radioButtons('quote1', 'Quote', inline=F,
                         c(None='',
                           'Double Quote'='"',
                           'Single Quote'="'"),
                         '"'),
            tags$hr(),
            #TODO: uncomment below for PINET
            # radioButtons("radio", label = h3("Chose your kinase-peptide association option"),
            #              choices = list("PiNET" = 1, "Upload" = 2), 
            #              selected = 2),
            
            #This takes a tab or comma seperated file as input, with the first column "Substrates" and
            #the second column "Kinases" being the only two columns and both required.
            conditionalPanel(
                #TODO: uncomment below for PINET and comment out "condition = TRUE"
                #condition = "input.radio==2",
                condition = TRUE,
                fileInput('file1', label = h5("Kinase-peptide association file"),
                          accept=c('text/csv', 'text/comma-separated-values','text/plain', '.csv')),
                h4("File options:"),
                checkboxInput('header2', 'Header', TRUE),
                radioButtons('sep2', 'Separator', inline=F,
                             c(Tab='\t',
                               Comma=',',
                               Semicolon=';'),
                             '\t'),
                radioButtons('quote2', 'Quote', inline=F,
                             c(None='',
                               'Double Quote'='"',
                               'Single Quote'="'"),
                             '"')
            ),
            tags$hr()
        ),
        
        #Main Panel
        mainPanel(
            tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: #7c7c7c;  color:#d6d6d6}
    .tabbable > .nav > li[class=active]    > a {background-color: #008E94; color:white}
  ")),
            tabsetPanel(
                tabPanel("Step 1: Options",
                         
                         h2("Data Selection:"),
                         #celltype
                         fluidRow(
                             column(4,
                                    dropdownButton(
                                        label = "Step 1: Experiment Barcode", status = "default", width = 80,
                                        actionButton(inputId = "celltypea2z", label = "Sort A to Z", icon = icon("sort-alpha-asc")),
                                        actionButton(inputId = "celltypez2a", label = "Sort Z to A", icon = icon("sort-alpha-desc")),
                                        br(),
                                        actionButton(inputId = "celltypeall", label = "(Un)select all"),
                                        checkboxGroupInput(inputId = "celltype", label = "Choose", choices = c()))),
                             column(4,
                                    dropdownButton(
                                        label = "Step 2: Control Samples", status = "default", width = 80,
                                        actionButton(inputId = "controlSamplesa2z", label = "Sort A to Z", icon = icon("sort-alpha-asc")),
                                        actionButton(inputId = "controlSamplesz2a", label = "Sort Z to A", icon = icon("sort-alpha-desc")),
                                        br(),
                                        actionButton(inputId = "controlSamplesall", label = "(Un)select all"),
                                        checkboxGroupInput(inputId = "controlSamples", label = "Choose", choices = c()))),
                             column(4,
                                    dropdownButton(
                                        label = "Step 3: Experimental Samples", status = "default", width = 80,
                                        actionButton(inputId = "experimentalSamplesa2z", label = "Sort A to Z", icon = icon("sort-alpha-asc")),
                                        actionButton(inputId = "experimentalSamplesz2a", label = "Sort Z to A", icon = icon("sort-alpha-desc")),
                                        br(),
                                        actionButton(inputId = "experimentalSamplesall", label = "(Un)select all"),
                                        checkboxGroupInput(inputId = "experimentalSamples", label = "Choose", choices = c())))
                         ),
                         tags$hr(),
                         h2("Preview Data:"),
                         
                         fluidRow(
                             h3("Preview of all data"),
                             shiny::dataTableOutput(outputId="contents0"),
                             #DT::dataTableOutput('contents0'),
                             #column(6,
                             #DT::dataTableOutput('contents1'),
                             h3("Preview of control data"),
                             shiny::dataTableOutput(outputId="contents1"),
                             #),
                             #column(6,
                             h3("Preview of experimental data"),
                             shiny::dataTableOutput(outputId="contents2")
                             #)
                         )
                         
                ),
                tabPanel("Step 2: Stringency",
                         sliderInput("FCslider", label = h3("Select Fold Change"), min = 0,  max = 2, value = c(0.85,1.15), step=0.05),
                         
                         h3("Select Stringency Level"),
                         fluidRow(
                             column(6,
                                    checkboxGroupInput("StrinRadio",label=NULL, choices = list("Max exposure >= 2" = 1, "R^2 >=0.9" = 2), selected = c(1,2)),
                                    actionButton("calculateData", "Find Significant Substrates")
                             ),
                             column(6,
                                    textOutput("qc_text"),
                                    textOutput("remain2")
                             )
                         ),
                         hr(),
                         shiny::dataTableOutput("fit_table"),
                         downloadButton('downloadDataFC', 'Download Substrate Results'),
                         downloadButton('downloadDataFC2', 'Download Substrate Results 2'),
                         selectInput("pep_sel", "Choose a Substrate:", choices=c()), 
                         plotOutput("pep_plot", width = "100%", height="600px"),
                         downloadButton('downloadData_LR', 'Download Graph'),
                         hr()
                         #dataTableOutput("subTable")
                         
                ),
                tabPanel("Step 3: Iteration", 
                         #This is the slider that allows you to choose the number of iterations of sampling you
                         #want to do, from 0 to 5000 (1000 default) and jumping in increments of 100.
                         sliderInput("itr_num", label = h4("Select the number of iterations"), min = 0, 
                                     max = 5000, value = 1000, step=100, width='100%'),
                         tags$hr(),
                         h4("Confirm and Run"),
                         textOutput("subnum"),
                         textOutput("subnum4"),
                         textOutput("subnum3"),
                         htmlOutput("subnum2"),
                         actionButton("krsa_btn", "Run KRSA")
                         
                ),
                tabPanel("Step 4: Results",
                         shiny::dataTableOutput("overallTable"),
                         downloadButton('downloadDataR', 'Download Results')
                ),
                tabPanel("Step 5: Histograms",
                         selectInput("kin_sel", "Choose a Kinase:", choices=c()), 
                         plotOutput("kin_plot", width = "100%", height="600px"),
                         downloadButton('downloadData', 'Download Histogram')
                ),
                tabPanel("Step 6: Heatmap",
                         fluidRow(
                             column(6,
                                    sliderInput("slider2", label = "Breaks", min = -10,
                                                max = 10, value = c(-4, 4)),
                                    selectInput("clust", "Clustering method:", choices=c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
                                    
                             ),
                             column(6,
                                    selectInput("low", "Choose low phosphorylation color:", choices=c("cyan", "orange", "yellow", "green", "darkgreen", "red", "blue", "purple", "magenta")),
                                    selectInput("mid", "Choose middle phosphorylation color:", choices=c("black", "white")),
                                    selectInput("high", "Choose high phosphorylation color:", choices=c("yellow", "orange", "red", "green", "darkgreen", "cyan", "blue", "purple", "magenta")),
                                    downloadButton('downloadDataH', 'Download Heatmap')
                             )
                         ),
                         
                         hr(),
                         plotOutput("heatmap", width = "100%", height="1200px")
                ),
                tabPanel("Step 7: Network",
                         fluidRow(
                             column(6,
                                    #sliderInput("sliderExpans", label = "Generations", min = 1, max = 3, value=1),
                                    sliderInput("nodeSize", label = "Node Size", min=1, max=10, value=3, step=1),
                                    sliderInput("nodeTextSize", label = "Node Text Size", min=1, max=10, value=6, step=1)
                                    #checkboxInput("low_conf", "Include Low Confidence Connections", FALSE),
                                    #selectInput("interType", "Interaction Type", choices=c("All Interactions", "Kinase Only"))
                             ),
                             column(6,
                                    selectInput("hitsColor", "Choose color for KRSA kinases:", choices=c("red","cyan", "orange", "yellow", "green", "darkgreen", "blue", "purple", "magenta", "white", "gray")),
                                    selectInput("interColor", "Choose color for interacting proteins:", choices=c("gray","red", "orange", "yellow", "green", "darkgreen", "cyan", "blue", "purple", "magenta", "white")),
                                    selectInput("layout", "Choose layout:", choices=c("Circle", "Fruchterman-Reingold", "Kamada Kawai", "LGL")),
                                    downloadButton('downloadDataN', 'Download Network')
                             )
                         ),
                         hr(),
                         withSpinner(plotOutput("networkPlot", width = "100%", height="1200px"), proxy.height = "10px") #TODO: trying to use proxy height to move the spinner up but no luck
                         
                         
                )
            )
            
        )
    )
))


###############################################################################
server <- function(input, output, session) {
    
    session$onSessionEnded(stopApp)
  
    theme_set(theme_bw())
    
    shared_values <- reactiveValues(
      sig_pep = "",
      kin_hits = "",
      ctl_samples = "",
      case_sample = "",
      barcodes = ""
    )
    
    
    
    #Reactive dataset for the input file 0 (raw data).
    data_set0 <- reactive({
        req(input$file0)
        inFile <- input$file0
        x<-read.csv(inFile$datapath, header=F, 
                            sep=input$sep1, quote=input$quote1, stringsAsFactors = F)
        reading_data(x)
         

    })
    
    
      observeEvent(input$upload1,{
        #z <- input$file0
        #data_set0() ->> xy
        data_set0() %>% pull(Barcode) %>% unique() %>% as.list() -> barcode_ids
        #cont <<- pull(data_set0(),Barcode) %>% unique()
        updateCheckboxGroupInput(session, "celltype",
                                 choices = barcode_ids
                                 #selected = NA #This change made it so all of the boxes were selected by default
        )
        
        
        data_set0() %>% pull(`Sample name`) %>% unique() %>% as.list() -> samples_ids
        
        updateCheckboxGroupInput(session, "controlSamples",
                                 choices = samples_ids
                                 #selected = NA #This change made it so all of the boxes were selected by default
        )
        
        updateCheckboxGroupInput(session, "experimentalSamples",
                                 choices = samples_ids
                                 #selected = NA #This change made it so all of the boxes were selected by default
        )
        
      })
      
      
      observeEvent(input$calculateData, {
        
        data_set0() %>% filter(`Sample name` %in% c(input$controlSamples, input$experimentalSamples)) -> filtered_samples
        
        filtered_samples %>% filter(`Exposure time` == 200) %>%
          select(Peptide, `Sample name`, Signal) %>% 
          spread(key = `Sample name`, value = Signal) %>% 
          dplyr::filter_at( vars(-Peptide) , all_vars(. >= as.numeric(2))) %>% 
          dplyr::filter(!Peptide %in% c("pVASP_150_164", "pTY3H_64_78", "ART_025_CXGLRRWSLGGLRRWSL")) %>%
          pull(Peptide) -> ppPassAll
        

        
        filtered_samples %>% dplyr::filter(Peptide %in% ppPassAll, Cycle == 124) %>% 
          select(Peptide, `Sample name`, `Exposure time`, Signal) %>% 
          nest(-`Sample name`, -Peptide) %>% 
          group_by(`Sample name`, Peptide) %>% 
          mutate(
            fit = map(data,~ lm(Signal ~ `Exposure time`, data = .x)),
            coefs = map(fit, tidy),
            summary = map(fit, glance),
            #.default argument to replace NA/NaN values
            slope = map_dbl(coefs,pluck,2,2,.default = 0),
            r.seq = map_dbl(summary, "r.squared")
          ) %>% 
          select(`Sample name`, Peptide, slope, r.seq) -> Dataset_PW_C1_Fit
        
        # Dataset_PW_C1_Fit %>% 
        #   select(`Sample name`, Peptide, r.seq) %>% spread(`Sample name`, r.seq) %>% 
        #   filter_at( vars(-Peptide) , all_vars(. >= 0.9)) %>% pull(Peptide) -> ppPassR2
        
        #pass_both <- intersect(ppPassAll, ppPassR2)
        
        Dataset_PW_C1_Fit %>% filter(Peptide %in% ppPassAll) %>% 
          mutate(Group = ifelse(`Sample name` %in% input$controlSamples, "Control", "Case")) %>% 
          group_by(Group, Peptide) %>% 
          summarise(SignalAvg = mean(slope), r.seqAvg = mean(r.seq)) -> grouped_fit
        
        grouped_fit %>%   
          select(Group, Peptide, r.seqAvg) %>% spread(Group, r.seqAvg) %>% 
          filter_at( vars(-Peptide) , all_vars(. >= 0.9)) %>% pull(Peptide) -> ppPassR2_avg
        
        #setdiff(ppPassR2_avg, ppPassR2)
        
        pass_both_2 <- intersect(ppPassAll, ppPassR2_avg)
        
        grouped_fit %>% filter(Peptide %in% pass_both_2) %>% 
          group_by(Peptide) %>% 
          mutate(FC = SignalAvg[Group == "Case"]/SignalAvg[Group == "Control"]) %>%
          filter(FC <= input$FCslider[1] | FC >= input$FCslider[2]) %>% pull(Peptide) %>% unique() -> sig_pep
        
      
        shared_values$sig_pep <- sig_pep
        
       output$qc_text <-renderText(paste0("removed", as.character(141 - length(ppPassAll)), " and ", as.character(141 - length(ppPassR2_avg))))
      
       
       updateSelectInput(session, "pep_sel", choices = sig_pep)
         
       output$pep_plot <- renderPlot({
         filtered_samples %>% dplyr::filter(Peptide %in% input$pep_sel, Cycle == 124) %>% 
           mutate(Group = ifelse(`Sample name` %in% input$controlSamples, "Control", "Case")) %>% 
           ggplot(aes(`Exposure time`, Signal, group = `Sample name`)) + 
           geom_smooth(method = lm, formula = y~x, se=F, aes(color = factor(Group)))
       })


          output$fit_table <- shiny::renderDataTable({
            Dataset_PW_C1_Fit
            #}, caption = "Preview of raw data file")
          }, options=list(scrollX=TRUE, pageLength = 10))
          
          
          output$heatmap <- renderPlot({
            
            Dataset_PW_C1_Fit %>% select(`Sample name`, Peptide, slope) %>% 
              pivot_wider(names_from = `Sample name`, values_from = slope) %>% 
              column_to_rownames("Peptide") %>% as.matrix() %>% log2() -> mm
            
            plotInputH(mm, input$slider2[1], input$slider2[2], input$low, input$mid, input$high, input$clust)
          })
          
        
      })
      
      output$subnum2 <- renderUI({
        str1 <- paste0("")
        #str2 <- as.character(resultsTable[,1])
        HTML(paste(str1, shared_values$sig_pep, sep = '<br/>'))
      })
      
      
    
    
    
    #Output for rendering the preview table in the main panel.
    #output$contents0 <- DT::renderDataTable({
    output$contents0 <- shiny::renderDataTable({
      data_set0() %>% filter(Barcode %in% input$celltype)
        #}, caption = "Preview of raw data file")
    }, options=list(scrollX=TRUE, pageLength = 10))
    
    output$contents1 <- shiny::renderDataTable({
      data_set0() %>% filter(`Sample name` %in% input$controlSamples)
      #}, caption = "Preview of raw data file")
    }, options=list(scrollX=TRUE, pageLength = 10))
    
    output$contents2 <- shiny::renderDataTable({
      data_set0() %>% filter(`Sample name` %in% input$experimentalSamples) 
      #}, caption = "Preview of raw data file")
    }, options=list(scrollX=TRUE, pageLength = 10))
    
    
    
    observeEvent(input$krsa_btn, {
      req(input$file1)
      
      withProgress(message = 'KRSA', value = 0, {
      
      inFile <- input$file1
      file <- read.csv(inFile$datapath, header=T, 
                  sep=input$sep1, quote=input$quote1, stringsAsFactors = F)
      
      file %>% separate_rows(Kinases) %>% rename(Kin = Kinases) -> cov_file
      
      
      incProgress(0.4)
      krsa_output <- krsa_2(input$itr_num, cov_file, shared_values$sig_pep)
      incProgress(0.8)
      
      output$overallTable <- shiny::renderDataTable({
        krsa_output$krsa_table
      }, options=list(scrollX=TRUE, pageLength = 10))
     
      })
      
      updateSelectInput(session, "kin_sel", choices = krsa_output$krsa_table$Kinase)
      
      output$kin_plot <- renderPlot({
        krsa_output$res %>% rename(Kinase = Kin) %>% 
          filter(Kinase == input$kin_sel) %>% 
          ggplot() +
          geom_histogram(aes(counts),binwidth = 1,fill= "gray30", color = "black") + 
          geom_rect(data=filter(krsa_output$krsa_table, Kinase == input$kin_sel),aes(xmin=SamplingAvg+(2*SD), xmax=SamplingAvg-(2*SD), ymin=0, ymax=Inf),
                    fill="gray", alpha=0.5
          ) +
          geom_vline(data=filter(krsa_output$krsa_table, Kinase == input$kin_sel), aes(xintercept = Observed), color = "red", size = 1,show.legend = F) +
          labs(x = "Hits", y = "Counts", title = input$kin_sel)
      })
      
      output$networkPlot <- renderPlot({
        filter(krsa_output$krsa_table, abs(Z) >= 2) %>% pull(Kinase) -> kin_hits
        plotInputN(kin_hits, input$hitsColor, input$interColor, input$nodeSize, input$nodeTextSize, input$layout)
      })
      
      
      
      
    })
    
    
    
    
    
}


###############################################################################
#Make the app!
shinyApp(ui, server)
