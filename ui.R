
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyWidgets)
library(reactable)

shinyUI(navbarPage(title = "KRSA",
                   header = tagList(
                     useShinydashboard(),
                     useShinyjs(),
                     tags$style(HTML("
                        .shiny-output-error-validation {
                         color: #ff0000;
                         font-weight: bold;
                         font-size: 20px;
                           }
                          "))
                   ),
                   id = "tabs",
                   theme = "style/style.css",
                   footer = includeHTML("footer.html"),
                   fluid = TRUE, 
                   collapsible = TRUE,
                   
                    

                   
                   
                   # Home Panel -----
                   tabPanel("Home",
                            includeHTML("home.html"),
                            tags$script(src = "plugins/scripts.js"),
                            tags$head(
                              tags$link(rel = "stylesheet", 
                                        type = "text/css", 
                                        href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css"),
                              tags$link(rel = "icon", 
                                        type = "image/png", 
                                        href = "images/logo_icon.png")
                            )
                   ),
                   
                   # Step1: Input Panel -----
                   tabPanel("Step1: Input",

                            fluidRow(
                              shinydashboard::box(title = "Signal",width = 6, status = "info",solidHeader=TRUE,
                                     fileInput("input_file", label = "Upload the Median_SigmBg crosstab file: ",
                                               accept = c("txt")
                                     ),
                                     actionButton("load_ex_data", "Use Example Dataset")
                              ),
                              # box(title = "Optional: Signal Saturation",width = 4,collapsible = T,collapsed = T,
                              #        fileInput("input_file2", label = "Optional: Upload the signal saturation crosstab file: ",
                              #                  accept = c("txt")
                              #        ),
                              #        actionButton("load_ex_data2", "Use Example Dataset")
                              # ),
                              shinydashboard::box(title = "Kinase Mapping", width = 6,status = "primary",solidHeader=TRUE,
                                     fileInput("map_file", label = "Upload a kinase-substrate association file: ",
                                               accept = c("txt")
                                     ),
                                     actionButton("load_map_file", "Use KRSA mapping")
                              )
                            ),
                            
                            fluidRow(
                              shinydashboard::box(title = "Preview Signal Table",width = 6, status = "info",solidHeader=TRUE,
                                  tableOutput("sig_tbl_preview")
                              ),
                              shinydashboard::box(title = "Preview Kinase Mapping", width = 6,status = "primary",solidHeader=TRUE,
                                  tableOutput("map_tbl_preview")
                              )
                            ),
                            
                            fluidRow(
                              column(12,align="center",
                                  actionButton("input_step_btn", "Go to Step 2", icon("paper-plane"), 
                                               style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                  tags$hr(),span(textOutput("err"), style="color:red")
                                  )
                            ),
                            
                            
                            
                   ),
                   
                   # Step2: Design Panel -----
                   tabPanel("Step2: Design Options",
                            
                            fluidRow(
                              shinydashboard::box(title = "Design",width = 6, status = "info",solidHeader=TRUE,
                                                  selectInput("group_col", label = "Select Column to define groups: ",
                                                            choices = list(""), selected = NULL
                                                  ),
                                                  selectInput("ctl_group", label = "Select the Control Group:",
                                                              choices = list("")
                                                  ),
                                                  selectInput("case_group", label = "Select the Case Group:",
                                                              choices = list("")
                                                  ),
                                                  selectInput("sampleName_col", label = "Select Columns to define uniques samples: ",
                                                              choices = list(""), selected = NULL, multiple = T
                                                  ),
                                                  
                              ),
                              
                              shinydashboard::box(title = "QC Options", width = 6,status = "primary",solidHeader=TRUE,
                                                  sliderInput("max_sig_qc", "Max Exposure Signal", min = 1, max = 100, value = 5, step = 1),
                                                  sliderInput("r2_qc", "Min R2", min = 0, max = 0.99, value = 0.90, step = 0.05)
                              )

                            ),
                            
                            fluidRow(
                              shinydashboard::box(title = "DP Options",width = 6, status = "info",solidHeader=TRUE,
                                                  sliderInput("lfc_thr", "Log2 Fold Change Cutoff", min = 0, max = 5, value = 0.2, step = 0.05),
                              ),
                              shinydashboard::box(title = "Sampling Options", width = 6,status = "primary",solidHeader=TRUE,
                                                  sliderInput("itr_num", "Number of Iterations", min = 100, max = 2000, value = 500, step = 100),
                                                  switchInput("use_seed", "Use Seed?", value = FALSE, onLabel = "Yes", offLabel = "No"),
                                                  numericInput("use_seed_num", "Input Seed Number: ", value = 123)
                              )
                            ),
                            
                            
                            fluidRow(
                              column(12,align="center",
                                     actionButton("start_krsa", "Run KRSA", icon("paper-plane"), 
                                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                     tags$hr(),span(textOutput("err2"), style="color:red")
                              )
                            )   
                            
                   ),
                   
                   # Step3: Results Overview Panel -----
                   tabPanel("Results: Overview", 
                            
                            tabsetPanel(#"output_panels",
                              tabPanel("Summary", 
                                       fluidRow(
                                         infoBoxOutput("summary_options")
                                       ),
                                       
                                       fluidRow(
                                         shinydashboard::box(title = "Peptides Selection",width = 12, status = "info",solidHeader=TRUE,
                                        valueBoxOutput("init_peps", width = 3),
                                        valueBoxOutput("qc_maxSig_peps", width = 3),
                                        valueBoxOutput("qc_r2_peps", width = 3),
                                        valueBoxOutput("lfc_peps", width = 3),
                                       )),
                                       
                                       fluidRow(
                                         shinydashboard::box(title = "LFC Table",width = 6, status = "info",solidHeader=TRUE,
                                                             dataTableOutput("lfc_table"),
                                                             downloadButton("lfc_table_download", label = "Save Table")
                                         ),
                                         shinydashboard::box(title = "Model Table", width = 6,status = "primary",solidHeader=TRUE,
                                                             dataTableOutput("model_table"),
                                                             downloadButton("model_table_download", label = "Save Table")
                                         )
                                       )
                                       
                                       
                                       ),
                              tabPanel("Heatmap", 
                                       selectInput("heatmap_op1", choices = c("Normalized", "Not"), label = "Option"),
                                       plotOutput("heatmap")
                                       ),
                              tabPanel("Boxplot", 
                                       
                                       ),
                              tabPanel("Cuvres Plot", 
                                       
                                       ),
                              tabPanel("WaterFall", 
                                       
                                       )
                            )
                            
                            ),
                   
                   # Kinase Panel -----
                   tabPanel("Results: Kinases", 
                            
                            tabsetPanel(
                              tabPanel("Table",
                                       ),
                              
                              tabPanel("Histogram",
                                       ),
                              tabPanel("Reverse KRSA",
                                       )
                            )
                            ),
                   
                   # Network Panel -----
                   tabPanel("Results: Network", 
                            
                            )
                   
                   
))