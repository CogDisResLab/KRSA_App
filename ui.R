## Libraries ----
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinydashboard)
library(shinycssloaders)
library(magrittr)
library(bsplus)
library(cicerone)

## Shiny UI ----
shinyUI(navbarPage(title = "KRSA",
                   header = tagList(
                     useShinydashboard(),
                     useShinyjs(),
                     use_cicerone(),
                     use_bs_popover(),
                     use_bs_tooltip(),
                     # GA Tag ------
                     HTML(
                       "<script async src='https://www.googletagmanager.com/gtag/js?id=G-0KRSJDY68Y'></script>
                         <script>
                         window.dataLayer = window.dataLayer || [];
                       function gtag(){dataLayer.push(arguments);}
                       gtag('js', new Date());
                       
                       gtag('config', 'G-0KRSJDY68Y');
                       </script>"
                     ),
                     ## css style ----
                     tags$style(HTML("
                        .shiny-output-error-validation {
                         color: #ff0000;
                         font-weight: bold;
                         font-size: 20px;
                        }
                           
                         .box.box-solid.box-primary>.box-header {
                          color: white;
                          background:#7190a7;
                                             }
                         
                         .box.box-solid.box-primary{
                         border-bottom-color:#bfbfbf;
                         border-left-color:#bfbfbf;
                         border-right-color:#bfbfbf;
                         border-top-color:#bfbfbf;
                         }
                         
                         .box.box-primary>.box-header {
                           color:#000000;
                           background:#fff
                                             }
                         
                         .box.box-primary{
                         border-bottom-color:#bfbfbf;
                         border-left-color:#bfbfbf;
                         border-right-color:#bfbfbf;
                         border-top-color:#bfbfbf;
                         }
                          "))
                   ),
                   ## ui elements -----
                   id = "tabs",
                   theme = "style/style.css",
                   footer = includeHTML("footer.html"),
                   fluid = TRUE, 
                   collapsible = TRUE,

                   # Home Panel -----
                   tabPanel("Home",
                            includeHTML("home.html")
                   ),
                   
                   # Step1: Input Panel -----
                   tabPanel("Step1: Input",
                            
                            fluidRow(column(offset = 11, width = 1, align="center",
                                            actionButton("help1", icon = icon("question"), label = "Help"),
                                            hr()
                                            )),
                            fluidRow(
                              shinydashboard::box(id = "box1_sig",title = "Signal",width = 6, status = "primary",solidHeader=TRUE,
                                     fileInput("input_file", label = "Upload the Median_SigmBg crosstab file: ",
                                               accept = c("txt")
                                     ) %>% 
                                       shinyInput_label_embed(
                                         htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>% 
                                         #shiny_iconlink(id = "hi", href= NULL) %>%
                                           bs_embed_popover(
                                             title = "File Upload",content  = "Upload your Median_SigmBg BioNavigator crosstab view file here. Must be a .txt file and tab delimited", 
                                             placement = "left", trigger = "focus"
                                           )
                                       )
                                     
                                     ,
                                     actionButton("load_ex_data", "Use Example Dataset")
                                     
                              ),
                              # box(title = "Optional: Signal Saturation",width = 4,collapsible = T,collapsed = T,
                              #        fileInput("input_file2", label = "Optional: Upload the signal saturation crosstab file: ",
                              #                  accept = c("txt")
                              #        ),
                              #        actionButton("load_ex_data2", "Use Example Dataset")
                              # ),
                              shinydashboard::box(id = "box2_map",title = "Kinase Mapping", width = 6,status = "primary",solidHeader=TRUE,
                                     fileInput("map_file", label = "Upload a kinase-substrate association file: ",
                                               accept = c("txt")
                                     ) %>% 
                                       shinyInput_label_embed(
                                         htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                           bs_embed_popover(
                                             title = "File Upload", content = "Upload your kinase-substrate association file here. Must be a .txt file and tab delimited. Two columns: Substrates (Peptide IDs on pamchip) and Kinases Kinases separated by spaces", 
                                             placement = "left", trigger = "focus"
                                           )
                                       )
                                     
                                     ,
                                     actionButton("load_map_file", "Use KRSA mapping")
                              ) 
                            ),
                            
                            fluidRow(
                              shinydashboard::box(id = "sig_prev_box", title = "Preview Signal Table",width = 6, status = "primary",solidHeader=TRUE,
                                  tableOutput("sig_tbl_preview")
                              ),
                              shinydashboard::box(id = "map_prev_box", title = "Preview Kinase Mapping", width = 6,status = "primary",solidHeader=TRUE,
                                  tableOutput("map_tbl_preview")
                              )
                            ),
                            
                            fluidRow(
                              column(12,align="center",
                                  actionButton("input_step_btn", "Go to Step 2", icon("paper-plane"), 
                                               style="color: #fff; background-color: #7190a7; border-color: #2e6da4"),
                                  tags$hr(),span(textOutput("err"), style="color:red")
                                  )
                            ),
    
                   ),
                   
                   # Step2: Design Panel -----
                   tabPanel("Step2: Design Options", id = "design_tab",
                            fluidRow(column(offset = 11, width = 1, align="center",
                                            actionButton("help2", icon = icon("question"), label = "Help"),
                                            hr()
                            )),
                            fluidRow(
                              shinydashboard::box(title = "Design", id = "design_box",width = 6, status = "primary",solidHeader=TRUE,
                                                  
                                                  div(id = "group_col_div",
                                                  
                                                  selectInput("group_col", label = "Select Column to define groups: ",
                                                            choices = list(""), selected = NULL
                                                  ) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Column Selection", content = "Choose the column the will define the different group samples", 
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                                                  ),
                                                  div(id = "ctl_group_div",
                                                  selectInput("ctl_group", label = "Select the Control Group:",
                                                              choices = list("")
                                                  )
                                                  ),
                                                  div(id = "case_group_div",
                                                  selectInput("case_group", label = "Select the Case Group:",
                                                              choices = list("")
                                                  )
                                                  ),
                                                  div(id = "sampleName_col_div",
                                                  selectInput("sampleName_col", label = "Select Columns to define unique samples: ",
                                                              choices = list(""), selected = NULL, multiple = T
                                                  ) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Columns Selection", content = "Choose the column or multiple columns that will uniquely define each individual sample", 
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                                                  ),
                                                  
                                                  
                              ),
                              
                              shinydashboard::box(id = "qc_box", title = "QC Options", width = 6,status = "primary",solidHeader=TRUE,
                                                  div(id = "max_qc_div",
                                                  sliderInput("max_sig_qc", "Max Exposure Signal", min = 1, max = 100, value = 5, step = 1) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Min Max Exposure", content = "Used to filter out peptides that have low signals. This will check the signals at max exposure and filter peptides that have lower signals that this value", 
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                                                  ),
                                                  div(id = "r2_qc_div",
                                                  sliderInput("r2_qc", "Min R2", min = 0, max = 0.99, value = 0.90, step = 0.05) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Min R2", content = "Used to filter out peptides that have nonlinear fit. This will check the R2 of the linear model and filter peptides that have R2 lower than this value", 
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                                                  )
                              )

                            ),
                            
                            fluidRow(
                              shinydashboard::box(id = "lfc_box",title = "Fold Change Options",width = 6, status = "primary",solidHeader=TRUE,
                                                  sliderInput("lfc_thr", "Log2 Fold Change Cutoff", min = 0, max = 5, value = 0.2, step = 0.05) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "LFC Cutoff", content = "Used as the cutoff of the log2 fold change values to determine the differentially phosphorylated peptides.", 
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                                                  ,
                                                  # byChip UI option
                                                  # switchInput("by_chip", "By Chip?", value = FALSE, onLabel = "Yes", offLabel = "No") %>% 
                                                  #   shinyInput_label_embed(
                                                  #     htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                  #       bs_embed_popover(
                                                  #         title = "By Chip", content = "Used to run the perumation test within each chip. Your data must have more than one chip and your group samples are found within each chip", 
                                                  #         placement = "left", trigger = "focus"
                                                  #       )
                                                  #   )
                              ),
                              shinydashboard::box(id = "sampling_box",title = "Sampling Options", width = 6,status = "primary",solidHeader=TRUE,
                                                  div(id = "itr_div",
                                                  sliderInput("itr_num", "Number of Iterations", min = 100, max = 2000, value = 500, step = 100) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Iteration", content = "The number of iterations for the permutation test", 
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                                                  
                                                  ),
                                                  div(id = "seed_div",
                                                  switchInput("use_seed", "Use Seed?", value = FALSE, onLabel = "Yes", offLabel = "No") %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Seed", content = "An option to run the analysis with a seed number to reproduce the results of the permutation test", 
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                                                    ),
                                                  div(id = "seed_num_div",
                                                  numericInput("use_seed_num", "Input Seed Number: ", value = 123) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Seed", content = "A seed number that can be used later to reproduce the results of the permutation test",
                                                          placement = "left", trigger = "focus"
                                                        )
                                                    )
                              )
                              )
                            ),
                            
                            
                            fluidRow(
                              column(12,align="center",
                                     actionButton("start_krsa", "Run KRSA", icon("paper-plane"), 
                                                  style="color: #fff; background-color: #7190a7; border-color: #2e6da4"),
                                     tags$hr(),span(textOutput("err2"), style="color:red")
                              )
                            )   
                            
                   ),
                   
                   # Step3: Results Overview Panel -----
                   tabPanel("Results: Overview", 
                            # fluidRow(column(offset = 11, width = 1, align="center",
                            #                 actionButton("help3", icon = icon("question"), label = "Help"),
                            #                 hr()
                            # )),
                            tabsetPanel(id = "res_tabs",
                              tabPanel("Summary", 
                                       fluidRow(
                                         column(width = 11,infoBoxOutput("summary_options")),
                                         column(width = 1, align="center",
                                                actionButton("help3", icon = icon("question"), label = "Help"),
                                                hr()
                                         )
                                       ),
                                       
                                       fluidRow(
                                        shinydashboard::box(id = "peps_box",title = "Peptides Selection",width = 12, status = "primary",solidHeader=TRUE,
                                        valueBoxOutput("init_peps", width = 3) %>% 
                                          bs_embed_tooltip(title = "Initial number of peptides in the input file"),
                                        valueBoxOutput("qc_maxSig_peps", width = 3) %>% 
                                          bs_embed_tooltip(title = "Number of peptides that passed QC 1: Min signal"),
                                        valueBoxOutput("qc_r2_peps", width = 3) %>% 
                                          bs_embed_tooltip(title = "Number of peptides that passed QC 2: R2 values"),
                                        valueBoxOutput("lfc_peps", width = 3) %>% 
                                          bs_embed_tooltip(title = "Number of peptides that passed LFC cutoff value"),
                                       )),
                                       
                                       fluidRow(
                                         shinydashboard::box(id = "lfc_tbl_box",title = "LFC Table",width = 6, status = "primary",solidHeader=TRUE,collapsible = T, collapsed = F,
                                                             dataTableOutput("lfc_table"),
                                                             downloadButton("lfc_table_download", label = "Save Table")
                                         ),
                                         shinydashboard::box(id = "model_tbl_box",title = "Model Table", width = 6,status = "primary",solidHeader=TRUE,collapsible = T, collapsed = F,
                                                             dataTableOutput("model_table"),
                                                             downloadButton("model_table_download", label = "Save Table")
                                         )
                                       )
                                       
                                       
                                       ),
                              tabPanel("Heatmap", 
                                       fluidRow(
                                         shinydashboard::box(id = "hm_options",title = "Options",width = 4, status = "primary",solidHeader=TRUE,
                                                             selectInput("heatmap_op1", "Data",
                                                                         choices = c("Normal", "Normalized"), label = "Data") %>% 
                                                               shinyInput_label_embed(
                                                                           htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                                             bs_embed_popover(
                                                                               title = "Data", content = "Here you can choose either the raw or normalized (by chip/barcode) data",
                                                                               placement = "top", trigger = "focus"
                                                                             )
                                                                         ),
                                                             selectInput("heatmap_op2", "Clustering Method",
                                                                         choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid"),
                                                                         selected = "ward.D2"
                                                             ) %>% 
                                                               shinyInput_label_embed(
                                                                 htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                                   bs_embed_popover(
                                                                     title = "Clustering Method", content = "Here you can choose the unsupervised clustering method that will be passed to hclust function to generate the heatmap",
                                                                     placement = "top", trigger = "focus"
                                                                   )
                                                               ),
                                                             selectInput("heatmap_op4", "Scale",
                                                                         choices = c("row", "column", "none"),
                                                                         selected = "row"
                                                             ) %>% 
                                                               shinyInput_label_embed(
                                                                 htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                                   bs_embed_popover(
                                                                     title = "Scale", content = "Here you can scale the values by row (peptide) or column (sample). 
                                                                     This will convert the values to z scores and will be displayed on the heatmap",
                                                                     placement = "top", trigger = "focus"
                                                                   )
                                                               ),
                                                             switchInput("heatmap_op3", "Cluster", onLabel = "Yes", offLabel = "No")
                                         ),
                                         shinydashboard::box(title = "Heatmap", width = 8,status = "primary",solidHeader=TRUE,
                                                             plotOutput("heatmap", height = "700px")
                                         )
                                       )

                                       
                                       ),
                              tabPanel("Boxplot", 
                                       plotOutput("boxplot") %>% withSpinner(color="#0dc5c1")
                                       ),
                              tabPanel("Cuvres Plot", 
                                       selectInput("crves_plot_opt1", "Peptide", choices = "") %>% 
                                         shinyInput_label_embed(
                                           htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                             bs_embed_popover(
                                               title = "Peptide", content = "Here you can select different peptides to plot",
                                               placement = "right", trigger = "focus"
                                             )
                                         ),
                                       plotOutput("cuvres_plot", height = "600px") %>% withSpinner(color="#0dc5c1")
                                       ),
                              tabPanel("WaterFall", 
                                       plotOutput("waterfall", height = "700px") %>% withSpinner(color="#0dc5c1")
                                       )
                            )
                            
                            ),
                   
                   # Kinase Panel -----
                   tabPanel("Results: Kinases", 
                            
                            tabsetPanel(
                              tabPanel("Table",
                                       dataTableOutput("krsa_table")
                                       ),
                              
                              tabPanel("Histogram",
                                       selectInput("histogram_opt1", "Kinase", choices = "") %>% 
                                         shinyInput_label_embed(
                                           htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                             bs_embed_popover(
                                               title = "Kinase", content = "Here you can choose different kinases to plot",
                                               placement = "right", trigger = "focus"
                                             )
                                         ),
                                       plotOutput("histogram")
                                       ),
                              tabPanel("Reverse KRSA",
                                       
                                       sliderInput("ReverseKRSA_opt1", label = "Z Score cuttof",min = 0.5, max = 5, step = 0.25, value = 2) %>% 
                                         shinyInput_label_embed(
                                           htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                             bs_embed_popover(
                                               title = "Z score Threshold", content = "Here you can control the Z score threshold. Higher values results in fewer kinases on the figure",
                                               placement = "right", trigger = "focus"
                                             )
                                         ),
                                       plotOutput("ReverseKRSA")
                                       )
                            )
                            ),
                   
                   # Network Panel -----
                   tabPanel("Results: Network", 
                            
                            fluidRow(
                              shinydashboard::box(title = "Network Options" ,width = 4, status = "primary",solidHeader=TRUE,
                                                  sliderInput("net_frq", label = "Frequency", min=1, max=10, value=4, step=1)  %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Frequency", content = "Here you can choose frequency threshold (minimum number of required connections between two kinases)",
                                                          placement = "bottom", trigger = "focus"
                                                        )
                                                    ), 
                                                  sliderInput("network_opt1", label = "Z score cutoff", min=0.5, max=10, value=2, step=0.25)  %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Z Score", content = "Here you can control the Z score threshold of kinases to be used an input in the network model",
                                                          placement = "top", trigger = "focus"
                                                        )
                                                    ), 
                                                  sliderInput("nodeSize", label = "Node Size", min=1, max=10, value=3, step=1) %>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Z Score", content = "Scale node (circle) size",
                                                          placement = "top", trigger = "focus"
                                                        )
                                                    ),
                                                  
                                                  sliderInput("nodeTextSize", label = "Node Text Size", min=1, max=10, value=6, step=1) %>% 
                                                    shinyInput_label_embed(
                                                    htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                      bs_embed_popover(
                                                        title = "Z Score", content = "Scale text size",
                                                        placement = "top", trigger = "focus"
                                                      )
                                                  ),
                                                  selectInput("layout", "Choose layout:", choices=c("Circle", "Fruchterman-Reingold", "Kamada Kawai", "LGL"))%>% 
                                                    shinyInput_label_embed(
                                                      htmltools::tags$a(shiny::icon(name = "info-circle"), href = "javascript:;") %>%
                                                        bs_embed_popover(
                                                          title = "Z Score", content = "Select layout method. This will be passed to the igraph function to generate the network",
                                                          placement = "top", trigger = "focus"
                                                        )
                                                    ),
                                     
                              ),
                              shinydashboard::box(title = "Network", width = 8, status = "primary",solidHeader=TRUE,
                                                  plotOutput("network", width = "100%", height="800px"),
                                                  downloadButton('downloadDataN', 'Download Network')
                              )
                            )
                            )
                   
))