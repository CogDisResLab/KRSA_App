library(shiny)
library(tidyverse)
library(KRSA)
library(shinydashboard)
library(broom)
library(igraph)

source("funcs/krsa_violin_plot_2.R")
source("funcs/krsa_ball_model_2.R")
source("ui-guide.R")

server <- function(input, output, session) {
  
  
  observeEvent(input$help1, {
  guide1$init()$start()
  })
  
  observeEvent(input$help2, {
    guide2$init()$start()
  })
  
  observeEvent(input$help3, {
    guide3$init()$start()
  })
  
  session$onSessionEnded(stopApp)
  
  # Define reactive values ----
  data_file <- reactiveValues()
  map_file <- reactiveValues()
  design_options <- reactiveValues()
  
  hideTab(inputId = "tabs", target = "Step2: Design Options")
  hideTab(inputId = "tabs", target = "Results: Overview")
  hideTab(inputId = "tabs", target = "Results: Kinases")
  hideTab(inputId = "tabs", target = "Results: Network")
  
  # Step1: Input -----
  
  ## load example files -----
  observeEvent(input$load_ex_data, {
    print("file_btn pressed")
    data_file$data <- krsa_read("data/datasets/DPLFC_MvsF_STK.txt")
  })
  
  observeEvent(input$load_map_file, {
    print("map_btn pressed")
    map_file$stk_v1<- KRSA_Mapping_STK_PamChip_87102_v1
    map_file$ptk_v2 <- KRSA_Mapping_PTK_PamChip_86402_v1
    
    map_file$final_map<- KRSA_Mapping_STK_PamChip_87102_v1
  })
  
  ## load actual files -----
  observe({
    req(input$input_file)
    inFile <- input$input_file
    
    data_file$data <- tryCatch(
      {
        
        krsa_read(inFile$datapath) 

      },
      error=function(cond) {
        message("Failed to parase the file")
        return(1)
      }
    )
    
  })

  observe({
    req(input$map_file)
    inFile <- input$map_file
    
    map_file$final_map <- tryCatch(
      {
        
        read_delim(inFile$datapath, delim = "\t") 
        
      },
      error=function(cond) {
        message("Failed to parase the kin-sub file")
        return(1)
      }
    )
    
  })
  
  ## preview files -----
  output$sig_tbl_preview <- renderTable({
    req(data_file$data)
      validate(
        need(data_file$data != 1, '- Failed to parse the input file'),
        need(colnames(data_file$data)[colnames(data_file$data) %in% c("ExposureTime", "Signal","Cycle", "Peptide", "SampleName")] != 5, "- One of these column is missing from the input file (ExposureTime, Signal,Cycle, Peptide, SampleName)")
        
      )

    data_file$data %>% head(6)
  })
  
  output$map_tbl_preview <- renderTable({
    req(map_file$final_map)
    validate(
      need(map_file$final_map != 1, '- Failed to parse the kinase-substrate file'),
      need(all(colnames(map_file$final_map) %in% c("Substrates", "Kinases")), "- One of these column is missing from the kinase-substrate file (Substrates, Kinases)")
    )
    
    map_file$final_map %>% head(5)
  })
  
  ## start input btn ----
  observeEvent(input$input_step_btn, {
    
    output$err <- renderText(
      validate(
        need(data_file$data, '- You need to load the signal data!'),
        need(data_file$data != 1, '- Failed to parse the input file'),
        need(colnames(data_file$data)[colnames(data_file$data) %in% c("ExposureTime", "Signal","Cycle", "Peptide", "SampleName")] != 5, "- One of these column is missing from the input file (ExposureTime, Signal,Cycle, Peptide, SampleName)"),
        need(map_file$final_map, '- You need to load the kinase-substarte mapping file!'),
        need(map_file$final_map != 1, '- Failed to parse the kinase-substarte file'),
        need(all(colnames(map_file$final_map) %in% c("Substrates", "Kinases")), "- One of these column is missing from the kinase-substrate file (Substrates, Kinases)")
      )
    )
    
    req(data_file$data, map_file$final_map, 
        !all(colnames(data_file$data) %in% c("ExposureTime", "Signal","Cycle", "Peptide", "SampleName")))
    
    group_colms <- colnames(data_file$data) %>% unique()
    group_colms <- group_colms[!group_colms %in% c("ExposureTime", "Signal","Cycle", "Peptide")]
    updateSelectInput(session, "group_col", choices = group_colms, 
                      selected = if(any("SampleName" %in% group_colms)){"SampleName"} else {group_colms[1]}
                        )
    
    observe({
      req(input$group_col)
      updateSelectInput(session, "ctl_group", choices = pull(data_file$data,input$group_col) %>% unique())
      updateSelectInput(session, "case_group", choices = pull(data_file$data,input$group_col) %>% unique())
      updateSelectInput(session, "sampleName_col", choices = group_colms[group_colms!=input$group_col])
    })
    
    observe({
      print(input$use_seed)
      toggleState(id = "use_seed_num", condition = input$use_seed == T)
    })

    showTab(session = session, inputId = "tabs", target = "Step2: Design Options")
    updateNavbarPage(session, "tabs", selected = "Step2: Design Options")
 
  })
  
  # Step2: starts KRSA -----
  observeEvent(input$start_krsa, {
    
    # TODO select samples 
    
    # Catch Errors for KRSA -----
    output$err2 <- renderText(
      validate(
        need(input$ctl_group != input$case_group, '- The Control and Case Groups need to be different!'),
        need(input$sampleName_col, '- Need unique sample name input')
      )
    )
    
    req(input$ctl_group != input$case_group, input$sampleName_col)
    
    disable("start_krsa")
    
    withProgress(message = "Loading Data", value = 0, {
      
    data_file$data$Peptide %>% unique() %>% {any(. %in% c("K6PL_766_778", "MPH6_140_152", "ACM4_456_468"))} -> check_chip
    
    if(check_chip) {
      chiptype <- "STK"
    } else {
      chiptype <- "PTK"
    }
      
    showTab(inputId = "tabs", target = "Results: Overview")
    showTab(inputId = "tabs", target = "Results: Kinases")
    showTab(inputId = "tabs", target = "Results: Network")
    
    output$summary_options <- renderInfoBox({
      infoBox(
        "Design Options Selected: ", 
        HTML(paste(
          paste0("Chip Type = ", chiptype),br(),
          paste0("Control Group = ",input$ctl_group),br(),
          paste0("Case Group = ",input$case_group),br(),
          paste0("Minimum Exposure Intensity Cutoff = ",input$max_sig_qc),br(),
          paste0("R2 Cutoff = ",input$r2_qc),br(),
          paste0("LFC Cutoff = ", input$lfc_thr),br(),
          paste0("Seed Used: ", 
                 ifelse(input$use_seed, "Yes", "No")
                 ),br(),
          ifelse(input$use_seed, paste0("Seed Number= ", round(input$use_seed_num)), "")
               )),
        
        icon = icon("list-ul"),
        color = "yellow"
      )
    })
    
    incProgress(2/10, message = "Filtering Data")
    data_file$filtered_data <- dplyr::filter(data_file$data, .data[[input$group_col]] %in% c(input$ctl_group, input$case_group))
    
    data_file$pep_init <- krsa_filter_ref_pep(data_file$filtered_data$Peptide %>% unique()) 
    
    data_file$filtered_data %>% 
      dplyr::filter(Peptide %in% data_file$pep_init) -> data_file$filtered_data

    data_file$filtered_data <- krsa_qc_steps(data_file$filtered_data, sat_qc = F)
    
    data_file$filtered_data %>% 
      mutate(Group = .data[[input$group_col]]) %>% 
      unite("SampleName",all_of(c(input$group_col,input$sampleName_col)), remove = F) -> data_file$filtered_data
    
    # group by sample
    data_file$filtered_data %>% 
      select(-Cycle, -Signal, -Peptide, -ExposureTime) %>% 
      distinct(.keep_all = T) %>% 
      group_by(SampleName) %>% 
      add_count() %>% 
      mutate(dup = if_else(n > 1, 1, 0)) %>% 
      select(-n) %>% 
      ungroup() %>% pull(dup) -> dup_req
    
    output$err2 <- renderText(
      validate(
        need(!any(dup_req == 1), '- Duplicated samples exist')
      )
    )
    
    enable("start_krsa")
    
    req(!any(dup_req == 1), input$sampleName_col)
    
    data_file$pw_200 <- krsa_extractEndPointMaxExp(data_file$filtered_data, chiptype)
    data_file$pw <- krsa_extractEndPoint(data_file$filtered_data, chiptype)
    data_file$pep_maxSig <- krsa_filter_lowPeps(data_file$pw_200, input$max_sig_qc)

    incProgress(3/10, message = "Fitting the linear model ...")
    data_file$data_modeled <- krsa_scaleModel(data_file$pw, data_file$pep_maxSig)

    data_file$pep_nonLinear <- krsa_filter_nonLinear(data_file$data_modeled$scaled, input$r2_qc)

    incProgress(4/10, message = "Calculating LFC")
    
    # data_file$data_modeled$scaled %>% 
    #   select(Barcode, Group) %>% 
    #   distinct(.keep_all = T) %>% 
    #   group_by(Barcode) %>% 
    #   count() %>% mutate(dup = if_else(n == 2, T, F)) %>% 
    #   pull(dup) %>% all(.) ->> found_inEachChip
    # 
    # print("foun_inChip")
    # print(found_inEachChip)
    # 
    # print("both_tests")
    # print(input$by_chip && found_inEachChip)
    # 
    # output$err2 <- renderText(
    # validate(
    #   # Check if more than one chip
    #   need(!input$by_chip && found_inEachChip, " - For byChip analysis, samples must be found in ech chip")
    # )
    # )
    # 
    # req(!all(input$by_chip && !found_inEachChip))
    
    data_file$lfc_table <- krsa_group_diff(data_file$data_modeled$scaled,
                                           c(input$case_group, input$ctl_group),
                                           data_file$pep_nonLinear, byChip = F)
    
    
    
    # TODO byChip option
    # if(input$by_chip) {
    #   data_file$pep_sig <- krsa_get_diff_byChip(data_file$lfc_table, LFC, lfc_thr = input$lfc_thr)
    # } else {
    #   data_file$lfc_table %>% 
    #       filter(LFC >= input$lfc_thr | LFC <= input$lfc_thr*-1) %>%
    #       pull(Peptide) -> data_file$pep_sig
    # }
    
    data_file$lfc_table %>%
      filter(LFC >= input$lfc_thr | LFC <= input$lfc_thr*-1) %>%
      pull(Peptide) -> data_file$pep_sig
    
    
    
    output$init_peps <- renderValueBox({
      valueBox(
        paste0(length(data_file$pep_init)), "Initial Peptides", icon = icon("list-ul"),
        color = "yellow"
      )
    })
    output$qc_maxSig_peps <- renderValueBox({
      valueBox(
        paste0(length(data_file$pep_maxSig)), "Max Sig Peptides", icon = icon("list-ul"),
        color = "yellow"
      )
    })
    output$qc_r2_peps <- renderValueBox({
      valueBox(
        paste0(length(data_file$pep_nonLinear)), "R2 Peptides", icon = icon("list-ul"),
        color = "yellow"
      )
    })
    output$lfc_peps <- renderValueBox({
      valueBox(
        paste0(length(data_file$pep_sig)), "LFC Peptides", icon = icon("list-ul"),
        color = "yellow"
      )
    })
    
    updateNavbarPage(session, "tabs", selected = "Results: Overview")
    
    output$lfc_table <- renderDataTable({
      data_file$lfc_table
    })
    
    output$lfc_table_download <- downloadHandler(
      filename = function() {
        paste("LFC_Table", ".txt", sep = "")
      },
      content = function(file) {
        write_delim(data_file$lfc_table, file, delim = "\t")
      }
    )
    
    output$model_table <- renderDataTable({
      data_file$data_modeled$scaled
    })
    
    output$model_table_download <- downloadHandler(
      filename = function() {
        paste("Model_Table", ".txt", sep = "")
      },
      content = function(file) {
        write_delim(data_file$data_modeled$scaled, file, delim = "\t")
      }
    )
    
    incProgress(5/10, message = "Genrating Figures")
    
    output$heatmap <- renderPlot({
      krsa_heatmap(

               data = if(input$heatmap_op1 == "Normalized") data_file$data_modeled$normalized
               else data_file$data_modeled$scaled, 
               peptides = data_file$pep_sig, 
               clustering_method = input$heatmap_op2,
               cluster_cols = input$heatmap_op3,
               scale = input$heatmap_op4
        )
    })
    
    output$boxplot <- renderPlot({
      krsa_violin_plot_2(data_file$data_modeled$scaled, 
                       data_file$pep_sig, 
                       facet = F
                       )
    })
    
    updateSelectInput(session, "crves_plot_opt1", choices = data_file$pep_sig,
                      selected = data_file$pep_sig[1]
                      )
    
    output$cuvres_plot <- renderPlot({
      krsa_curve_plot(data_file$pw, input$crves_plot_opt1)
    })
    
    output$waterfall <- renderPlot({
      krsa_waterfall(data_file$lfc_table, input$lfc_thr, byChip = F)
    })
    
    # Step4: Kinase Analysis -----
    
    incProgress(8/10, message = "Resampling analysis")
    
    #chipCov <<- KRSA_coverage_STK_PamChip_87102_v1
    #KRSA_file <<- KRSA_Mapping_STK_PamChip_87102_v1
    
    map_file$final_cov <- separate_rows(map_file$final_map, Kinases, sep = "\\s+") %>% 
      rename(Kin = Kinases)
    
    krsa(data_file$pep_sig, return_count = T, itr = input$itr_num, 
         seed = if(input$use_seed) input$use_seed_num else runif(1, 1, 100) %>% round(),
         map_file = map_file$final_map,
         cov_file = map_file$final_cov
         ) -> data_file$krsa
    
    output$krsa_table <- renderDataTable({
      data_file$krsa$KRSA_Table %>% arrange(desc(abs(Z)))
    })
    
    updateSelectInput(session, "histogram_opt1", choices = data_file$krsa$KRSA_Table$Kinase)
    
    output$histogram <- renderPlot({
      krsa_histogram_plot(data_file$krsa$KRSA_Table, data_file$krsa$count_mtx, input$histogram_opt1) +
        theme(strip.text = element_text(size = 20))
    })
    
    output$ReverseKRSA <- renderPlot({
      krsa_reverse_krsa_plot(map_file$final_cov, data_file$lfc_table,
                             filter(data_file$krsa$KRSA_Table, abs(Z) >= input$ReverseKRSA_opt1) %>% pull(Kinase)
                             ,lfc_thr = input$lfc_thr, byChip = F)
    })
    
    # Step5: Network -----
    

    output$network <- renderPlot({
      krsa_ball_model_2(
        kinase_hits = filter(data_file$krsa$KRSA_Table, abs(Z) >= input$network_opt1) %>% pull(Kinase),
        data_file$krsa$KRSA_Table,
        input$net_frq,
        input$nodeSize,
        input$nodeTextSize,
        input$layout
      )
    })
    #   
    #   
    # })
    

    })
  })
  
  
  
}