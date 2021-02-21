library(shiny)
library(tidyverse)
library(reactable)
library(KRSA)


server <- function(input, output, session) {
  
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
    data_file$data <- krsa_read(inFile$datapath)
    
  })

  observe({
    req(input$map_file)
    inFile <- input$map_file
    map_file$final_map <- read_delim(inFile$datapath, delim = "\t")
    
  })
  
  ## preview files -----
  output$sig_tbl_preview <- renderTable({
    req(data_file$data)
    data_file$data %>% head(6)
  })
  
  output$map_tbl_preview <- renderTable({
    req(map_file$final_map)
    map_file$final_map %>% head(5)
  })
  
  ## start input btn ----
  observeEvent(input$input_step_btn, {
    
    output$err <- renderText(
      validate(
        need(data_file$data, '- You need to load the signal data!'),
        need(map_file$final_map, '- You need to load the kinase-substarte mapping file!')
      )
    )
    
    req(data_file$data, map_file$final_map)
    
    # TODO: check if is STK or PTK
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
    
    
    
    
    showTab(inputId = "tabs", target = "Step2: Design Options")
    updateNavbarPage(session, "tabs", selected = "Step2: Design Options")
    
    
    
    
    
  })
  
  # Step2: starts KRSA -----
  observeEvent(input$start_krsa, {
    
    # Catch Errors for KRSA -----
    output$err2 <- renderText(
      validate(
        need(input$ctl_group != input$case_group, '- The Control and Case Groups need to be different!')
      )
    )
    
    req(input$ctl_group != input$case_group)
    
    chiptype <- "STK"

    showTab(inputId = "tabs", target = "Results: Overview")
    showTab(inputId = "tabs", target = "Results: Kinases")
    showTab(inputId = "tabs", target = "Results: Network")
    
    updateNavbarPage(session, "tabs", selected = "Results: Overview")
    
    output$summary_options <- renderInfoBox({
      infoBox(
        "Design Options Selected: ", 
        HTML(paste(
          # TODO update chip type based on input
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
    
    data_file$filtered_data <- dplyr::filter(data_file$data, get(input$group_col) %in% c(input$ctl_group, input$case_group))
    
    data_file$pep_init <- krsa_filter_ref_pep(data_file$filtered_data$Peptide %>% unique()) 
    
    data_file$filtered_data %>% 
      dplyr::filter(Peptide %in% data_file$pep_init) -> data_file$filtered_data
    
    
    data_file$filtered_data <- krsa_qc_steps(data_file$filtered_data, sat_qc = F)
    
    #TODO need to dynamically create new names
    data_file$filtered_data %>% 
      mutate(Group = get(input$group_col), SampleName = paste(Group,get(input$sampleName_col), sep = "_")) -> data_file$filtered_data
    
    hgf <<- data_file$filtered_data
    
    data_file$pw_200 <- krsa_extractEndPointMaxExp(data_file$filtered_data, chiptype)
    data_file$pw <- krsa_extractEndPoint(data_file$filtered_data, chiptype)
  
    
    data_file$pep_maxSig <- krsa_filter_lowPeps(data_file$pw_200, input$max_sig_qc)
    
    
    data_file$data_modeled <- krsa_scaleModel(data_file$pw, data_file$pep_maxSig)

    data_file$pep_nonLinear <- krsa_filter_nonLinear(data_file$data_modeled$scaled, input$r2_qc)

    data_file$lfc_table <- krsa_group_diff(data_file$data_modeled$scaled,
                                           c(input$case_group, input$ctl_group),
                                           data_file$pep_nonLinear, byChip = F)
    
    data_file$pep_sig <- krsa_get_diff(data_file$lfc_table, LFC, lfc_thr = input$lfc_thr)
    
    jhg <<- data_file$data_modeled$normalized
    
    jhg$Peptide %>% head(10) %>% unique() -> ppp
    
    krsa_heatmap(jhg, ppp)

    
    output$init_peps <- renderValueBox({
      valueBox(
        paste0(length(data_file$pep_init)), "Initial Peptides", icon = icon("list-ul")
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
    
    output$heatmap <- renderPlot({
      krsa_heatmap(
        ifelse(input$heatmap_op1 == "Normalized", data_file$data_modeled$normalized, 
               data_file$data_modeled$scaled
               ), 
        data_file$pep_sig)
    })
    
  })
  
  ## Figures -----
  
  ### Heatmap -----
  
  # Step4: Kinase Analysis -----
  # Step5: Network -----
  
  
    


  
}