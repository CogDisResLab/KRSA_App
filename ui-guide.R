library(cicerone)

guide1 <- Cicerone$
  new()$
  step(
    el = "box1_sig",
    title = "Signal File",
    description = "This is the main input file. Here you can upload your BioNavigator Median_SigmBg crosstab file"
  )$
  step(
    el = "load_ex_data",
    title = "Use Example Data",
    description = "Click this button to use an example dataset"
  )$
  step(
    el = "box2_map",
    title = "Kinase-Substrate File",
    description = "This is the kinase-substrate input file. 
    Here you can upload your own mapping file or use the built-in mapping file in KRSA.
    This mapping file will be used for the upstream kinase analysis"
  )$
  step(
    el = "load_map_file",
    title = "Use Built-in Kinase-Substrate Mapping File",
    description = "Click this button to use the built-in kinase-substrate mapping file"
  )$
  # step(
  #   el = "sig_prev_box",
  #   title = "Use Built-in Kinase-Substrate Mapping File",
  #   description = "Click this button to use the built-in kinase-substrate mapping file"
  # )$
  # step(
  #   el = "map_prev_box",
  #   title = "Use Built-in Kinase-Substrate Mapping File",
  #   description = "Click this button to use the built-in kinase-substrate mapping file"
  # )$

  step(
    el = "input_step_btn",
    title = "Next Step",
    position = "top",
    description = "Click this button to go to the next step of the analysis. 
    The next step will is choosing the analysis paramters"
  )
  # step(
  #   el = "[data-value='Step2: Design Options']",
  #   title = "Next Step",
  #   position = "bottom",
  #   description = "Click this button to go to the next step of the analysis. 
  #   The next step will is choosing the analysis paramters",
  #   is_id = FALSE
  # )$
  # step(
  #   el = "design_box",
  #   title = "Design Options",
  #   position = "right",
  #   description = HTML("<p>This box contains the options to select your samples and groups for the analysis.</p>
  #   <p>These options include:
  #   Grouping your samples, selecting the 'Control' and 'Case' groups, and assign unique ids to the individual samples")
  # )


guide2 <- Cicerone$
  new()$
  step(
    el = "design_box",
    title = "Design Options",
    position = "right",
    description = HTML("<p>This box contains the options to select your samples and groups for the analysis.</p>
    <p>These options include:
    Grouping your samples, selecting the 'Control' and 'Case' groups, and assign unique ids to the individual samples")
  )$
  step(
    el = "group_col_div",
    title = "Grouppig",
    position = "right",
    description = HTML("<p>Grouping samples based on the different factors that exist in the input file. 
    Use a groupping factor that will results in defing the two groups you would like analyze</p>
    If you're using the example dataset, choose 'Comment 1'
    ")
  )$
  step(
    el = "ctl_group_div",
    title = "Control Group",
    position = "right",
    description = HTML("<p>Select the 'control' group that will be used as the baseline/reference group</p>
                       If you're using the <u>example dataset</u>, choose 'M' here")
  )$
  step(
    el = "case_group_div",
    title = "Control Group",
    position = "right",
    description = HTML("<p>Select the 'case' group that will be used as the experimental group</p>
                       If you're using the <u>example dataset</u>, choose 'F' here
                       ")
  )$
  step(
    el = "sampleName_col_div",
    title = "Sample IDs",
    position = "right",
    description = HTML("<p>Choose the column or multiple columns that will uniquely define each individual sample.</p>
                      If you're using the <u>example dataset</u>, choose 'SampleName' here
                      ")
  )$
  step(
      el = "lfc_box",
      title = "LFC Threshold",
      position = "top",
      description = "Here you can adjust the log2 fold chnage (LFC) threshold. 
      Selecting a higher value is a more conservative approach and will results in fewer peptides 
      deemed to be differentially phosphorylated"
    )$
  step(
    el = "qc_box",
    title = "QC Options",
    position = "left",
    description = HTML("<p>Here you can adjust the quality control options. 
    These options will filter out peptides if they didn't pass these qc checks.</p>")
  )$
  step(
    el = "max_qc_div",
    title = "Min Signal Threshold",
    position = "left",
    description = "Select the value to be used to filter out peptides that have low signals to minimize noise. 
    This will check the signals at max exposure (usually 200ms) and filter out peptides that have lower signals that this value"
  )$
  step(
    el = "r2_qc_div",
    title = "Min R2 Threshold",
    position = "left",
    description = "Select the value to be used used to filter out peptides that have nonlinear fit. 
    This will check the R2 of the linear model and filter peptides that have R2 lower than this value"
  )$
  step(
    el = "sampling_box",
    title = "Sampling Options",
    position = "left",
    description = "Here you can adjust the sampling options. 
    These options icludes number of iterations and optional seed number"
  )$
  step(
    el = "itr_div",
    title = "Iterations",
    position = "left",
    description = "Here you can adjust the number of iterations of the permutation test. 
    The default value (2000 iterations) is a good value to select, but this could be adjusted if necessary"
  )$
  step(
    el = "seed_div",
    title = "Seed",
    position = "left",
    description = "Here is an option to run the analysis with a seed number to reproduce the results of the permutation test"
  )$
  step(
    el = "seed_num_div",
    title = "Sampling Options",
    position = "left",
    description = "If you selected the seed option, here oyou can enter seed number (must be a positive integer) that can be used later to reproduce the results of the permutation test"
  )$
  step(
    el = "start_krsa",
    title = "Run RKSA",
    position = "top",
    description = "Click this button to run KRSA! (This step will take time to process, see progress bar in the bottom right)"
  )

guide3 <- Cicerone$
  new()$
  step(
    el = "summary_options",
    title = "Parameters Summary",
    position = "right",
    description = "Here are parameters summary that were selected from the previous step"
  )$
  step(
    el = "peps_box",
    title = "Peptides Summary",
    position = "top",
    description = HTML("<p>Here are peptides selection summary (Peptides that passed QC and LFC threshold)</p>
    
    Hover over each box for more info")
  )$
  step(
    el = "lfc_tbl_box",
    title = "LFC Table",
    position = "right",
    description = "Here is the log2 fold chnage (LFC) table. 
    This table shows the LFC for each peptide based on the two group comparison chosen in the design step"
  )$
  step(
    el = "lfc_table_download",
    title = "Iniital pepetides",
    position = "right",
    description = "You can click this button to save the table as a csv file",
  )$
  step(
    el = "model_tbl_box",
    title = "Model Table",
    position = "left",
    description = "Here is the linear model fit table. 
    This table shows the results of the linear model fit for each peptide (and for each sample). It shows the slope and R2 values",
  )$
  step(
    el = "res_tabs",
    title = "Iniital pepetides",
    position = "bottom",
    description = "Click on the other tabs for more visualizations options",
  )

