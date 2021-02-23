# input files -----
input_file_modal <- bs_modal(
  id = "input_file",
  title = "Input File",
  size = "large",
  body = includeMarkdown(system.file("markdown", "modal.md", package = "bsplus"))
)