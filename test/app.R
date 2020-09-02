#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            fileInput('file0', label = h3("Upload raw data file"),
                      accept=c('text/csv', 'text/comma-separated-values','text/plain', '.csv'))
        ),

        # Show a plot of the generated distribution
        mainPanel(
            shiny::dataTableOutput(outputId="contents0")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    
    
    
    data_set0 <- reactive({
        req(input$file0)
        inFile <- input$file0
        read.csv(inFile$datapath, header=T, 
                            sep="\t", quote='"', stringsAsFactors = F)
    })
    
    output$contents0 <- shiny::renderDataTable({
        data_set0()
        #}, caption = "Preview of raw data file")
    }, options=list(scrollX=TRUE, pageLength = 10))

}

# Run the application 
shinyApp(ui = ui, server = server)
