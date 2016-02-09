library(shiny)
library(shinyTable)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Proliferation Assay Plot"),
  
  sidebarPanel(
      textInput('experiment',label = h5("Enter experiment name"), value = "Experiment"),
      textInput('gene',label = h5("Enter gene name"), value = "Gene"),
       numericInput('replicate',label = h5("Enter number of replicates"), value = 3),
       numericInput('timepoint',label = h5("Enter number of time points"), value = 4),
       actionButton("submit", "Update Table")
  ),
  
  # Show the simple table
  mainPanel(
    helpText(HTML("Insert data")),
    textOutput("ogene"),
    htable("tbl", colHeaders="provided"), # clickId="tblClick",
    downloadButton('downloadData', 'Download Data'),
    
    plotOutput("assay"),
    downloadButton('downloadPlot', 'Download Plot')
    #     h3("Current Selection"),
#     verbatimTextOutput("clickText")
  )
))
