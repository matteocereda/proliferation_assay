library(shiny)
library(shinyTable)
library(matrixStats)
library(dplyr)
library(tidyr)
library(ggplot2)
options(scipen=3)


# which fields get saved 
fieldsAll <- c("experiment", "num_genes", "num_reps", "num_time_points")

# which fields are mandatory
fieldsMandatory <- c("experiment", "num_genes", "num_reps", "num_time_points")

# add an asterisk to an input label
labelMandatory <- function(label) {
  tagList(label,span("*", class = "mandatory_star"))
}

# get current Epoch time
epochTime <- function() { return(as.integer(Sys.time())) }

# get a formatted string of the timestamp (exclude colons as they are invalid
# characters in Windows filenames)
humanTime <- function() {  format(Sys.time(), "%Y%m%d-%H%M%OS")}

# save the results to a file
saveData <- function(data) {
  fileName <- sprintf("%s_%s.csv",humanTime(), digest::digest(data))
  write.csv(x = data, file = file.path(responsesDir, fileName), row.names = FALSE, quote = TRUE)
}

# load all responses into a data.frame
loadData <- function() {
  files <- list.files(file.path(responsesDir), full.names = TRUE)
  data <- lapply(files, read.csv, stringsAsFactors = FALSE)
  #data <- dplyr::rbind_all(data)
  data <- do.call(rbind, data)
  data
}

# directory where responses get stored
responsesDir <- file.path("responses")

# CSS to use in the app
appCSS <-
".mandatory_star { color: red; }
.shiny-input-container { margin-top: 25px; }
#submit_msg { margin-left: 15px; }
#error { color: red; }
#adminPanel { border: 4px solid #aaa; padding: 0 20px 20px; }"

# usernames that are admins
adminUsers <- c("admin", "prof")

panel_width = 2
result_with = 10
shinyApp(
  ui = fluidPage(
    shinyjs::useShinyjs(),
    shinyjs::inlineCSS(appCSS),

    titlePanel("Proliferation assay curve"),
    p("This app is develop to help people to draw curves with R"),
    
    fluidRow(
      column(panel_width,
             div(
               id = "form",
               textInput("experiment",         labelMandatory("Experiment"), ""),
               numericInput("num_genes",       labelMandatory("Number of conditions"), 2, min=2),
               numericInput("num_reps",        labelMandatory("Number of replicates"), 1, min=1),
               numericInput("num_time_points", labelMandatory("Number of time points"), 2, min = 2),
               numericInput("len_time_point",  labelMandatory("Time point length"), 24, min = 1),
               actionButton("submit", "Submit", class = "btn-primary")
               # ,
#                shinyjs::hidden(
#                  span(id = "submit_msg", "Submitting..."),
#                  div(id = "error",
#                      div(br(), tags$b("Error: "), span(id = "error_msg"))
#                  )
#                )
             ),
             shinyjs::hidden(
               div(
                 id = "form_genes",
                 br(), 
                 uiOutput("genePanel"),
                 actionButton("submit_gene", "GO!", class = "btn-primary")               #         DT::dataTableOutput("responsesTable") 
               )
      )
      ),
      column(result_with,
             
             shinyjs::hidden(
               div(
                 id = "form_table",
                 h4("Insert your data here"),
                 htable("dataTable", colHeaders="provided", rowNames = 'provided'),
                 br(), 
                 actionButton("submit_table", "Plot data", class = "btn-primary"),
                 downloadButton('downloadData', 'Download Data')

               )
             ),

            shinyjs::hidden(
              div(
                id = "form_plot",
                br(), 
                plotOutput("assay"),
                downloadButton('downloadPlot', 'Download Plot')
              )
            )
      )
    )
    
    
  ),
  
  server = function(input, output, session) {
    
    current = reactiveValues(tbl =NULL, num_genes=NULL, num_reps=NULL, num_time_points=NULL)
    
    # Enable the Submit button when all mandatory fields are filled out
    observe({
      mandatoryFilled <- vapply(fieldsMandatory,   function(x) { !is.null(input[[x]]) && input[[x]] != "" && input[[x]] != 0  }, logical(1))
      mandatoryFilled <- all(mandatoryFilled)
      shinyjs::toggleState(id = "submit", condition = mandatoryFilled)
    })
    
     
    # When the Submit button is clicked, submit the response
    observeEvent(input$submit, {
      
      # User-experience stuff
      shinyjs::disable("submit")
      shinyjs::hide("error")

      # Save the data (show an error message in case of error)
      tryCatch({
        shinyjs::show("form_genes")
        shinyjs::disable("submit")
      # },
#       error = function(err) {
#         shinyjs::html("error_msg", err$message)
#         shinyjs::show(id = "error", anim = TRUE, animType = "fade")
#       },
#       finally = {
#         shinyjs::enable("submit")
#         shinyjs::hide("submit_msg")
      })
    })
    
    # submit another response
    observeEvent(input$submit_gene, {
      tryCatch({
        genes = c()
        for(i in grep("Gene_", names(input)) ) genes = c(genes, input[[ names(input)[i] ]][1] )
        current$tbl <- data.frame(matrix(0, nr=input$num_time_points, nc=input$num_reps*input$num_genes))
        colnames(current$tbl) = rep(genes, each=input$num_reps)
        rownames(current$tbl) = as.character( seq(input$len_time_point, input$len_time_point*input$num_time_points, by=input$len_time_point) )
        shinyjs::show("form_table")
        shinyjs::show("form_plot")
        shinyjs::disable("submit_genes")
#       },
#       error = function(err) {
#         shinyjs::html("error_msg", err$message)
#         shinyjs::show(id = "error", anim = TRUE, animType = "fade")
#       },
#       finally = {
#         shinyjs::enable("submit")
#         shinyjs::hide("submit_msg")
      })
    })
    
    # render the form_gene
    output$genePanel <- renderUI({
      ifelse(input$num_genes>0,
       return(lapply(1:input$num_genes, function(i) textInput( paste0("Gene_",i), "Insert condition (or gene) name", ""))),
       return()
      )
      })
    
    output$dataTable = renderHtable({
      if(length(grep("Gene_", names(input))>0)){
        if(!is.null(input$dataTable)){
          current$tbl<<-input$dataTable
        }
        return(current$tbl)
      }else{
        return()
      }
    })
    
    observeEvent(input$submit_table, {
      tryCatch({
        shinyjs::show("form_plot")
        shinyjs::disable("submit_table")
      },
      error = function(err) {
        shinyjs::html("error_msg", err$message)
        shinyjs::show(id = "error", anim = TRUE, animType = "fade")
      },
      finally = {
        shinyjs::enable("submit")
        shinyjs::hide("submit_msg")
      })
    })
    
    output$assay = renderPlot({
      if(input$submit_table>0){

      data = current$tbl
      for(i in 1:ncol(data)) data[,i] = as.numeric(data[,i])
      print(str(data))
      plot_data = sapply(unique(colnames(data)), function(x) rowMeans(data[, colnames(data) == x, drop = FALSE], na.rm = TRUE))
      
      # Calculate standard deviations for the replicates which will be used to draw the error bars in the plot
      data = as.matrix(data)
      st_devs = sapply(unique(colnames(data)), function(x) rowSds(data[, colnames(data) == x, drop = FALSE], na.rm = TRUE))
      st_devs_long = st_devs %>% as.data.frame() %>% gather(key = Condition, value = sd, 1:ncol(st_devs))
      
      # Draw the plot, which will be produced at the same directory
      # Each row does a job which is self-explicatory, and options for each job can be changed as desired
      p = plot_data %>% data.frame() %>%
        mutate(time_point = rownames(plot_data)) %>%
        gather(key = Condition, value = RFU, 1:ncol(plot_data)) %>%
        mutate(sd = st_devs_long %>% .$sd) %>%
        ggplot(aes(x = time_point, y = RFU, group = Condition, colour = Condition, shape = Condition, ymin = RFU - sd, ymax = RFU + sd)) +
        geom_line() +
        geom_point(size = 5) +
        geom_errorbar(width = 0.1) +
        labs(x = "hours", y = "Normalized RFU", title = input$experiment) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
        theme(legend.justification = c(0,1), legend.position = c(0,1))
      current$p<<-p
      return(p)  
      }else{
        return()
      }
      
    })
    
    # Allow user to download responses
    output$downloadData <- downloadHandler(
      filename = function() { paste(input$experiment, '.tsv', sep='') },
      content = function(file) {
        write.table(current$tbl, file, sep="\t",quote = F, row.names = F)
      }
    )
    output$downloadPlot <- downloadHandler(
      filename = function() { paste(input$experiment, '.pdf', sep='') },
      content = function(file) {
        pdf(file,w=12,h=8)
        print(current$p)
        dev.off()
      })
  }
)