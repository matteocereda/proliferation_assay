library(shiny)
library(shinyTable)
library(matrixStats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
options(scipen=3)


# which fields get saved 
fieldsAll <- c("experiment", "num_genes", "num_reps", "num_time_points")

# which fields are mandatory
fieldsMandatory <- c("experiment", "num_genes", "num_reps", "num_time_points")

# add an asterisk to an input label
labelMandatory <- function(label) {
  tagList(label,span("*", class = "mandatory_star"))
}


# CSS to use in the app
appCSS <-
".mandatory_star { color: red; }
#error { color: red; }"


panel_width = 2
result_with = 10

shinyApp(
  ui = fluidPage(
    navbarPage("MC",

    tabPanel("Proliferation_assay",
             titlePanel("Proliferation assay curve"),
             p("This app is develop to help people to draw curves with R"),
             shinyjs::useShinyjs(),
             shinyjs::inlineCSS(appCSS),
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
             ),
             shinyjs::hidden(
               div(
                 id = "form_genes",
                 br(), 
                 uiOutput("genePanel"),
                 actionButton("submit_gene", "GO!", class = "btn-primary") 
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
            )
 ,
 tabPanel("Expression",
          titlePanel("Expression bar plot"),
          p("This app is develop to help people to draw bars with R"),
          shinyjs::useShinyjs(),
          shinyjs::inlineCSS(appCSS),
          fluidRow(
            column(panel_width,
                   div(
                     id = "expr_form",
                     textInput("expr_experiment",         labelMandatory("Experiment"), ""),
                     numericInput("num_expr_cond",       labelMandatory("Number of conditions"), 2, min=2),
                     numericInput("num_expr_genes",        labelMandatory("Number of genes"), 1, min=1),
                     actionButton("expr_submit", "Submit", class = "btn-primary")
                   ),
                   shinyjs::hidden(
                     div(
                       id = "expr_form_genes",
                       br(), 
                       uiOutput("expr_condPanel"),
                       uiOutput("expr_genePanel"),
                       actionButton("expr_submit_gene", "GO!", class = "btn-primary") 
                     )
                   )
            )
             ,
             column(result_with,
                   shinyjs::hidden(
                     div(
                       id = "expr_form_table",
                        h4("Insert your data here"),
                        htable("expr_dataTable", colHeaders="provided", rowNames = 'provided'), br(), 
                        actionButton("expr_submit_table", "Plot data", class = "btn-primary"),
                        downloadButton('expr_downloadData', 'Download Data')
                       
                     )
                   ),
                   shinyjs::hidden(
                     div(
                       id = "expr_form_plot", br(), 
                       plotOutput("expr_assay"),
                       downloadButton('expr_downloadPlot', 'Download Plot')
                     )
                   )
             )
          )
          )
 )
 ),
  
  server = function(input, output, session) {
    
    current = reactiveValues(tbl =NULL, etbl =NULL, num_genes=NULL, num_reps=NULL, num_time_points=NULL)
    
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
    
    observeEvent(input$expr_submit, {
      print("STICAZZI")
      shinyjs::disable("expr_submit")
      shinyjs::show("expr_form_genes")
    })
    
    output$expr_condPanel <- renderUI({
      ifelse(input$num_expr_cond>0,
             return(lapply(1:input$num_expr_cond, function(i) textInput( paste0("exprCond_",i), "Insert CONDITION", ""))),
             return()
      )
    })
    output$expr_genePanel <- renderUI({
      ifelse(input$num_expr_genes>0,
             return(lapply(1:input$num_expr_genes, function(i) textInput( paste0("exprGene_",i), "Insert GENE", ""))),
             return()
      )
    })
    
    observeEvent(input$expr_submit_gene, {
        genes = c()
        for(i in grep("exprGene_", names(input)) ) genes = c(genes, input[[ names(input)[i] ]][1] )
        cond = c()
        for(i in grep("exprCond_", names(input)) ) cond = c(cond, input[[ names(input)[i] ]][1] )
        print(cond)
        current$etbl <- data.frame(matrix(0, nr=input$num_expr_cond, nc=input$num_expr_genes))
        colnames(current$etbl) = genes
        rownames(current$etbl) = cond
        print(current$etbl)
        shinyjs::show("expr_form_table")
        shinyjs::show("expr_form_plot")
        shinyjs::disable("expr_submit_genes")
    })
    
    output$expr_dataTable = renderHtable({
      if(length(grep("exprGene_", names(input))>0)){
        if(!is.null(input$expr_dataTable)){
          current$etbl<<-input$expr_dataTable
        }
        return(current$etbl)
      }else{
        return()
      }
    })
    
    output$expr_assay = renderPlot({
      if(input$expr_submit_table>0){
        print("ACHI")
        data = current$etbl
        for(i in 1:ncol(data)) data[,i] = as.numeric(data[,i])
        print(str(data))
        cn = colnames(data)
        data = data %>% mutate(treatment = rownames(data))
        data = melt(data, measure.vars = cn, id.var="treatment")
        colnames(data)[2] = "gene"
        data$value = as.numeric(data$value)
        pe = 
          ggplot(data, aes(x = treatment, y = signif(value, 1), group = gene, fill = gene, label = value, ymax = max(value) * 1.05)) + 
          geom_bar(position = "dodge", stat = "identity") + 
          geom_text(position = position_dodge(width = 1), vjust = -1) +
          labs(x = "", y = "Normalized Ratio", title = "") + 
          guides(fill = guide_legend(title = "")) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) + 
          theme(legend.justification = c(1,1), legend.position = c(1,1))
        
        current$pe<<-pe
        return(pe)  
      }else{
        return()
      }
      
    })
    
    output$expr_downloadData <- downloadHandler(
      filename = function() { paste(input$expr_experiment, '.tsv', sep='') },
      content = function(file) {
        write.table(current$etbl, file, sep="\t",quote = F, row.names = F)
      }
    )
    
    output$expr_downloadPlot <- downloadHandler(
      filename = function() { paste(input$expr_experiment, '.pdf', sep='') },
      content = function(file) {
        pdf(file,w=12,h=8)
        print(current$pe)
        dev.off()
      })
    
    
    
  }

)