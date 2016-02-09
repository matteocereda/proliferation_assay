rea
shinyServer(function(input, output, session) {
  
 current = reactiveValues(gene=NULL, replicate=NULL, experiment=NULL, tbl=NULL, p=NULL)
#   observe({
#     # if(input$update > 0) {
#       current$replicate = input$replicate
#       current$gene = input$gene
#     # }
#   })
  
 validate <- function(tbl){
   updateTableStyle(session, "tbl", "valid", which(as.numeric(tbl$num2) < 50), 2)
   updateTableStyle(session, "tbl", "warning",  which(as.numeric(tbl$num2) >= 50 & as.numeric(tbl$num2) < 100), 2)
   updateTableStyle(session, "tbl", "invalid",  which(as.numeric(tbl$num2) >= 100), 2)    
 }
 
 
 
  output$ogene = renderText({return(input$submit)})
  
  Data = reactive({
    if (input$submit > 0) {
      tbl <- data.frame(matrix(0, nr=input$timepoint, nc=input$replicate*2))
      colnames(tbl) = c(rep("Negative",input$replicate), rep(input$gene,input$replicate))
      current$tbl <<- tbl
      return(current$tbl)
    }else{
      if (is.null(input$tbl)){
        tbl <- data.frame(matrix(0, nr=input$timepoint, nc=input$replicate*2))
        colnames(tbl) = c(rep("Negative",input$replicate), rep(input$gene,input$replicate))
        current$tbl <<- tbl
      } else{
        current$tbl <<- input$tbl
      }
      print(current$tbl)
      return(current$tbl)
    }
  })
 

  output$tbl <- renderHtable({Data()})
#   output$tbl <- renderHtable({
#     if (is.null(input$tbl)){
#       tbl <- data.frame(matrix(0, nr=input$timepoint, nc=input$replicate*2))
#       colnames(tbl) = c(rep("Negative",input$replicate), rep(input$gene,input$replicate))
#       current$tbl <<- tbl
#     } else{
#       current$tbl <<- input$tbl
#     }
#     print(current$tbl)
#     return(current$tbl)
#   })  
  
  output$assay = renderPlot({
    
    print(current$tbl)
    data = current$tbl
    rownames(data) = as.character(c(24,48,76,92))
    for(i in 1:ncol(data)) data[,i] = as.numeric(data[,i])
    
    print("data")
    print(str(data))
    plot_data = sapply(unique(colnames(data)), function(x) rowMeans(data[, colnames(data) == x, drop = FALSE], na.rm = TRUE))
    
    # Calculate standard deviations for the replicates which will be used to draw the error bars in the plot
    data = as.matrix(data)
    st_devs = sapply(unique(colnames(data)), function(x) rowSds(data[, colnames(data) == x, drop = FALSE], na.rm = TRUE))
    st_devs_long = st_devs %>% as.data.frame() %>% gather(key = treatment, value = sd, 1:ncol(st_devs))
    
    # Draw the plot, which will be produced at the same directory
    # Each row does a job which is self-explicatory, and options for each job can be changed as desired
    p = plot_data %>% data.frame() %>%
      mutate(time_point = rownames(plot_data)) %>%
      gather(key = treatment, value = RFU, 1:ncol(plot_data)) %>%
      mutate(sd = st_devs_long %>% .$sd) %>%
      ggplot(aes(x = time_point, y = RFU, group = treatment, colour = treatment, shape = treatment, ymin = RFU - sd, ymax = RFU + sd)) +
      geom_line() +
      geom_point(size = 5) +
      geom_errorbar(width = 0.1) +
      labs(x = "hours", y = "Normalized RFU", title = input$experiment) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
      theme(legend.justification = c(0,1), legend.position = c(0,1))
    current$p<<-p
      return(p)  
      
  })
  
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
    }
  )
 
  
#   output$tbl <- renderHtable({
#      cols <- input$replicate*2
#      tbl<- data.frame(matrix(0, nr=4, nc=cols))
#      colnames(tbl) = c(rep("Negative",input$replicate), rep(input$gene,input$replicate)) 
#      return(tbl)
#    })
#   
#    output$clickText <- renderPrint({
#     input$tblClick
#   })
})

# Define server logic required to draw a histogram
# shinyServer(function(input, output) {
#   
#   # Expression that generates a histogram. The expression is
#   # wrapped in a call to renderPlot to indicate that:
#   #
#   #  1) It is "reactive" and therefore should re-execute automatically
#   #     when inputs change
#   #  2) Its output type is a plot
# 
#   observe({
#     if (input$test == 0) 
#       return()
#     isolate({ 
#       output$value <-renderTable( 
#         mdat <- matrix(NA, nrow = input$test, ncol = 2, byrow = TRUE) 
#       )
#       
# #   output$distPlot <- renderPlot({
# #     x    <- faithful[, 2]  # Old Faithful Geyser data
# #     # bins <- seq(min(x), max(x), length.out = input$bins + 1)
# #     
# #     # draw the histogram with the specified number of bins
# #     hist(x, breaks = 10, col = 'darkgray', border = 'white')
# #   }
#   )
# })