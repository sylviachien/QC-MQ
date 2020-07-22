# app_10

options(shiny.maxRequestSize=10000*1024^2)
library(BiocManager)
options(repos = BiocManager::repositories())


###################################################################################
# load libraries
library(GenomicAlignments)
library(knitr)
library(kableExtra)
library(colorspace)
library(shinyalert)
library(shinyWidgets)
library(shinyBS)
library(tools)
library(shinyjs)
library(pvca)
library(Biobase)
library(ggthemes)
library(ggplot2)
library(scales)
library(reshape2)
library(stringr)
library(quantro)
library(plyr)
library(tidyr)
library(plotly)
library(dplyr)

###################################################################################

server <- function(input, output, session){
  source("functions.R")
  shinyjs::useShinyjs()
  
  observeEvent(input$doe_demo_table, {
    demo_doe <- read.table("demo_run_file/doe_demo", sep = "\t", header = TRUE)
    demo_doe <- demo_doe[1:11,]
    output$doe_example <- renderTable(demo_doe, striped = TRUE, bordered = TRUE)
  })
  
  ##################################
  # quick check for uploaded files #
  ##################################
  qc_allow <- 0
  
  ### 1 proteinGroups.txt
  observeEvent(input$pg_file, {
    pg_file <- input$pg_file$datapath
    pg_name <- input$pg_file$name
    if(pg_name!="proteinGroups.txt"){
      createAlert(session, anchorId = "alert1", alertId = "a1", style = "danger", content =  paste0("The uploaded file '", pg_name, "' is NOT 'proteinGroups.txt'!"))
      qc_allow_pg = qc_allow + 1
    }else{
      shinyBS::closeAlert(session, "a1")
      qc_allow_pg = 0
    }
  })
  ### 2 peptides.txt
  observeEvent(input$pt_file, {
    pt_file <- input$pt_file$datapath
    pt_name <- input$pt_file$name
    if(!is.null(pt_file) && pt_name!="peptides.txt"){
      createAlert(session, anchorId = "alert2", alertId = "a2", style = "danger", content =  paste0("The uploaded file '", pt_name, "' is NOT 'peptides.txt'!"))
    }else{
      shinyBS::closeAlert(session, "a2")
    }
  })
  ### 3 msms.txt
  observeEvent(input$ms_file, {
    ms_file <- input$ms_file$datapath
    ms_name <- input$ms_file$name
    if(!is.null(ms_file) && ms_name!="msms.txt"){
      createAlert(session, anchorId = "alert3", alertId = "a3", style = "danger", content =  paste0("The uploaded file '", ms_name, "' is NOT 'msms.txt'!"))
    }else{
      shinyBS::closeAlert(session, "a3")
    }
  })
  ### 4 summary.txt
  observeEvent(input$sum_file, {
    sum_file <- input$sum_file$datapath
    sum_name <- input$sum_file$name
    if(!is.null(sum_file) && sum_name!="summary.txt"){
      createAlert(session, anchorId = "alert4", alertId = "a4", style = "danger", content =  paste0("The uploaded file '", sum_name, "' is NOT 'summary.txt'!"))
    }else{
      shinyBS::closeAlert(session, "a4")
    }
  })
  ### 5 doe.txt
  observeEvent(input$doe_file, {
    if(isTruthy(input$doe_file)){
      doe_file <- input$doe_file$datapath
      doe <- read.table(doe_file, sep = input$sep, header = T)
      
      colnames(doe) <- gsub(" ", ".", colnames(doe))
      colnames(doe) <- tolower(colnames(doe))
      doe$type <- tolower(doe$type)
      # modify the sample.id if necessary
      doe$sample.id <- gsub(" ", ".", doe$sample.id)
      doe$sample.id <- gsub("-", ".", doe$sample.id)
      
      # check required columns in doe
      requiredColList <- c("sample.id", "run.order","type", "sample.type")
      colnames(doe) <- tolower(colnames(doe))
      doe.factors <- colnames(doe)
      noMatchCol <- 0
      notMatchCol_list <- ""
      for(i in requiredColList){
        if(!(i %in% doe.factors)){
          notMatchCol_list <- c(notMatchCol_list, i)
          
        }else{
          noMatchCol <- noMatchCol+1
        }
      }
      notMatchCol_list <- notMatchCol_list[-1]
      if(length(notMatchCol_list)!=0){
        createAlert(session, "alert5a", alertId = "a5a",
                    style = "danger",
                    content = paste0(notMatchCol_list," column is NOT in DOE! "),
                    append = TRUE)
      }else{
        if(noMatchCol==4){
          shinyBS::closeAlert(session, "a5a")
        }
      }

      # check if samples listed in DOE are existing in proteinGroups file
      if(isTruthy(input$pg_file)){
        shinyBS::closeAlert(session, "a10")
        pg_file <- input$pg_file$datapath
        pg_data <- read.table(pg_file, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
        pg_data_colname <- colnames(pg_data)
        pg_data_intensity_col <- pg_data %>% 
          select(starts_with("Intensity."))
        pg_data_intensity_colname <- colnames(pg_data_intensity_col)
        output$pg_data_intensity_colname <- renderText({
          pg_data_intensity_colname
        })
        
        doe_intensity_sampleid <- paste0("Intensity.", sort(doe$sample.id))
        doe_sample_count = length(doe$sample.id)
        output$doe_sample_name <- renderText({
          doe_intensity_sampleid
        })
        
        # 1. if DOE has more samples than proteinGroups
        notMatchSample_list1 <- ""
        for(i in doe_intensity_sampleid){
          if(!(i %in% pg_data_intensity_colname)){
            sample = gsub("Intensity.","",i)
            #notMatchSample_list1 <- c(notMatchSample_list1, i)
            notMatchSample_list1 <- c(notMatchSample_list1, sample)
          }
        }
        notMatchSample_list1 <- notMatchSample_list1[-1]
        if(length(notMatchSample_list1)!=0){
          #notMatchSample_list1 < gsub("Intensity.", "", notMatchSample_list1)
          createAlert(session, "alert5b", alertId = "a5b",
                      style = "warning",
                      content = paste0(notMatchSample_list1, " in DOE is not found in proteinGroups data file! "),
                      append = TRUE)
        }else{
          shinyBS::closeAlert(session, "a5b")
        }
        
        # 2. if DOE has less samples than proteinGroups
        notMatchSample_list2 <- ""
        for(i in pg_data_intensity_colname){
          if(!(i %in% doe_intensity_sampleid)){
            sample = gsub("Intensity.","",i)
            #notMatchSample_list2 <- c(notMatchSample_list2, i)
            notMatchSample_list2 <- c(notMatchSample_list2, sample)
          }
        }
        notMatchSample_list2 <- notMatchSample_list2[-1]
        if(length(notMatchSample_list2)!=0){
          #notMatchSample_list2 < gsub("Intensity.", "", notMatchSample_list2)
          createAlert(session, "alert5b", alertId = "a5b",
                      style = "warning",
                      content = paste0(notMatchSample_list2, " in proteinGroups data is not found in DOE file! "),
                      append = TRUE)
        }else{
          shinyBS::closeAlert(session, "a5c")
        }
        
        # 1+2 if DOE and proteinGroups are not consistent
        pg_data_sample_count = length(pg_data_intensity_colname)
        if(doe_sample_count!=pg_data_sample_count){
          createAlert(session, "alert5d", alertId = "a5d",
                      style = "danger",
                      content = "Data in DOE is not consistent with data in proteinGroups.txt")
        }else{
          shinyBS::closeAlert(session, "a5d")
        }
      }else{
        createAlert(session, "alert10", alertId = "a10",
                    style = "warning",
                    content = "The 'proteinGroups.txt' is required!")
      }
    }else{
      output$doe_descrp <- renderText({
        paste0("No DOE is provided in this analysis.")
      })
    }
  })
  
  ###################################################################################

  
  ####################################
  # quick check for input parameters #
  ####################################
  observeEvent(input$custom_spikeinRun,{
    if(isTruthy(input$pg_file)){
      pg_file <- input$pg_file$datapath
      pg <- read.table(pg_file, header = T, sep = "\t",
                       stringsAsFactors = F, quote = "\"",
                       comment.char = "")
      if(!custom_spike %in% pg$Protein.IDs){
        createAlert(session, "alert11", alertId = "a11",
                    style = "danger",
                    content = "The Protein.ID(s) is not found in the data!")
      }else{
        shinyBS::closeAlert(session, "a11")
      }
      
    }
  })
  
  parameter <- reactive({
    ### LFQRun
    LFQRun <- input$LFQRun
    if(isTruthy(input$pg_file)){
      pg_file <- input$pg_file$datapath
      pg_data <- read.table(pg_file, header = T, sep = "\t",
                            stringsAsFactors = F, quote = "\"",
                            comment.char = "")
      pg_data_lfq_col <- pg_data %>% 
        select(starts_with("LFQ"))
      if(length(colnames(pg_data_lfq_col))>0){
        LFQRun <- "Y"
      }else{
        LFQRun <- "N"
        output$lfq_descrp <- renderText({
          paste0("No LFQ is performed in this data.")
        })
      }
    }
    ### spikeinRun
    spikeinRun <- input$spikeinRun
    if(isTruthy(input$pg_file)){
      pg_file <- input$pg_file$datapath
      pg <- read.table(pg_file, header = T, sep = "\t",
                            stringsAsFactors = F, quote = "\"",
                            comment.char = "")
      if(length(pg$Protein.IDs[pg$Protein.IDs %in% "P00000"])==0){
        spikeinRun <- "N"
        output$spike_descrp <- renderText({
          paste0("No spike-in is added in this experiment.")
        })
      }else{
        spikeinRun <- "Y"
        output$spike_descrp <- renderText({
          paste0("")
        })
      }
    }
    
    ### prefix
    prefix <- input$prefix
    
    ### exp_title
    if(isTruthy(input$exp_title)){
      exp_title <- input$exp_title
      exp_title <- gsub(" ","_", exp_title)
    }else{
      exp_title <- "QC-MQ"
    }
  })
  
  
  ###################################################################################
  
  
  ###########################
  # response upon the files #
  ###########################
  
  ### 1 proteinGroups.txt
  pg_response <- reactive({
    pg_file <- input$pg_file$datapath
    pg <- read.table(pg_file, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
    
    parsed_pg <- process_pg(pg)
    
    # check spikein
    if(length(pg$Protein.IDs[pg$Protein.IDs %in% "P00000"])==0){
      spikeinRun <- "N"
      output$spike_descrp <- renderText({
        paste0("No spike-in is added in this experiment.")
      })
    }else{
      spikeinRun <- "Y"
      output$spike_descrp <- renderText({
        paste0("")
      })
    }
    
    # check LFQ
    pg_lfq_col <- parsed_pg %>% 
      select(starts_with("LFQ"))
    if(length(colnames(pg_lfq_col))>0){
      LFQRun <- "Y"
    }else{
      LFQRun <- "N"
    }
    ################################################################################################
    # datasets for download
    ## raw
    table_log2_raw <- getLog2Table(parsed_pg, "raw")
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    table_contam_log2_raw <- getLog2ContamTable(parsed_pg, "raw")
    
    output$dtable_log2_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_log2_raw_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_log2_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dtable_log2_clean_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_log2_raw_clean_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_log2_clean_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dtable_contam_log2_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_log2_raw_contaminant_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_contam_log2_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    ## lfq
    if(LFQRun=="Y"){
      table_log2_lfq <- getLog2Table(parsed_pg, "lfq")
      table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
      table_contam_log2_lfq <- getLog2ContamTable(parsed_pg, "lfq")
      
      output$dtable_log2_lfq <- downloadHandler(
        filename = function(){
          paste0(input$prefix,"_log2_lfq_data_",Sys.Date(),".tsv")
        },
        content = function(file){
          write.table(table_log2_lfq, file, sep = "\t", row.names = FALSE)
        }
      )
      
      output$dtable_log2_clean_lfq <- downloadHandler(
        filename = function(){
          paste0(input$prefix,"_log2_lfq_clean_data_",Sys.Date(),".tsv")
        },
        content = function(file){
          write.table(table_log2_clean_lfq, file, sep = "\t", row.names = FALSE)
        }
      )
      
      output$dtable_contam_log2_lfq <- downloadHandler(
        filename = function(){
          paste0(input$prefix,"_log2_lfq_contaminant_data_",Sys.Date(),".tsv")
        },
        content = function(file){
          write.table(table_contam_log2_lfq, file, sep = "\t", row.names = FALSE)
        }
      )
      
    }
    

    ##########################################################################################################
    # expCheck - spikein
    
    if(length(parsed_pg$Protein.IDs[parsed_pg$Protein.IDs %in% "P00000"])!=0){
      ## raw
      spike_dat_raw <- getSpikeInfo(parsed_pg, "raw")
      spike_raw_log2_mean <- mean(spike_dat_raw$log2.intensity)
      spike_raw_log2_sd <- sd(spike_dat_raw$log2.intensity)
      cutline_below_raw = spike_raw_log2_mean - 2*spike_raw_log2_sd
      cutline_above_raw = spike_raw_log2_mean + 2*spike_raw_log2_sd
      below_raw <- getSpikeBelow2sdSample(spike_dat_raw)
      above_raw <- getSpikeAbove2sdSample(spike_dat_raw)
      spike_dat_raw <- spike_dat_raw %>% 
        mutate(outlier = ifelse(sample %in% below_raw$sample | sample %in% above_raw$sample, "#FF6666", "#00CCCC"))
      
      output$spike_plot_raw <- renderPlotly({
        n = nrow(spike_dat_raw)
        plot_ly(spike_dat_raw, x=~sample, y=~log2.intensity, type = "scatter", color =~ I(outlier), name="Spikein") %>% 
          add_segments(y=cutline_above_raw, yend=cutline_above_raw, x=spike_dat_raw$sample[1], xend=spike_dat_raw$sample[n], color=I("gray"), name="Mean+2SD") %>% 
          add_segments(y=spike_raw_log2_mean, yend=spike_raw_log2_mean, x=spike_dat_raw$sample[1], xend=spike_dat_raw$sample[n], color=I("black"), name="Mean") %>% 
          add_segments(y=cutline_below_raw, yend=cutline_below_raw, x=spike_dat_raw$sample[1], xend=spike_dat_raw$sample[n], color=I("gray"), name="Mean-2SD") 
      })
      
      output$spikelow_raw <- renderTable({
        below_raw %>% select(sample, log2.intensity)
      })
      
      output$spikehigh_raw <- renderTable({
        above_raw %>% select(sample, log2.intensity)
      })
      
      ## lfq
      if(LFQRun=="Y"){
        spike_dat_lfq <- getSpikeInfo(parsed_pg, "lfq")
        spike_lfq_log2_mean <- mean(spike_dat_lfq$log2.intensity)
        spike_lfq_log2_sd <- sd(spike_dat_lfq$log2.intensity)
        cutline_below_lfq = spike_lfq_log2_mean - 2*spike_lfq_log2_sd
        cutline_above_lfq = spike_lfq_log2_mean + 2*spike_lfq_log2_sd
        below_lfq <- getSpikeBelow2sdSample(spike_dat_lfq)
        above_lfq <- getSpikeAbove2sdSample(spike_dat_lfq)
        spike_dat_lfq <- spike_dat_lfq %>% 
          mutate(outlier = ifelse(sample %in% below_lfq$sample | sample %in% above_lfq$sample, "#FF6666", "#00CCCC"))
        
        output$spike_plot_lfq <- renderPlotly({
          n = nrow(spike_dat_lfq)
          plot_ly(spike_dat_lfq, x=~sample, y=~log2.intensity, type = "scatter", color =~ I(outlier), name="Spikein") %>% 
            add_segments(y=cutline_above_lfq, yend=cutline_above_lfq, x=spike_dat_lfq$sample[1], xend=spike_dat_lfq$sample[n], color=I("gray"), name="Mean+2SD") %>% 
            add_segments(y=spike_lfq_log2_mean, yend=spike_lfq_log2_mean, x=spike_dat_lfq$sample[1], xend=spike_dat_lfq$sample[n], color=I("black"), name="Mean") %>% 
            add_segments(y=cutline_below_lfq, yend=cutline_below_lfq, x=spike_dat_lfq$sample[1], xend=spike_dat_lfq$sample[n], color=I("gray"), name="Mean-2SD") 
        })
        
        output$spikelow_lfq <- renderTable({
          below_lfq %>% select(sample, log2.intensity)
        })
        
        output$spikehigh_lfq <- renderTable({
          above_lfq %>% select(sample, log2.intensity)
        })
      }
      
    }

    
    

    ##################################################################################################################################################################################
    # contaminants
    
    ## raw
    raw_con_info <- getContamInfo(parsed_pg, "raw")
    raw_con_t1 <- getContamTable1(raw_con_info)
    
    raw_con_t1g <- gather(raw_con_t1, key = "Contaminant.type", value = "Total.intensity", -sample)
    
    output$contam_raw <- renderPlotly(
      ggplotly(
        ggplot(raw_con_t1g, aes(x=sample, y=Total.intensity, color=Contaminant.type))+
          geom_point()+
          labs(y="Total intensity (1og10 scale)")+
          theme_minimal()+
          theme(axis.text.x = element_text(angle=90))+
          scale_y_log10()
      )
    )
    
    ## lfq
    if(LFQRun=="Y"){
      lfq_con_info <- getContamInfo(parsed_pg, "lfq")
      lfq_con_t1 <- getContamTable1(lfq_con_info)
      
      lfq_con_t1g <- gather(lfq_con_t1, key = "Contaminant.type", value = "Total.intensity", -sample)
      
      output$contam_lfq <- renderPlotly(
        ggplotly(
          ggplot(lfq_con_t1g, aes(x=sample, y=Total.intensity, color=Contaminant.type))+
            geom_point()+
            labs(y="Total intensity (1og10 scale)")+
            theme_minimal()+
            theme(axis.text.x = element_text(angle=90))+
            scale_y_log10()
        )
      )
    }
    
    ###############################################################################################################
    # protein counts
    ## raw
    raw_prog_count <- getProteinCount(parsed_pg, "raw")
    output$pg_count_raw <- renderPlotly(
      plot_ly(raw_prog_count, x=~sample, y=~count, type = "scatter") %>% 
        layout(yaxis=list(title="Number of Protein Group"))
    )
    
    ## lfq
    if(LFQRun=="Y"){
      lfq_prog_count <- getProteinCount(parsed_pg, "lfq")
      output$pg_count_lfq <- renderPlotly(
        plot_ly(lfq_prog_count, x=~sample, y=~count, type = "scatter", color=I("#FF6666")) %>% 
          layout(yaxis=list(title="Number of Protein Group"))
      )
    }

    ###############################################################################################################
    # protein intensity
    
    # total by individual
    ## raw
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    prog_int_sum_raw <- getTotalIntensity(table_log2_clean_raw)
    
    output$pg_int_raw <- renderPlotly(
      plot_ly(prog_int_sum_raw, x=~sample, y=~Total.Log2.Intensity, type="scatter")
    )
    
    ## lfq
    if(LFQRun=="Y"){
      table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
      prog_int_sum_lfq <- getTotalIntensity(table_log2_clean_lfq)
      
      output$pg_int_lfq <- renderPlotly(
        plot_ly(prog_int_sum_lfq, x=~sample, y=~Total.Log2.Intensity, type="scatter", color=I("#FF6666"))
      )
    }

    ###############################################################################################################
    # protein intensity
    
    # violin distribution
    ## raw
    pg_ints_raw <- getIntensity(table_log2_clean_raw)
    pg_ints_raw_melt <- melt(pg_ints_raw, id.vars="sample")
    pg_ints_raw_melt$value = as.numeric(pg_ints_raw_melt$value)
    
    output$pg_int_vio_raw <- renderPlotly(
      plot_ly(pg_ints_raw_melt, x=~sample, y=~value, type="violin", split =~sample, color = I("#3366CC"), showlegend =FALSE) %>% 
        layout(yaxis=list(title="Log2 Intensity"), xaxis=list(title="Sample"))
    )
    
    ## lfq
    if(LFQRun=="Y"){
      pg_ints_lfq <- getIntensity(table_log2_clean_lfq)
      pg_ints_lfq_melt <- melt(pg_ints_lfq, id.vars="sample")
      pg_ints_lfq_melt$value = as.numeric(pg_ints_lfq_melt$value)
      
      output$pg_int_vio_lfq <- renderPlotly(
        plot_ly(pg_ints_lfq_melt, x=~sample, y=~value, type="violin", split =~sample, color = I("#FF6666"), showlegend =FALSE) %>% 
          layout(yaxis=list(title="Log2 Intensity"), xaxis=list(title="Sample"))
      )
    }
    
    ###############################################################################################################
    # protein intensity
    
    # Raw vs LFQ overll
    # get the melt data from raw and lfq: pg_ints_raw/lfq_melt
    if(LFQRun=="Y"){
      melt_table <- data.frame(raw=pg_ints_raw_melt$value,
                               lfq=pg_ints_lfq_melt$value)
      
      output$pg_int_rawlfq <- renderPlotly(
        plot_ly(melt_table, x=~raw, y=~lfq, type="scatter", marker=list(opacity=0.3)) %>% 
          layout(xaxis=list(title="Raw Intensity"), yaxis=list(title="LFQ Intensity"))
      )
    }
    
    
    ###############################################################################################################
    # CV
    ## raw
    cv_raw <- ddply(pg_ints_raw_melt, .(variable), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_raw <- renderPlotly(
      ggplotly(
        ggplot(cv_raw, aes(x=cv))+
          geom_density(alpha=0.5, fill="#3366CC")+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    
    ## lfq
    if(LFQRun=="Y"){
      cv_lfq <- ddply(pg_ints_lfq_melt, .(variable), summarise, cv=sd(value)/mean(value)*100)
      
      output$pg_int_cv_lfq <- renderPlotly(
        ggplotly(
          ggplot(cv_lfq, aes(x=cv))+
            geom_density(alpha=0.5, fill="#FF6666")+
            theme_minimal()+
            labs(x="Coefficient of Variation (%)", y="Density")
        )
      )
    }
    
    
    ###############################################################################################################
    # top20
    # get the mean intensity of all the proteins from table_log2_clean_raw/lfq
    
    ## raw
    melt_table_name_raw <- getMeltwithNames(table_log2_clean_raw)
    annotated_mean_table_raw <- countMeanbyProtein(melt_table_name_raw)
    top20.raw <- annotated_mean_table_raw[1:20,]
    output$top20_raw <- renderTable(
      top20.raw, striped = TRUE
    )
    
    output$dannotated_mean_table_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_mean_log2_raw_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(annotated_mean_table_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    ## lfq
    if(LFQRun=="Y"){
      melt_table_name_lfq <- getMeltwithNames(table_log2_clean_lfq)
      annotated_mean_table_lfq <- countMeanbyProtein(melt_table_name_lfq)
      top20.lfq <- annotated_mean_table_lfq[1:20,]
      output$top20_lfq <- renderTable(
        top20.lfq, striped = TRUE
      )
      
      output$dannotated_mean_table_lfq <- downloadHandler(
        filename = function(){
          paste0(input$prefix,"_mean_log2_lfq_",Sys.Date(),".tsv")
        },
        content = function(file){
          write.table(annotated_mean_table_lfq, file, sep = "\t", row.names = FALSE)
        }
      )
    }
    
    
    ###############################################################################################################
    # pca
    ## raw
    pca_table_raw <- getPCAtable(table_log2_clean_raw)
    pc12_raw <- getPC12Percentage(table_log2_clean_raw)
    
    output$pca_raw <- renderPlotly(
      #ggplotly(
        #ggplot(pca_table_raw, aes(x=PC1, y=PC2, color=sample) )+ 
        #  geom_point(alpha=0.5, size=3)+
        #  theme_bw()+
        #  labs(x=paste0("PC1 (", pc12_raw[1], "%)"),
        #       y=paste0("PC2 (", pc12_raw[2], "%)"))
      #)
      plot_ly(pca_table_raw, x=~PC1, y=~PC2, type = "scatter", color=~sample, text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE),
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
    )
    ## lfq
    if(LFQRun=="Y"){
      pca_table_lfq <- getPCAtable(table_log2_clean_lfq)
      pc12_lfq <- getPC12Percentage(table_log2_clean_lfq)
      
      
      output$pca_lfq <- renderPlotly(
        #ggplotly(
        #  ggplot(pca_table_lfq, aes(x=PC1, y=PC2, color=sample) )+ 
        #    geom_point(alpha=0.5, size=3)+
        #    theme_bw()+
        #    labs(x=paste0("PC1 (", pc12_lfq[1], "%)"),
        #         y=paste0("PC2 (", pc12_lfq[2], "%)"))
        #)
        plot_ly(pca_table_lfq, x=~PC1, y=~PC2, type = "scatter", color=~sample, text=~sample, marker=list(size=16, opacity=0.5)) %>% 
          layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE),
                 yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
        
        
        
      )
    }
    
    
    ###############################################################################################################
    # end
    
  })
  
  ### 2 peptides.txt
  pt_response <- reactive({
    pt_file <- input$pt_file$datapath
    pt <- read.table(pt_file, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
    
    # get the average of Cysteine
    cys_avg = sum(pt$C.Count)/nrow(pt)
    # get the average of missed cleavage
    ms_clvg_avg = sum(pt$Missed.cleavage)/nrow(pt)
    # make checkpoint table
    avg_table = data.frame("Feature" = c("Cysteine","Missed Cleavage"),
                           "Percentage" = c(paste0(format(round(cys_avg*100, 2), nsmall = 2), " %"), paste0(format(round(ms_clvg_avg*100, 2), nsmall = 2), " %")))
    output$expck_table <- renderTable(
      avg_table, striped = T, bordered = T
    )
    
    # get peptide count by individual sample
    
    pep_count <- getPepCount(pt)
    
    output$pt_count <- renderPlotly(
      plot_ly(pep_count, x=~sample, y=~count, type = "scatter") %>% 
        layout(yaxis=list(title="Number of Peptide"))
    )
    
    # get peptide total log2 intensity
    
    pep_int_sum <- getTotalIntensityPep(pt)
    
    output$pt_int <- renderPlotly(
      plot_ly(pep_int_sum, x=~sample, y=~Total.Log2.Intensity, type="scatter", color = I("#3366CC"))
    )
    
  })
  
  ### 3 msms.txt
  ### 4 summary.txt
  ms_sum_response <- reactive({
    ms_file <- input$ms_file$datapath
    ms <- read.table(ms_file, header = T, sep = "\t",
                     stringsAsFactors = F, quote = "\"",
                     comment.char = "")
    sum_file <- input$sum_file$datapath
    sumtxt <- read.table(sum_file, header = T, sep = "\t",
                         stringsAsFactors = F, quote = "\"",
                         comment.char = "")
    
    ms_clvg_table <- ms %>% 
      select(Raw.file, Missed.cleavages) %>% 
      left_join(sumtxt %>% select(Raw.file, Experiment)) %>% 
      dplyr::count(Experiment, Missed.cleavages) %>% 
      spread(Missed.cleavages, n) %>% 
      mutate(sum = `0`+`1`+`2`) %>% 
      mutate(`percentage.0` = `0`/sum*100,
             `percentage.1` = `1`/sum*100,
             `percentage.2` = `2`/sum*100)
    
    output$msclvg_persample <- renderPlotly(
      plot_ly(ms_clvg_table, x= ~Experiment, y=~`percentage.0`, type = 'bar', name = "percentage.0") %>% 
        add_trace(y = ~`percentage.1`, name = "percentage.1") %>%
        add_trace(y = ~`percentage.2`, name = "percentage.2") %>%
        layout(yaxis=list(title = 'Percentage (%)'), barmode = 'group')
    )
    
  })
  
  ### 5 doe.txt
  doe_response <- reactive({
    doe_file <- input$doe_file$datapath
    doe <- read.table(doe_file, sep = input$sep, header = T)
    
    colnames(doe) <- gsub(" ", ".", colnames(doe))
    colnames(doe) <- tolower(colnames(doe))
    doe$type <- tolower(doe$type)
    
    doe$sample.id <- gsub(" ", ".", doe$sample.id)
    doe$sample.id <- gsub("-", ".", doe$sample.id)
    
    
    
    ################################################################################################
    # read pg as well
    pg <- read.table(input$pg_file$datapath, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
    parsed_pg <- process_pg(pg)
    
    # check LFQ
    pg_lfq_col <- parsed_pg %>% 
      select(starts_with("LFQ"))
    if(length(colnames(pg_lfq_col))>0){
      LFQRun <- "Y"
    }else{
      LFQRun <- "N"
    }
    ################################################################################################
    # data structure
    int_col <- pg %>% 
      select(starts_with("Intensity."))
    pg_sample_count <- ncol(int_col)
    
    output$doe_descrp <- renderText({
      paste0("There are ", pg_sample_count, " samples in the data")
    })
    
    sampletype_table <- dplyr::count(doe, sample.type)
    
    output$datastructure_table <- renderTable(
      sampletype_table, striped = T, bordered = T
    )
    
    output$datastructure_plot <- renderPlotly(
      plot_ly(sampletype_table, x= ~sample.type, y=~n, type = 'bar', color = ~sample.type) %>% 
        layout(yaxis=list(title = 'Count'))
    )
    
    ################################################################################################
    # contaminant type with sample type #pg
    ## raw
    raw_con_info <- getContamInfo(parsed_pg, "raw")
    raw_con_t2 <- getContamTable2(raw_con_info, doe)
    
    output$contam_raw_st <- renderPlotly(
      plot_ly(raw_con_t2, y=~log2(intensity), x=~sample.type, color=~kind_contam, type = "box") %>% 
        layout(yaxis=list(title="Log2 Raw Intensity" ),boxmode = "group")
    )
    
    ## lfq
    if(LFQRun=="Y"){
      lfq_con_info <- getContamInfo(parsed_pg, "lfq")
      lfq_con_t2 <- getContamTable2(lfq_con_info, doe)
      
      output$contam_lfq_st <- renderPlotly(
        plot_ly(lfq_con_t2, y=~log2(intensity), x=~sample.type, color=~kind_contam, type = "box") %>% 
          layout(yaxis=list(title="Log2 LFQ Intensity" ),boxmode = "group")
      )
    }
    
    ################################################################################################
    # protein counts 
    # sample type
    ## raw
    raw_prog_count_doe <- getProteinCountDoe(parsed_pg, "raw", doe)
    output$pg_count_st_raw <- renderPlotly(
      plot_ly(raw_prog_count_doe, x=~sample.type, y=~count, type = "box") %>% 
        layout(yaxis=list(title="Number of Protein Groups"))
    )
    
    ## lfq
    if(LFQRun=="Y"){
      lfq_prog_count_doe <- getProteinCountDoe(parsed_pg, "lfq", doe)
      output$pg_count_st_lfq <- renderPlotly(
        plot_ly(lfq_prog_count_doe, x=~sample.type, y=~count, type = "box", color=I("#FF6666")) %>% 
          layout(yaxis=list(title="Number of Protein Groups"))
      )
    }
    
    
    # type cdf plot
    ## raw
    per_data_melt_raw <- getExpressedSamplePercentPerProteinEachType(parsed_pg, "raw", doe)
    output$pg_count_t_raw <- renderPlot(
      ggplot(per_data_melt_raw, aes(x=Value, color=Type))+
        stat_ecdf(geom = "step")+
        scale_x_reverse()+
        labs(x="Percent of samples expressed per protein (%)", y="Proportion of expressed proteins")+
        theme_minimal()
    )
    
    ## lfq
    if(LFQRun=="Y"){
      per_data_melt_lfq <- getExpressedSamplePercentPerProteinEachType(parsed_pg, "lfq", doe)
      output$pg_count_t_lfq <- renderPlot(
        ggplot(per_data_melt_lfq, aes(x=Value, color=Type))+
          stat_ecdf(geom = "step")+
          scale_x_reverse()+
          labs(x="Percent of samples expressed per protein (%)", y="Proportion of expressed proteins")+
          theme_minimal()
      )
    }
    
    
    ################################################################################################
    # peptide counts
    # get peptide count by sample.type and type
    
    if(isTruthy(input$pt_file)){
      pt <- read.table(input$pt_file$datapath, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
      
      pep_count_doe <- getPepCountDoe(pt, doe)
      
      output$pt_count_st <- renderPlotly(
        plot_ly(pep_count_doe, x=~sample.type, y=~count, type = "box") %>% 
          layout(yaxis=list(title="Number of Peptide"), xaxis=list(title="Sample Type"))
      )
      
      output$pt_count_t <- renderPlotly(
        plot_ly(pep_count_doe, x=~type, y=~count, type = "box", color=~type) %>% 
          layout(yaxis=list(title="Number of Peptide"), xaxis=list(title="Type"))
      )
    }
    
    ################################################################################################
    # protein intensity - mean proportion plot
    ## raw
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    pg_ints_raw <- getIntensity(table_log2_clean_raw)
    pg_ints_raw_melt <- melt(pg_ints_raw, id.vars="sample")
    pg_ints_raw_melt$value = as.numeric(pg_ints_raw_melt$value)
    
    mp_raw <- getMedianIntPlotData(pg_ints_raw_melt, doe)
    
    output$pg_int_mpp_raw <- renderPlotly(
      plot_ly(mp_raw, type = "scatter", mode = "markers") %>% 
        add_trace(y = ~All, name = 'All sample', marker = list(opacity=0.2, color="#3366CC")) %>% 
        layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
               xaxis=list(title="Protein Group"))
    )
    
    output$pg_int_mpp_t_raw <- renderPlotly(
      plot_ly(mp_raw, type = "scatter", mode = "markers") %>% 
        add_trace(y = ~Control, name = 'Control', marker = list(opacity=0.2)) %>%
        add_trace(y = ~Sample, name = 'Sample', marker = list(opacity=0.2)) %>% 
        layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
               xaxis=list(title="Protein Group"))
    )
    ## lfq
    if(LFQRun=="Y"){
      table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
      pg_ints_lfq <- getIntensity(table_log2_clean_lfq)
      pg_ints_lfq_melt <- melt(pg_ints_lfq, id.vars="sample")
      pg_ints_lfq_melt$value = as.numeric(pg_ints_lfq_melt$value)
      
      mp_lfq <- getMedianIntPlotData(pg_ints_lfq_melt, doe)
      
      output$pg_int_mpp_lfq <- renderPlotly(
        plot_ly(mp_lfq, type = "scatter", mode = "markers") %>% 
          add_trace(y = ~All, name = 'All sample', marker = list(opacity=0.2, color="#FF6666")) %>% 
          layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
                 xaxis=list(title="Protein Group"))
      )
      
      output$pg_int_mpp_t_lfq <- renderPlotly(
        plot_ly(mp_lfq, type = "scatter", mode = "markers") %>% 
          add_trace(y = ~Control, name = 'Control', marker = list(opacity=0.2)) %>%
          add_trace(y = ~Sample, name = 'Sample', marker = list(opacity=0.2)) %>% 
          layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
                 xaxis=list(title="Protein Group"))
      )
    }
    
    
    ################################################################################################
    # protein intensity - cdf plot and violin plot
    ## raw
    pg_ints_raw_melt_doe = pg_ints_raw_melt %>% 
      left_join(doe, by=c("sample"="sample.id"))
    
    output$pg_int_cdf_st_raw <- renderPlot(
      ggplot(pg_ints_raw_melt_doe, aes(x=value+0.000001, color=sample))+
        stat_ecdf()+
        labs(x = "Log2 Intensities by Protein", 
             y = "Fraction of Library") + 
        theme_minimal() +
        facet_wrap(~sample.type, ncol = 2)
      
    )
    
    output$pg_int_vio_raw_doe <- renderPlotly(
      ggplotly(
        ggplot(pg_ints_raw_melt_doe, aes(x=reorder(sample,run.order), y=value, color=sample.type, fill=sample.type))+
          geom_violin()+
          theme_minimal()+
          labs(x="Sample by Run order", y="Log2 Intensity")+
          theme(axis.text.x = element_text(angle = 90))
      )
    )
    
    ## lfq
    if(LFQRun=="Y"){
      pg_ints_lfq_melt_doe = pg_ints_lfq_melt %>% 
        left_join(doe, by=c("sample"="sample.id"))
      
      output$pg_int_cdf_st_lfq <- renderPlot(
        ggplot(pg_ints_lfq_melt_doe, aes(x=value+0.000001, color=sample))+
          stat_ecdf()+
          labs(x = "Log2 Intensities by Protein", 
               y = "Fraction of Library") + 
          theme_minimal() +
          facet_wrap(~sample.type, ncol = 2)
      )
      
      output$pg_int_vio_lfq_doe <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_lfq_melt_doe, aes(x=reorder(sample,run.order), y=value, color=sample.type, fill=sample.type))+
            geom_violin()+
            theme_minimal()+
            labs(x="Sample by Run order", y="Log2 Intensity")+
            theme(axis.text.x = element_text(angle = 90))
        )
      ) 
    }
     
    
    ################################################################################################
    # protein intensity - Raw vs LFQ sample.type
    # get the melt data from raw and lfq: pg_ints_raw/lfq_melt_doe
    if(LFQRun=="Y"){
      melt_table_doe <- data.frame(sample=pg_ints_raw_melt_doe$sample,
                                   type=pg_ints_raw_melt_doe$type,
                                   sample.type=pg_ints_raw_melt_doe$sample.type,
                                   raw=pg_ints_raw_melt_doe$value,
                                   lfq=pg_ints_lfq_melt_doe$value)
      
      st_num <- length(unique(melt_table_doe$sample.type))
      if(st_num<=3){
        n=2
      }else{
        n=3
      }
      n <- as.numeric(n)
      
      output$pg_int_st_rawlfq <- renderPlotly(
        
        melt_table_doe %>% 
          group_by(sample.type) %>% 
          group_map(~ plot_ly(data = ., x=~raw, y=~lfq, color=~sample.type, type="scatter", marker=list(opacity=0.3)), keep = TRUE) %>% 
          subplot(nrows = n, shareX = TRUE, shareY = TRUE) %>% 
          layout(xaxis=list(title="Raw Intensity"), yaxis=list(title="LFQ Intensity"))
      )
    }
    
    ################################################################################################
    # CV with doe, sample type
    ## raw
    cv_raw_st <- ddply(pg_ints_raw_melt_doe, .(variable, sample.type), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_st_raw <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_st, aes(x=cv, fill=sample.type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    output$pg_int_cv_st_raw_split <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_st, aes(x=cv, fill=sample.type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")+
          facet_grid(~sample.type)
      )
    )
    
    ## lfq
    if(LFQRun=="Y"){
      cv_lfq_st <- ddply(pg_ints_lfq_melt_doe, .(variable, sample.type), summarise, cv=sd(value)/mean(value)*100)
      
      output$pg_int_cv_st_lfq <- renderPlotly(
        ggplotly(
          ggplot(cv_lfq_st, aes(x=cv, fill=sample.type))+
            geom_density(alpha=0.5)+
            theme_minimal()+
            labs(x="Coefficient of Variation (%)", y="Density")
        )
      )
      output$pg_int_cv_st_lfq_split <- renderPlotly(
        ggplotly(
          ggplot(cv_lfq_st, aes(x=cv, fill=sample.type))+
            geom_density(alpha=0.5)+
            theme_minimal()+
            labs(x="Coefficient of Variation (%)", y="Density")+
            facet_grid(~sample.type)
        )
      )
    }
    
    
    # CV with doe, type
    ## raw
    cv_raw_t <- ddply(pg_ints_raw_melt_doe, .(variable, type), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_t_raw <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_t, aes(x=cv, fill=type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    
    output$pg_int_cv_t_raw_split <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_t, aes(x=cv, fill=type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")+
          facet_grid(~type)
      )
    )
    
    ## lfq
    if(LFQRun=="Y"){
      cv_lfq_t <- ddply(pg_ints_lfq_melt_doe, .(variable, type), summarise, cv=sd(value)/mean(value)*100)
      
      output$pg_int_cv_t_lfq <- renderPlotly(
        ggplotly(
          ggplot(cv_lfq_t, aes(x=cv, fill=type))+
            geom_density(alpha=0.5)+
            theme_minimal()+
            labs(x="Coefficient of Variation (%)", y="Density")
        )
      )
      
      output$pg_int_cv_t_lfq_split <- renderPlotly(
        ggplotly(
          ggplot(cv_lfq_t, aes(x=cv, fill=type))+
            geom_density(alpha=0.5)+
            theme_minimal()+
            labs(x="Coefficient of Variation (%)", y="Density")+
            facet_grid(~type)
        )
      )
    }
    
    
    ################################################################################################
    # top20 sample.type
    ## raw
    melt_table_name_raw <- getMeltwithNames(table_log2_clean_raw)
    melt_table_name_raw_doe <- melt_table_name_raw %>% 
      left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
    
 
    wide_sampletype_table_raw <- countMeanbyProteinSampletypeDoe(melt_table_name_raw_doe)
    name_col_st_raw <- wide_sampletype_table_raw[,c(1:3)]
    measure_col_st_raw <- wide_sampletype_table_raw[,c(4:ncol(wide_sampletype_table_raw))]
    
    output$sampletype_raw <- renderUI({
      st_list <- sort(unique(doe$sample.type))
      selectInput("st_selected_raw", h4("Choose Sample Type"),
                  choices = st_list, selected = 1)
    })
    
    
    output$dwide_sampletype_sort_table_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_meanBySampletype_log2_raw_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(wide_sampletype_table_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    ## lfq
    if(LFQRun=="Y"){
      melt_table_name_lfq <- getMeltwithNames(table_log2_clean_lfq)
      melt_table_name_lfq_doe <- melt_table_name_lfq %>% 
        left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
      
      wide_sampletype_table_lfq <- countMeanbyProteinSampletypeDoe(melt_table_name_lfq_doe)
      name_col_st_lfq <- wide_sampletype_table_lfq[,c(1:3)]
      measure_col_st_lfq <- wide_sampletype_table_lfq[,c(4:ncol(wide_sampletype_table_lfq))]
      
      output$sampletype_lfq <- renderUI({
        st_list <- sort(unique(doe$sample.type))
        selectInput("st_selected_lfq", h4("Choose Sample Type"),
                    choices = st_list, selected = 1)
      })
      
      
      output$dwide_sampletype_sort_table_lfq <- downloadHandler(
        filename = function(){
          paste0(input$prefix,"_meanBySampletype_log2_lfq_",Sys.Date(),".tsv")
        },
        content = function(file){
          write.table(wide_sampletype_sort_table_lfq, file, sep = "\t", row.names = FALSE)
        }
      )
    }
    

    # top20 type
    ## raw
    wide_type_table_raw <- countMeanbyProteinTypeDoe(melt_table_name_raw_doe)
    name_col_t_raw <- wide_type_table_raw[,c(1:3)]
    measure_col_t_raw <- wide_type_table_raw[,c(4:ncol(wide_type_table_raw))]
    
    output$type_raw <- renderUI({
      t_list <- sort(unique(doe$type))
      selectInput("t_selected_raw", h4("Choose Type"),
                  choices = t_list, selected = 1)
    })
    
    
    output$dwide_type_sort_table_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_meanByType_log2_raw_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(wide_type_table_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    ## lfq
    if(LFQRun=="Y"){
      
      wide_type_table_lfq <- countMeanbyProteinTypeDoe(melt_table_name_lfq_doe)
      name_col_t_lfq <- wide_type_table_lfq[,c(1:3)]
      measure_col_t_lfq <- wide_type_table_lfq[,c(4:ncol(wide_type_table_lfq))]
      
      output$type_lfq <- renderUI({
        t_list <- sort(unique(doe$type))
        selectInput("t_selected_lfq", h4("Choose Type"),
                    choices = t_list, selected = 1)
      })
      
      
      output$dwide_type_sort_table_lfq <- downloadHandler(
        filename = function(){
          paste0(input$prefix,"_meanByType_log2_lfq_",Sys.Date(),".tsv")
        },
        content = function(file){
          write.table(wide_type_table_lfq, file, sep = "\t", row.names = FALSE)
        }
      )
    }
    
    
    ################################################################################################
    # pca sample    sample.type, run.order
    ## raw
    pca_table_raw_doe <- getPCAtableDOE(table_log2_clean_raw, doe)
    pc12_raw <- getPC12Percentage(table_log2_clean_raw)
    
    output$pca_st_raw <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_raw_doe, aes(x=PC1, y=PC2, color=sample.type, shape=type) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_raw[1], "%)"),
      #         y=paste0("PC2 (", pc12_raw[2], "%)"))
      #)
      plot_ly(pca_table_raw_doe, x=~PC1, y=~PC2, type = "scatter", color=~sample.type, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE), 
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
    )
    
    output$pca_ro_raw <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_raw_doe, aes(x=PC1, y=PC2, color=run.order, shape=type) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_raw[1], "%)"),
      #         y=paste0("PC2 (", pc12_raw[2], "%)"))
      #)
      plot_ly(pca_table_raw_doe, x=~PC1, y=~PC2, type = "scatter", color=~run.order, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE), 
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
    )
    
    ## lfq
    if(LFQRun=="Y"){
      pca_table_lfq_doe <- getPCAtableDOE(table_log2_clean_lfq, doe)
      pc12_lfq <- getPC12Percentage(table_log2_clean_lfq)
      
      output$pca_st_lfq <- renderPlotly(
        #ggplotly(
        #  ggplot(pca_table_lfq_doe, aes(x=PC1, y=PC2, color=sample.type, shape=type) )+ 
        #    geom_point(alpha=0.5, size=3)+
        #    theme_bw()+
        #    labs(x=paste0("PC1 (", pc12_lfq[1], "%)"),
        #         y=paste0("PC2 (", pc12_lfq[2], "%)"))
        #)
        plot_ly(pca_table_lfq_doe, x=~PC1, y=~PC2, type = "scatter", color=~sample.type, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
          layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE), 
                 yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
      )
      
      output$pca_ro_lfq <- renderPlotly(
        #ggplotly(
        #  ggplot(pca_table_lfq_doe, aes(x=PC1, y=PC2, color=run.order, shape=type) )+ 
        #    geom_point(alpha=0.5, size=3)+
        #    theme_bw()+
        #    labs(x=paste0("PC1 (", pc12_lfq[1], "%)"),
        #         y=paste0("PC2 (", pc12_lfq[2], "%)"))
        #)
        plot_ly(pca_table_lfq_doe, x=~PC1, y=~PC2, type = "scatter", color=~run.order, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
          layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE), 
                 yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
      )
    }
    
    
    ################################################################################################
    # contribution of the variance
    ## raw 
    pvca_raw <- getPVCAresult(table_log2_clean_raw, doe)
    
    output$ctbv_raw <- renderPlotly(
      plot_ly(pvca_raw, x=~Variable, y=~Proportion.Variance, type="bar", color=~Variable) %>% 
        layout(yaxis=list(title="Weighted average proportion variance"), xaxis=list(title="Variables"))
    )
    
    ## lfq
    if(LFQRun=="Y"){
      pvca_lfq <- getPVCAresult(table_log2_clean_lfq, doe)
      
      output$ctbv_lfq <- renderPlotly(
        plot_ly(pvca_lfq, x=~Variable, y=~Proportion.Variance, type="bar", color=~Variable) %>% 
          layout(yaxis=list(title="Weighted average proportion variance"), xaxis=list(title="Variables"))
      )
    }
    
    ################################################################################################
    # quantro
    # check the number of sample.type is >= 2
    if(length(unique(doe$sample.type))<2){
      output$qt_notallow <- renderText(
        "The number of 'sample.type' is less than 2, Quantro can not be applied!"
      )
    }else{
      output$qt_notallow <- renderText(
        ""
      )
      ## raw
      qt_raw <- doQuantro(table_log2_clean_raw, doe)
      output$qt_qt_raw <- renderTable(
        qt_raw
      )
      qt_raw_a <- doQuantroAnova(table_log2_clean_raw, doe)
      output$qt_anov_raw <- renderTable(
        qt_raw_a, striped = TRUE
      )
      
      ## lfq
      if(LFQRun=="Y"){
        qt_lfq <- doQuantro(table_log2_clean_lfq, doe)
        output$qt_qt_lfq <- renderTable(
          qt_lfq
        )
        qt_lfq_a <- doQuantroAnova(table_log2_clean_lfq, doe)
        output$qt_anov_lfq <- renderTable(
          qt_lfq_a, striped = TRUE
        )
      }
      
      # use pg_ints_raw_melt_doe or pg_ints_lfq_melt_doe
      # boxplot
      ## raw
      pg_ints_raw_melt_doe$value = as.numeric(pg_ints_raw_melt_doe$value)
      output$qt_box_raw <- renderPlotly(
        plot_ly(pg_ints_raw_melt_doe %>% arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
          layout(yaxis=list(title="Log2 Intensity "),
                 xaxis=list(title="Sample", 
                            categoryorder= "array",
                            categoryarray=~sample.type))
      )
      
      ## lfq
      if(LFQRun=="Y"){
        pg_ints_lfq_melt_doe$value = as.numeric(pg_ints_lfq_melt_doe$value)
        output$qt_box_lfq <- renderPlotly(
          plot_ly(pg_ints_lfq_melt_doe %>% arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
            layout(yaxis=list(title="Log2 Intensity "),
                   xaxis=list(title="Sample", 
                              categoryorder= "array",
                              categoryarray=~sample.type))
        )
      }
      
      # density
      ## raw
      output$qt_den_raw <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe, aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()
        )
      )
      
      output$qt_den_raw_split <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe, aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()+
            facet_wrap(~ sample.type)
        )
      )
      ## lfq
      if(LFQRun=="Y"){
        output$qt_den_lfq <- renderPlotly(
          ggplotly(
            ggplot(pg_ints_raw_melt_doe, aes(x=value))+
              stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
              labs(x="Log2 Intensity", y="Density")+
              theme_minimal()
          )
        )
        output$qt_den_lfq_split <- renderPlotly(
          ggplotly(
            ggplot(pg_ints_raw_melt_doe, aes(x=value))+
              stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
              labs(x="Log2 Intensity", y="Density")+
              theme_minimal()+
              facet_wrap(~ sample.type)
          )
        )
      }
      
    }
    
    
    # quantro --- No control
    nc_doe <- doe %>% filter(type=="sample")
    if(length(unique(nc_doe$sample.type))<2){
      output$nc_qt_notallow <- renderText(
        "The number of 'sample.type' after filtering out the control samples is less than 2, Quantro can not be applied!"
      )
    }else{
      ## raw
      nc_qt_raw <- doQuantroNoControl(table_log2_clean_raw, doe)
      output$nc_qt_qt_raw <- renderTable(
        nc_qt_raw
      )
      nc_qt_raw_a <- doQuantroNoControlAnova(table_log2_clean_raw, doe)
      output$nc_qt_anov_raw <- renderTable(
        nc_qt_raw_a, striped = TRUE
      )
      
      ## lfq
      if(LFQRun=="Y"){
        nc_qt_lfq <- doQuantroNoControl(table_log2_clean_lfq, doe)
        output$nc_qt_qt_lfq <- renderTable(
          nc_qt_lfq
        )
        nc_qt_lfq_a <- doQuantroNoControlAnova(table_log2_clean_lfq, doe)
        output$nc_qt_anov_lfq <- renderTable(
          nc_qt_lfq_a, striped = TRUE
        )
      }
      
      # use pg_ints_raw_melt_doe or pg_ints_lfq_melt_doe
      # boxplot
      ## raw
      pg_ints_raw_melt_doe$value = as.numeric(pg_ints_raw_melt_doe$value)
      output$nc_qt_box_raw <- renderPlotly(
        plot_ly(pg_ints_raw_melt_doe %>% filter(type=="sample") %>%  arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
          layout(yaxis=list(title="Log2 Intensity "),
                 xaxis=list(title="Sample", 
                            categoryorder= "array",
                            categoryarray=~sample.type))
      )
      
      ## lfq
      if(LFQRun=="Y"){
        pg_ints_lfq_melt_doe$value = as.numeric(pg_ints_lfq_melt_doe$value)
        output$nc_qt_box_lfq <- renderPlotly(
          plot_ly(pg_ints_lfq_melt_doe %>% filter(type=="sample") %>% arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
            layout(yaxis=list(title="Log2 Intensity "),
                   xaxis=list(title="Sample", 
                              categoryorder= "array",
                              categoryarray=~sample.type))
        )
      }
      
      
      # density
      ## raw
      output$nc_qt_den_raw <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()
        )
      )
      
      output$nc_qt_den_raw_split <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()+
            facet_wrap(~ sample.type)
        )
      )
      ## lfq
      if(LFQRun=="Y"){
        output$nc_qt_den_lfq <- renderPlotly(
          ggplotly(
            ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
              stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
              labs(x="Log2 Intensity", y="Density")+
              theme_minimal()
          )
        )
        output$nc_qt_den_lfq_split <- renderPlotly(
          ggplotly(
            ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
              stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
              labs(x="Log2 Intensity", y="Density")+
              theme_minimal()+
              facet_wrap(~ sample.type)
          )
        )
      }
      
    }
    

    ################################################################################################
    # end
    
  })

  
  ###################################################################################
  
  ############
  # start QC #
  ############
  
  observeEvent(input$start,{
    qc_allow <- 0
    useSweetAlert()
    
    progressSweetAlert(
      session = session,
      id = "pb1",
      title = "In progress",
      striped = TRUE,
      value = 0
    )
    Sys.sleep(2)
    updateProgressBar(session = session, id = "pb1", value = 10)
    ### 1 proteinGroups.txt
    if(!isTruthy(input$pg_file)){
      closeSweetAlert()
      shinyalert(title = "The 'proteinGroups.txt' is required!!", 
                 text = "Please upload the file and try again.", 
                 type = "error")
      
    }else{
      
      pg_response()
      updateProgressBar(session = session, id = "pb1", value = 20)
      
      ### 6 prefix
      if(!isTruthy(input$prefix)){
        closeSweetAlert()
        shinyalert(title = "The 'prefix' is required!!", 
                   text = "Please type the 'prefix' and click 'Start QC' to try again.", 
                   type = "input",
                   inputType = "text",
                   inputId = "prefix")
        
        qc_allow = qc_allow + 1
      }else{
        shinyBS::closeAlert(session, "a6")
        updateProgressBar(session = session, id = "pb1", value = 30)
        
        ### 2 peptides.txt
        if(isTruthy(input$pt_file)){
          pt_response()
        }
        updateProgressBar(session = session, id = "pb1", value = 40)
        
        ### 3 msms.txt
        ### 4 summary.txt
        if(isTruthy(input$ms_file)){
          if(isTruthy(input$sum_file)){
            ms_sum_response()
          }
        }
        
        updateProgressBar(session = session, id = "pb1", value = 60)
        
        ### 5 doe.txt
        if(isTruthy(input$doe_file)){
          doe_response()
        }
        updateProgressBar(session = session, id = "pb1", value = 80)
        Sys.sleep(2)
        updateProgressBar(session = session, id = "pb1", value = 100)
        Sys.sleep(2)
        closeSweetAlert(session = session)
        shinyalert(title = "QC is finished!", text = "Please check the results by clicking the tabs in the sidebar.", type = "success")
        
      }
    }
    

    ### confirm the required inputs are there
    #if(qc_allow != 0){
    #  createAlert(session, "alert7", alertId = "a7",
    #              style = "danger",
    #              content = "Please provide required 'proteinGroups.txt' or prefix to perform the QC! ")
    #  updateProgressBar(session = session, id = "pb1", status = "danger", value = 100)
    #}else{
    #  shinyBS::closeAlert(session, "a7")
      
    #}
    
        
  })
  
  ###################################################################################
  
  observeEvent(input$st_selected_raw,{
    
    # read pg and doe first
    if(isTruthy(input$pg_file$datapath)){
      pg <- read.table(input$pg_file$datapath, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
    }else{
      pg <- readRDS("demo_run_file/pg_example.rds")
    }
    
    if(isTruthy(input$doe_file$datapath)){
      doe_file <- input$doe_file$datapath
      doe <- read.table(doe_file, sep = input$sep, header = T)
    }else{
      doe <- readRDS("demo_run_file/doe_example.rds")
    }
    
    parsed_pg <- process_pg(pg)
    
    # check LFQ
    pg_lfq_col <- parsed_pg %>% 
      select(starts_with("LFQ"))
    if(length(colnames(pg_lfq_col))>0){
      LFQRun <- "Y"
    }else{
      LFQRun <- "N"
    }
    
    colnames(doe) <- gsub(" ", ".", colnames(doe))
    colnames(doe) <- tolower(colnames(doe))
    doe$type <- tolower(doe$type)
    
    doe$sample.id <- gsub(" ", ".", doe$sample.id)
    doe$sample.id <- gsub("-", ".", doe$sample.id)
    
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    
    melt_table_name_raw <- getMeltwithNames(table_log2_clean_raw)
    melt_table_name_raw_doe <- melt_table_name_raw %>% 
      left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
    
    wide_sampletype_table_raw <- countMeanbyProteinSampletypeDoe(melt_table_name_raw_doe)
    name_col_st_raw <- wide_sampletype_table_raw[,c(1:3)]
    
    st <- input$st_selected_raw
    display_table_raw <- data.frame(name_col_st_raw, wide_sampletype_table_raw[,st])
    colnames(display_table_raw)[4] <- "Log2.Intensity"
    display_table_raw <- display_table_raw %>% arrange(desc(Log2.Intensity))
    display_table_raw <- display_table_raw[c(1:20),]
    display_table_raw <- rbind(c(st,"","",""), display_table_raw)
    
    output$top20_st_raw <- renderTable(
      #long_sampletype_top20_raw, striped = TRUE
      display_table_raw, striped=TRUE
    )
    shinyalert(title = "Table is updated!", text = "", type = "success")
  })
  
  
  ###################################################################################
  
  observeEvent(input$st_selected_lfq,{
    
    if(isTruthy(input$pg_file$datapath)){
      pg <- read.table(input$pg_file$datapath, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
    }else{
      pg <- readRDS("demo_run_file/pg_example.rds")
    }
    
    if(isTruthy(input$doe_file$datapath)){
      doe_file <- input$doe_file$datapath
      doe <- read.table(doe_file, sep = input$sep, header = T)
    }else{
      doe <- readRDS("demo_run_file/doe_example.rds")
    }
    
    parsed_pg <- process_pg(pg)
    
    # check LFQ
    pg_lfq_col <- parsed_pg %>% 
      select(starts_with("LFQ"))
    if(length(colnames(pg_lfq_col))>0){
      LFQRun <- "Y"
    }else{
      LFQRun <- "N"
    }
    
  
    colnames(doe) <- gsub(" ", ".", colnames(doe))
    colnames(doe) <- tolower(colnames(doe))
    doe$type <- tolower(doe$type)
    
    doe$sample.id <- gsub(" ", ".", doe$sample.id)
    doe$sample.id <- gsub("-", ".", doe$sample.id)
    
    table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
    
    melt_table_name_lfq <- getMeltwithNames(table_log2_clean_lfq)
    melt_table_name_lfq_doe <- melt_table_name_lfq %>% 
      left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
    
    wide_sampletype_table_lfq <- countMeanbyProteinSampletypeDoe(melt_table_name_lfq_doe)
    name_col_st_lfq <- wide_sampletype_table_lfq[,c(1:3)]
    
    st <- input$st_selected_lfq
    display_table_lfq <- data.frame(name_col_st_lfq, wide_sampletype_table_lfq[,st])
    colnames(display_table_lfq)[4] <- "Log2.LFQ.Intensity"
    display_table_lfq <- display_table_lfq %>% arrange(desc(Log2.LFQ.Intensity))
    display_table_lfq <- display_table_lfq[c(1:20),]
    display_table_lfq <- rbind(c(st,"","",""), display_table_lfq)
    
    output$top20_st_lfq <- renderTable(
      #long_sampletype_top20_lfq, striped = TRUE
      display_table_lfq, striped=TRUE
    )
    shinyalert(title = "Table is updated!", text = "", type = "success")
  })
  
  
  ###################################################################################
  observeEvent(input$t_selected_raw,{
    
    if(isTruthy(input$pg_file$datapath)){
      pg <- read.table(input$pg_file$datapath, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
    }else{
      pg <- readRDS("demo_run_file/pg_example.rds")
    }
    
    if(isTruthy(input$doe_file$datapath)){
      doe_file <- input$doe_file$datapath
      doe <- read.table(doe_file, sep = input$sep, header = T)
    }else{
      doe <- readRDS("demo_run_file/doe_example.rds")
    }
    
    parsed_pg <- process_pg(pg)
    
    # check LFQ
    pg_lfq_col <- parsed_pg %>% 
      select(starts_with("LFQ"))
    if(length(colnames(pg_lfq_col))>0){
      LFQRun <- "Y"
    }else{
      LFQRun <- "N"
    }
    
    
    colnames(doe) <- gsub(" ", ".", colnames(doe))
    colnames(doe) <- tolower(colnames(doe))
    doe$type <- tolower(doe$type)
    
    doe$sample.id <- gsub(" ", ".", doe$sample.id)
    doe$sample.id <- gsub("-", ".", doe$sample.id)
    
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    
    melt_table_name_raw <- getMeltwithNames(table_log2_clean_raw)
    melt_table_name_raw_doe <- melt_table_name_raw %>% 
      left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
    
    wide_type_table_raw <- countMeanbyProteinTypeDoe(melt_table_name_raw_doe)
    name_col_t_raw <- wide_type_table_raw[,c(1:3)]
    
    t <- input$t_selected_raw
    display_table_raw <- data.frame(name_col_t_raw, wide_type_table_raw[,t])
    colnames(display_table_raw)[4] <- "Log2.Intensity"
    display_table_raw <- display_table_raw %>% arrange(desc(Log2.Intensity))
    display_table_raw <- display_table_raw[c(1:20),]
    display_table_raw <- rbind(c(t,"","",""), display_table_raw)
    
    output$top20_t_raw <- renderTable(
      #long_sampletype_top20_raw, striped = TRUE
      display_table_raw, striped=TRUE
    )
    shinyalert(title = "Table is updated!", text = "", type = "success")
  })
  
  ###################################################################################
  observeEvent(input$t_selected_lfq,{
    # read pg and doe first

    if(isTruthy(input$pg_file$datapath)){
      pg <- read.table(input$pg_file$datapath, header = T, sep = "\t", stringsAsFactors = F, quote = "\"", comment.char = "")
    }else{
      pg <- readRDS("demo_run_file/pg_example.rds")
    }
    
    if(isTruthy(input$doe_file$datapath)){
      doe_file <- input$doe_file$datapath
      doe <- read.table(doe_file, sep = input$sep, header = T)
    }else{
      doe <- readRDS("demo_run_file/doe_example.rds")
    }
    
    parsed_pg <- process_pg(pg)
    
    # check LFQ
    pg_lfq_col <- parsed_pg %>% 
      select(starts_with("LFQ"))
    if(length(colnames(pg_lfq_col))>0){
      LFQRun <- "Y"
    }else{
      LFQRun <- "N"
    }
    
    colnames(doe) <- gsub(" ", ".", colnames(doe))
    colnames(doe) <- tolower(colnames(doe))
    doe$type <- tolower(doe$type)
    
    doe$sample.id <- gsub(" ", ".", doe$sample.id)
    doe$sample.id <- gsub("-", ".", doe$sample.id)
    
    table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
    
    melt_table_name_lfq <- getMeltwithNames(table_log2_clean_lfq)
    melt_table_name_lfq_doe <- melt_table_name_lfq %>% 
      left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
    
    wide_type_table_lfq <- countMeanbyProteinTypeDoe(melt_table_name_lfq_doe)
    name_col_t_lfq <- wide_type_table_lfq[,c(1:3)]
    
    t <- input$t_selected_lfq
    display_table_lfq <- data.frame(name_col_t_lfq, wide_type_table_lfq[,t])
    colnames(display_table_lfq)[4] <- "Log2.Intensity"
    display_table_lfq <- display_table_lfq %>% arrange(desc(Log2.Intensity))
    display_table_lfq <- display_table_lfq[c(1:20),]
    display_table_lfq <- rbind(c(t,"","",""), display_table_lfq)
    
    output$top20_t_lfq <- renderTable(
      #long_sampletype_top20_raw, striped = TRUE
      display_table_lfq, striped=TRUE
    )
    shinyalert(title = "Table is updated!", text = "", type = "success")
  })
  
  
  
  
  
  ###################################################################################
  
  ##############################################
  # demo run - responses upon the example file #
  ##############################################
  
  #### 1 proteinGroups.txt
  pg_ex_response <- reactive({
    pg_ex <- readRDS("demo_run_file/pg_example.rds")
    parsed_pg <- process_pg(pg_ex)
    
    ################################################################################################
    # datasets for download

    table_log2_raw <- getLog2Table(parsed_pg, "raw")
    table_log2_lfq <- getLog2Table(parsed_pg, "lfq")
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
    table_contam_log2_raw <- getLog2ContamTable(parsed_pg, "raw")
    table_contam_log2_lfq <- getLog2ContamTable(parsed_pg, "lfq")
    
    output$dtable_log2_raw <- downloadHandler(
      filename = function(){
        paste0("demo","_log2_raw_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_log2_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dtable_log2_lfq <- downloadHandler(
      filename = function(){
        paste0("demo","_log2_lfq_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_log2_lfq, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dtable_log2_clean_raw <- downloadHandler(
      filename = function(){
        paste0("demo","_log2_raw_clean_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_log2_clean_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dtable_log2_clean_lfq <- downloadHandler(
      filename = function(){
        paste0("demo","_log2_lfq_clean_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_log2_clean_lfq, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dtable_contam_log2_raw <- downloadHandler(
      filename = function(){
        paste0("demo","_log2_raw_contaminant_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_contam_log2_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dtable_contam_log2_lfq <- downloadHandler(
      filename = function(){
        paste0("demo","_log2_lfq_contaminant_data_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(table_contam_log2_lfq, file, sep = "\t", row.names = FALSE)
      }
    )
    
    ##########################################################################################################
    # expCheck - spikein
    
    spike_dat_raw <- getSpikeInfo(parsed_pg, "raw")
    spike_dat_lfq <- getSpikeInfo(parsed_pg, "lfq")
    
    spike_raw_log2_mean <- mean(spike_dat_raw$log2.intensity)
    spike_raw_log2_sd <- sd(spike_dat_raw$log2.intensity)
    
    spike_lfq_log2_mean <- mean(spike_dat_lfq$log2.intensity)
    spike_lfq_log2_sd <- sd(spike_dat_lfq$log2.intensity)
    
    cutline_below_raw = spike_raw_log2_mean - 2*spike_raw_log2_sd
    cutline_above_raw = spike_raw_log2_mean + 2*spike_raw_log2_sd
    
    cutline_below_lfq = spike_lfq_log2_mean - 2*spike_lfq_log2_sd
    cutline_above_lfq = spike_lfq_log2_mean + 2*spike_lfq_log2_sd
    
    below_raw <- getSpikeBelow2sdSample(spike_dat_raw)
    below_lfq <- getSpikeBelow2sdSample(spike_dat_lfq)
    
    above_raw <- getSpikeAbove2sdSample(spike_dat_raw)
    above_lfq <- getSpikeAbove2sdSample(spike_dat_lfq)
    
    spike_dat_raw <- spike_dat_raw %>% 
      mutate(outlier = ifelse(sample %in% below_raw$sample | sample %in% above_raw$sample, "#FF6666", "#00CCCC"))
    
    spike_dat_lfq <- spike_dat_lfq %>% 
      mutate(outlier = ifelse(sample %in% below_lfq$sample | sample %in% above_lfq$sample, "#FF6666", "#00CCCC"))
    
    output$spike_plot_raw <- renderPlotly({
      n = nrow(spike_dat_raw)
      plot_ly(spike_dat_raw, x=~sample, y=~log2.intensity, type = "scatter", color =~ I(outlier), name="Spikein") %>% 
        add_segments(y=cutline_above_raw, yend=cutline_above_raw, x=spike_dat_raw$sample[1], xend=spike_dat_raw$sample[n], color=I("gray"), name="Mean+2SD") %>% 
        add_segments(y=spike_raw_log2_mean, yend=spike_raw_log2_mean, x=spike_dat_raw$sample[1], xend=spike_dat_raw$sample[n], color=I("black"), name="Mean") %>% 
        add_segments(y=cutline_below_raw, yend=cutline_below_raw, x=spike_dat_raw$sample[1], xend=spike_dat_raw$sample[n], color=I("gray"), name="Mean-2SD") 
    })
    
    output$spike_plot_lfq <- renderPlotly({
      n = nrow(spike_dat_lfq)
      plot_ly(spike_dat_lfq, x=~sample, y=~log2.intensity, type = "scatter", color =~ I(outlier), name="Spikein") %>% 
        add_segments(y=cutline_above_lfq, yend=cutline_above_lfq, x=spike_dat_lfq$sample[1], xend=spike_dat_lfq$sample[n], color=I("gray"), name="Mean+2SD") %>% 
        add_segments(y=spike_lfq_log2_mean, yend=spike_lfq_log2_mean, x=spike_dat_lfq$sample[1], xend=spike_dat_lfq$sample[n], color=I("black"), name="Mean") %>% 
        add_segments(y=cutline_below_lfq, yend=cutline_below_lfq, x=spike_dat_lfq$sample[1], xend=spike_dat_lfq$sample[n], color=I("gray"), name="Mean-2SD") 
    })
    
    output$spikelow_raw <- renderTable({
      below_raw %>% select(sample, log2.intensity)
    })
    
    output$spikehigh_raw <- renderTable({
      above_raw %>% select(sample, log2.intensity)
    })
    
    output$spikelow_lfq <- renderTable({
      below_lfq %>% select(sample, log2.intensity)
    })
    
    output$spikehigh_lfq <- renderTable({
      above_lfq %>% select(sample, log2.intensity)
    })
    
    
    ##################################################################################################################################################################################
    # contaminants
    
    ## raw
    raw_con_info <- getContamInfo(parsed_pg, "raw")
    raw_con_t1 <- getContamTable1(raw_con_info)
    
    raw_con_t1g <- gather(raw_con_t1, key = "Contaminant.type", value = "Total.intensity", -sample)
    
    output$contam_raw <- renderPlotly(
      ggplotly(
        ggplot(raw_con_t1g, aes(x=sample, y=Total.intensity, color=Contaminant.type))+
          geom_point()+
          labs(y="Total intensity (1og10 scale)")+
          theme_minimal()+
          theme(axis.text.x = element_text(angle=90))+
          scale_y_log10()
      )
    )
    
    ## lfq
    lfq_con_info <- getContamInfo(parsed_pg, "lfq")
    lfq_con_t1 <- getContamTable1(lfq_con_info)
    
    lfq_con_t1g <- gather(lfq_con_t1, key = "Contaminant.type", value = "Total.intensity", -sample)
    
    output$contam_lfq <- renderPlotly(
      ggplotly(
        ggplot(lfq_con_t1g, aes(x=sample, y=Total.intensity, color=Contaminant.type))+
          geom_point()+
          labs(y="Total intensity (1og10 scale)")+
          theme_minimal()+
          theme(axis.text.x = element_text(angle=90))+
          scale_y_log10()
      )
    )
    
    ###############################################################################################################
    # protein counts
    ## raw
    raw_prog_count <- getProteinCount(parsed_pg, "raw")
    output$pg_count_raw <- renderPlotly(
      plot_ly(raw_prog_count, x=~sample, y=~count, type = "scatter") %>% 
        layout(yaxis=list(title="Number of Protein Group"))
    )
    
    ## lfq
    lfq_prog_count <- getProteinCount(parsed_pg, "lfq")
    output$pg_count_lfq <- renderPlotly(
      plot_ly(lfq_prog_count, x=~sample, y=~count, type = "scatter", color=I("#FF6666")) %>% 
        layout(yaxis=list(title="Number of Protein Group"))
    )
    
    ###############################################################################################################
    # protein intensity
    
    # total by individual
    ## raw
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    prog_int_sum_raw <- getTotalIntensity(table_log2_clean_raw)
    
    output$pg_int_raw <- renderPlotly(
      plot_ly(prog_int_sum_raw, x=~sample, y=~Total.Log2.Intensity, type="scatter")
    )
    
    ## lfq
    table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
    prog_int_sum_lfq <- getTotalIntensity(table_log2_clean_lfq)
    
    output$pg_int_lfq <- renderPlotly(
      plot_ly(prog_int_sum_lfq, x=~sample, y=~Total.Log2.Intensity, type="scatter", color=I("#FF6666"))
    )
    
    ###############################################################################################################
    # protein intensity
    
    # violin distribution
    ## raw
    pg_ints_raw <- getIntensity(table_log2_clean_raw)
    pg_ints_raw_melt <- melt(pg_ints_raw, id.vars="sample")
    pg_ints_raw_melt$value = as.numeric(pg_ints_raw_melt$value)
    
    output$pg_int_vio_raw <- renderPlotly(
      plot_ly(pg_ints_raw_melt, x=~sample, y=~value, type="violin", split =~sample, color = I("#3366CC"), showlegend =FALSE) %>% 
        layout(yaxis=list(title="Log2 Intensity"), xaxis=list(title="Sample"))
      
    )
    ## lfq
    pg_ints_lfq <- getIntensity(table_log2_clean_lfq)
    pg_ints_lfq_melt <- melt(pg_ints_lfq, id.vars="sample")
    pg_ints_lfq_melt$value = as.numeric(pg_ints_lfq_melt$value)
    
    output$pg_int_vio_lfq <- renderPlotly(
      plot_ly(pg_ints_lfq_melt, x=~sample, y=~value, type="violin", split =~sample, color = I("#FF6666"), showlegend =FALSE) %>% 
        layout(yaxis=list(title="Log2 Intensity"), xaxis=list(title="Sample"))
      
    )
    
    ###############################################################################################################
    # protein intensity
    
    # Raw vs LFQ overll
    # get the melt data from raw and lfq: pg_ints_raw/lfq_melt
    melt_table <- data.frame(raw=pg_ints_raw_melt$value,
                             lfq=pg_ints_lfq_melt$value)
    
    output$pg_int_rawlfq <- renderPlotly(
      plot_ly(melt_table, x=~raw, y=~lfq, type="scatter", marker=list(opacity=0.3)) %>% 
        layout(xaxis=list(title="Raw Intensity"), yaxis=list(title="LFQ Intensity"))
    )
    
    ###############################################################################################################
    # CV
    ## raw
    cv_raw <- ddply(pg_ints_raw_melt, .(variable), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_raw <- renderPlotly(
      ggplotly(
        ggplot(cv_raw, aes(x=cv))+
          geom_density(alpha=0.5, fill="#3366CC")+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    
    ## lfq
    cv_lfq <- ddply(pg_ints_lfq_melt, .(variable), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_lfq <- renderPlotly(
      ggplotly(
        ggplot(cv_lfq, aes(x=cv))+
          geom_density(alpha=0.5, fill="#FF6666")+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    
    ###############################################################################################################
    # top20
    # get the mean intensity of all the proteins from table_log2_clean_raw/lfq
    
    ## raw
    melt_table_name_raw <- getMeltwithNames(table_log2_clean_raw)
    annotated_mean_table_raw <- countMeanbyProtein(melt_table_name_raw)
    top20.raw <- annotated_mean_table_raw[1:20,]
    output$top20_raw <- renderTable(
      top20.raw, striped = TRUE
    )
    
    ## lfq
    melt_table_name_lfq <- getMeltwithNames(table_log2_clean_lfq)
    annotated_mean_table_lfq <- countMeanbyProtein(melt_table_name_lfq)
    top20.lfq <- annotated_mean_table_lfq[1:20,]
    output$top20_lfq <- renderTable(
      top20.lfq, striped = TRUE
    )
    
    output$dannotated_mean_table_raw <- downloadHandler(
      filename = function(){
        paste0("demo","_mean_log2_raw_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(annotated_mean_table_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dannotated_mean_table_lfq <- downloadHandler(
      filename = function(){
        paste0("demo","_mean_log2_lfq_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(annotated_mean_table_lfq, file, sep = "\t", row.names = FALSE)
      }
    )
    
    ###############################################################################################################
    # pca
    ## raw
    pca_table_raw <- getPCAtable(table_log2_clean_raw)
    pc12_raw <- getPC12Percentage(table_log2_clean_raw)
    
    output$pca_raw <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_raw, aes(x=PC1, y=PC2, color=sample) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_raw[1], "%)"),
      #         y=paste0("PC2 (", pc12_raw[2], "%)"))
      #)
      plot_ly(pca_table_raw, x=~PC1, y=~PC2, type = "scatter", color=~sample, text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE), 
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
    )
    ## lfq
    pca_table_lfq <- getPCAtable(table_log2_clean_lfq)
    pc12_lfq <- getPC12Percentage(table_log2_clean_lfq)
    
    
    output$pca_lfq <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_lfq, aes(x=PC1, y=PC2, color=sample) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_lfq[1], "%)"),
      #         y=paste0("PC2 (", pc12_lfq[2], "%)"))
      #)
      plot_ly(pca_table_lfq, x=~PC1, y=~PC2, type = "scatter", color=~sample, text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE), 
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
      
    )
    
    ###############################################################################################################
    # end

  })
  
  #### 2 peptides.txt
  pt_ex_response <- reactive({
    pt <- readRDS("demo_run_file/pt_example.rds")
    
    # get the average of Cysteine
    cys_avg = sum(pt$C.Count)/nrow(pt)
    # get the average of missed cleavage
    ms_clvg_avg = sum(pt$Missed.cleavage)/nrow(pt)
    # make checkpoint table
    avg_table = data.frame("Feature" = c("Cysteine","Missed Cleavage"),
                           "Percentage" = c(paste0(format(round(cys_avg*100, 2), nsmall = 2), " %"), paste0(format(round(ms_clvg_avg*100, 2), nsmall = 2), " %")))
    output$expck_table <- renderTable(
      avg_table, striped = T, bordered = T
    )
    
    # get peptide count by individual sample
    
    pep_count <- getPepCount(pt)
    
    output$pt_count <- renderPlotly(
      plot_ly(pep_count, x=~sample, y=~count, type = "scatter") %>% 
        layout(yaxis=list(title="Number of Peptide"))
    )
    
    # get peptide total log2 intensity
    
    pep_int_sum <- getTotalIntensityPep(pt)
    
    output$pt_int <- renderPlotly(
      plot_ly(pep_int_sum, x=~sample, y=~Total.Log2.Intensity, type="scatter", color = I("#3366CC"))
    )
    
  })
  
  #### 3 msms.txt
  #### 4 summary.txt
  ms_sum_ex_response <- reactive({
    ms <- readRDS("demo_run_file/ms_example.rds")
    sumtxt <- readRDS("demo_run_file/sm_example.rds")
    
    ms_clvg_table <- ms %>% 
      select(Raw.file, Missed.cleavages) %>% 
      left_join(sumtxt %>% select(Raw.file, Experiment)) %>% 
      dplyr::count(Experiment, Missed.cleavages) %>% 
      spread(Missed.cleavages, n) %>% 
      mutate(sum = `0`+`1`+`2`) %>% 
      mutate(`percentage.0` = `0`/sum*100,
             `percentage.1` = `1`/sum*100,
             `percentage.2` = `2`/sum*100)
    
    output$msclvg_persample <- renderPlotly(
      plot_ly(ms_clvg_table, x= ~Experiment, y=~`percentage.0`, type = 'bar', name = "percentage.0") %>% 
        add_trace(y = ~`percentage.1`, name = "percentage.1") %>%
        add_trace(y = ~`percentage.2`, name = "percentage.2") %>%
        layout(yaxis=list(title = 'Percentage (%)'), barmode = 'group')
    )
    
  })
  
  
  #### 5 doe.txt
  doe_ex_response <- reactive({
    doe <- readRDS("demo_run_file/doe_example.rds")
    doe$sample.id <- gsub(" ", ".", doe$sample.id)
    doe$sample.id <- gsub("-", ".", doe$sample.id)
    
    colnames(doe) <- tolower(colnames(doe))
    doe$type <- tolower(doe$type)

    ################################################################################################
    # read pg as well
    pg <- readRDS("demo_run_file/pg_example.rds")
    parsed_pg <- process_pg(pg)
    
    ################################################################################################
    # data structure
    int_col <- pg %>% 
      select(starts_with("Intensity."))
    pg_sample_count <- ncol(int_col)
    
    output$doe_descrp <- renderText({
      paste0("There are ", pg_sample_count, " samples in the data")
    })
    
    sampletype_table <- dplyr::count(doe, sample.type)
    
    output$datastructure_table <- renderTable(
      sampletype_table, striped = T, bordered = T
    )
    
    output$datastructure_plot <- renderPlotly(
      plot_ly(sampletype_table, x= ~sample.type, y=~n, type = 'bar', color = ~sample.type) %>% 
        layout(yaxis=list(title = 'Count'))
    )
    
    ################################################################################################
    # contaminant type with sample type #pg
    ## raw
    raw_con_info <- getContamInfo(parsed_pg, "raw")
    raw_con_t2 <- getContamTable2(raw_con_info, doe)
    
    output$contam_raw_st <- renderPlotly(
      plot_ly(raw_con_t2, y=~log2(intensity), x=~sample.type, color=~kind_contam, type = "box") %>% 
        layout(yaxis=list(title="Log2 Raw Intensity" ),boxmode = "group")
    )
    
    ## lfq
    lfq_con_info <- getContamInfo(parsed_pg, "raw")
    lfq_con_t2 <- getContamTable2(lfq_con_info, doe)
    
    output$contam_lfq_st <- renderPlotly(
      plot_ly(lfq_con_t2, y=~log2(intensity), x=~sample.type, color=~kind_contam, type = "box") %>% 
        layout(yaxis=list(title="Log2 LFQ Intensity" ),boxmode = "group")
    )
    
    ################################################################################################
    # protein counts 
    # sample type
    ## raw
    raw_prog_count_doe <- getProteinCountDoe(parsed_pg, "raw", doe)
    output$pg_count_st_raw <- renderPlotly(
      plot_ly(raw_prog_count_doe, x=~sample.type, y=~count, type = "box") %>% 
        layout(yaxis=list(title="Number of Protein Groups"))
    )
    
    ## lfq
    lfq_prog_count_doe <- getProteinCountDoe(parsed_pg, "lfq", doe)
    output$pg_count_st_lfq <- renderPlotly(
      plot_ly(lfq_prog_count_doe, x=~sample.type, y=~count, type = "box", color=I("#FF6666")) %>% 
        layout(yaxis=list(title="Number of Protein Groups"))
    )
    
    # type cdf plot
    ## raw
    per_data_melt_raw <- getExpressedSamplePercentPerProteinEachType(parsed_pg, "raw", doe)
    output$pg_count_t_raw <- renderPlot(
      ggplot(per_data_melt_raw, aes(x=Value, color=Type))+
        stat_ecdf(geom = "step")+
        scale_x_reverse()+
        labs(x="Percent of samples expressed per protein (%)", y="Proportion of expressed proteins")+
        theme_minimal()
    )
    
    ## lfq
    per_data_melt_lfq <- getExpressedSamplePercentPerProteinEachType(parsed_pg, "lfq", doe)
    output$pg_count_t_lfq <- renderPlot(
      ggplot(per_data_melt_lfq, aes(x=Value, color=Type))+
        stat_ecdf(geom = "step")+
        scale_x_reverse()+
        labs(x="Percent of samples expressed per protein (%)", y="Proportion of expressed proteins")+
        theme_minimal()
    )
    
    ################################################################################################
    # peptide counts
    # get peptide count by sample.type and type
    
    pt <- readRDS("demo_run_file/pt_example.rds")
    pep_count_doe <- getPepCountDoe(pt, doe)
    
    output$pt_count_st <- renderPlotly(
      plot_ly(pep_count_doe, x=~sample.type, y=~count, type = "box") %>% 
        layout(yaxis=list(title="Number of Peptide"), xaxis=list(title="Sample Type"))
    )
    
    output$pt_count_t <- renderPlotly(
      plot_ly(pep_count_doe, x=~type, y=~count, type = "box", color=~type) %>% 
        layout(yaxis=list(title="Number of Peptide"), xaxis=list(title="Type"))
    )
    
    ################################################################################################
    # protein intensity - mean proportion plot
    ## raw
    table_log2_clean_raw <- getLog2CleanTable(parsed_pg, "raw")
    pg_ints_raw <- getIntensity(table_log2_clean_raw)
    pg_ints_raw_melt <- melt(pg_ints_raw, id.vars="sample")
    pg_ints_raw_melt$value = as.numeric(pg_ints_raw_melt$value)
    
    mp_raw <- getMedianIntPlotData(pg_ints_raw_melt, doe)
    
    output$pg_int_mpp_raw <- renderPlotly(
      plot_ly(mp_raw, type = "scatter", mode = "markers") %>% 
        add_trace(y = ~All, name = 'All sample', marker = list(opacity=0.2, color="#3366CC")) %>% 
        layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
               xaxis=list(title="Protein Group"))
    )
    output$pg_int_mpp_t_raw <- renderPlotly(
      plot_ly(mp_raw, type = "scatter", mode = "markers") %>% 
        add_trace(y = ~Control, name = 'Control', marker = list(opacity=0.2)) %>%
        add_trace(y = ~Sample, name = 'Sample', marker = list(opacity=0.2)) %>% 
        layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
               xaxis=list(title="Protein Group"))
    )
    ## lfq
    table_log2_clean_lfq <- getLog2CleanTable(parsed_pg, "lfq")
    pg_ints_lfq <- getIntensity(table_log2_clean_lfq)
    pg_ints_lfq_melt <- melt(pg_ints_lfq, id.vars="sample")
    pg_ints_lfq_melt$value = as.numeric(pg_ints_lfq_melt$value)
    
    mp_lfq <- getMedianIntPlotData(pg_ints_lfq_melt, doe)
    
    output$pg_int_mpp_lfq <- renderPlotly(
      plot_ly(mp_lfq, type = "scatter", mode = "markers") %>% 
        add_trace(y = ~All, name = 'All sample', marker = list(opacity=0.2, color="#FF6666")) %>% 
        layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
               xaxis=list(title="Protein Group"))
    )
    output$pg_int_mpp_t_lfq <- renderPlotly(
      plot_ly(mp_lfq, type = "scatter", mode = "markers") %>% 
        add_trace(y = ~Control, name = 'Control', marker = list(opacity=0.2)) %>%
        add_trace(y = ~Sample, name = 'Sample', marker = list(opacity=0.2)) %>% 
        layout(yaxis=list(title="Mean Intensity/\nSum of Mean Intensities (log scale)", type='log'), 
               xaxis=list(title="Protein Group"))
    )
    
    ################################################################################################
    # protein intensity - cdf plot
    ## raw
    pg_ints_raw_melt_doe = pg_ints_raw_melt %>% 
      left_join(doe, by=c("sample"="sample.id"))
    
    output$pg_int_cdf_st_raw <- renderPlot(
      ggplot(pg_ints_raw_melt_doe, aes(x=value+0.000001, color=sample))+
        stat_ecdf()+
        labs(x = "Log2 Intensities by Protein", 
             y = "Fraction of Library") + 
        theme_minimal() +
        facet_wrap(~sample.type, ncol=2)
    )
    
    output$pg_int_vio_raw_doe <- renderPlotly(
      ggplotly(
        ggplot(pg_ints_raw_melt_doe, aes(x=reorder(sample,run.order), y=value, color=sample.type, fill=sample.type))+
          geom_violin()+
          theme_minimal()+
          labs(x="Sample by Run order", y="Log2 Intensity")+
          theme(axis.text.x = element_text(angle = 90))
      )
    )
    
    ## lfq
    pg_ints_lfq_melt_doe = pg_ints_lfq_melt %>% 
      left_join(doe, by=c("sample"="sample.id"))
    
    output$pg_int_cdf_st_lfq <- renderPlot(
      ggplot(pg_ints_lfq_melt_doe, aes(x=value+0.000001, color=sample))+
        stat_ecdf()+
        labs(x = "Log2 Intensities by Protein", 
             y = "Fraction of Library") + 
        theme_minimal() +
        facet_wrap(~sample.type, ncol=2)
    )
    
    output$pg_int_vio_lfq_doe <- renderPlotly(
      ggplotly(
        ggplot(pg_ints_lfq_melt_doe, aes(x=reorder(sample,run.order), y=value, color=sample.type, fill=sample.type))+
          geom_violin()+
          theme_minimal()+
          labs(x="Sample by Run order", y="Log2 Intensity")+
          theme(axis.text.x = element_text(angle = 90))
      )
    )
    
    ################################################################################################
    # protein intensity - Raw vs LFQ sample.type
    # get the melt data from raw and lfq: pg_ints_raw/lfq_melt_doe
    melt_table_doe <- data.frame(sample=pg_ints_raw_melt_doe$sample,
                                 type=pg_ints_raw_melt_doe$type,
                                 sample.type=pg_ints_raw_melt_doe$sample.type,
                                 raw=pg_ints_raw_melt_doe$value,
                                 lfq=pg_ints_lfq_melt_doe$value)
    
    st_num <- length(unique(melt_table_doe$sample.type))
    if(st_num<=3){
      n=2
    }else{
      n=3
    }
    n <- as.numeric(n)
    
    output$pg_int_st_rawlfq <- renderPlotly(
      melt_table_doe %>% 
        group_by(sample.type) %>% 
        group_map(~ plot_ly(data = ., x=~raw, y=~lfq, color=~sample.type, type="scatter", marker=list(opacity=0.3)), keep = TRUE) %>% 
        subplot(nrows = n, shareX = TRUE, shareY = TRUE) %>% 
        layout(xaxis=list(title="Raw Intensity"), yaxis=list(title="LFQ Intensity"))
    )
    
    ################################################################################################
    # CV with doe, sample type
    ## raw
    cv_raw_st <- ddply(pg_ints_raw_melt_doe, .(variable, sample.type), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_st_raw <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_st, aes(x=cv, fill=sample.type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    output$pg_int_cv_st_raw_split <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_st, aes(x=cv, fill=sample.type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")+
          facet_grid(~sample.type)
      )
    )
    
    ## lfq
    cv_lfq_st <- ddply(pg_ints_lfq_melt_doe, .(variable, sample.type), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_st_lfq <- renderPlotly(
      ggplotly(
        ggplot(cv_lfq_st, aes(x=cv, fill=sample.type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    output$pg_int_cv_st_lfq_split <- renderPlotly(
      ggplotly(
        ggplot(cv_lfq_st, aes(x=cv, fill=sample.type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")+
          facet_grid(~sample.type)
      )
    )
    
    # CV with doe, type
    ## raw
    cv_raw_t <- ddply(pg_ints_raw_melt_doe, .(variable, type), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_t_raw <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_t, aes(x=cv, fill=type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    
    output$pg_int_cv_t_raw_split <- renderPlotly(
      ggplotly(
        ggplot(cv_raw_t, aes(x=cv, fill=type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")+
          facet_grid(~type)
      )
    )
    
    ## lfq
    cv_lfq_t <- ddply(pg_ints_lfq_melt_doe, .(variable, type), summarise, cv=sd(value)/mean(value)*100)
    
    output$pg_int_cv_t_lfq <- renderPlotly(
      ggplotly(
        ggplot(cv_lfq_t, aes(x=cv, fill=type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")
      )
    )
    
    output$pg_int_cv_t_lfq_split <- renderPlotly(
      ggplotly(
        ggplot(cv_lfq_t, aes(x=cv, fill=type))+
          geom_density(alpha=0.5)+
          theme_minimal()+
          labs(x="Coefficient of Variation (%)", y="Density")+
          facet_grid(~type)
      )
    )
    
    ################################################################################################
    # top20 sample.type
    ## raw
    melt_table_name_raw <- getMeltwithNames(table_log2_clean_raw)
    melt_table_name_raw_doe <- melt_table_name_raw %>% 
      left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
    
    
    wide_sampletype_table_raw <- countMeanbyProteinSampletypeDoe(melt_table_name_raw_doe)
    name_col_st_raw <- wide_sampletype_table_raw[,c(1:3)]
    measure_col_st_raw <- wide_sampletype_table_raw[,c(4:ncol(wide_sampletype_table_raw))]
    
    output$sampletype_raw <- renderUI({
      st_list <- sort(unique(doe$sample.type))
      selectInput("st_selected_raw", h4("Choose Sample Type"),
                  choices = st_list, selected = 1)
    })
  
    ## lfq
    melt_table_name_lfq <- getMeltwithNames(table_log2_clean_lfq)
    melt_table_name_lfq_doe <- melt_table_name_lfq %>% 
      left_join(doe %>% select(sample.id, sample.type, type), by=c("sample"="sample.id"))
    
    
    wide_sampletype_table_lfq <- countMeanbyProteinSampletypeDoe(melt_table_name_lfq_doe)
    name_col_st_lfq <- wide_sampletype_table_lfq[,c(1:3)]
    measure_col_st_lfq <- wide_sampletype_table_lfq[,c(4:ncol(wide_sampletype_table_lfq))]
    
    output$sampletype_lfq <- renderUI({
      st_list <- sort(unique(doe$sample.type))
      selectInput("st_selected_lfq", h4("Choose Sample Type"),
                  choices = st_list, selected = 1)
    })
    
    
    output$dwide_sampletype_sort_table_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_meanBySampletype_log2_raw_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(wide_sampletype_table_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dwide_sampletype_sort_table_lfq <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_meanBySampletype_log2_lfq_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(wide_sampletype_table_lfq, file, sep = "\t", row.names = FALSE)
      }
    )
    
    # top20 type
    ## raw
    
    wide_type_table_raw <- countMeanbyProteinTypeDoe(melt_table_name_raw_doe)
    name_col_t_raw <- wide_type_table_raw[,c(1:3)]
    measure_col_t_raw <- wide_type_table_raw[,c(4:ncol(wide_type_table_raw))]
    
    output$type_raw <- renderUI({
      t_list <- sort(unique(doe$type))
      selectInput("t_selected_raw", h4("Choose Type"),
                  choices = t_list, selected = 1)
    })
    
    ## lfq
    
    wide_type_table_lfq <- countMeanbyProteinTypeDoe(melt_table_name_lfq_doe)
    name_col_t_lfq <- wide_type_table_lfq[,c(1:3)]
    measure_col_t_lfq <- wide_type_table_lfq[,c(4:ncol(wide_type_table_lfq))]
    
    output$type_lfq <- renderUI({
      t_list <- sort(unique(doe$type))
      selectInput("t_selected_lfq", h4("Choose Type"),
                  choices = t_list, selected = 1)
    })  
    
    
    output$dwide_type_sort_table_raw <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_meanByType_log2_raw_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(wide_type_table_raw, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$dwide_type_sort_table_lfq <- downloadHandler(
      filename = function(){
        paste0(input$prefix,"_meanByType_log2_lfq_",Sys.Date(),".tsv")
      },
      content = function(file){
        write.table(wide_type_table_lfq, file, sep = "\t", row.names = FALSE)
      }
    )
    
    ################################################################################################
    # pca sample    sample.type, run.order
    ## raw
    pca_table_raw_doe <- getPCAtableDOE(table_log2_clean_raw, doe)
    pc12_raw <- getPC12Percentage(table_log2_clean_raw)
    
    output$pca_st_raw <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_raw_doe, aes(x=PC1, y=PC2, color=sample.type, shape=type) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_raw[1], "%)"),
      #         y=paste0("PC2 (", pc12_raw[2], "%)"))
      #)
      plot_ly(pca_table_raw_doe, x=~PC1, y=~PC2, type = "scatter", color=~sample.type, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE),
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
      
      
    )
    
    output$pca_ro_raw <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_raw_doe, aes(x=PC1, y=PC2, color=run.order, shape=type) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_raw[1], "%)"),
      #         y=paste0("PC2 (", pc12_raw[2], "%)"))
      #)
      plot_ly(pca_table_raw_doe, x=~PC1, y=~PC2, type = "scatter", color=~run.order, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE),
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
    )
    
    ## lfq
    pca_table_lfq_doe <- getPCAtableDOE(table_log2_clean_lfq, doe)
    pc12_lfq <- getPC12Percentage(table_log2_clean_lfq)
    
    output$pca_st_lfq <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_lfq_doe, aes(x=PC1, y=PC2, color=sample.type, shape=type) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_lfq[1], "%)"),
      #         y=paste0("PC2 (", pc12_lfq[2], "%)"))
      #)
      plot_ly(pca_table_lfq_doe, x=~PC1, y=~PC2, type = "scatter", color=~sample.type, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)"), zeroline=FALSE),
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)"), zeroline=FALSE))
    )
    
    output$pca_ro_lfq <- renderPlotly(
      #ggplotly(
      #  ggplot(pca_table_lfq_doe, aes(x=PC1, y=PC2, color=run.order, shape=type) )+ 
      #    geom_point(alpha=0.5, size=3)+
      #    theme_bw()+
      #    labs(x=paste0("PC1 (", pc12_lfq[1], "%)"),
      #         y=paste0("PC2 (", pc12_lfq[2], "%)"))
      #)
      plot_ly(pca_table_lfq_doe, x=~PC1, y=~PC2, type = "scatter", color=~run.order, symbol=~type , text=~sample, marker=list(size=16, opacity=0.5)) %>% 
        layout(xaxis=list(title=paste0("PC1 (", pc12_raw[1], "%)")),
               yaxis=list(title=paste0("PC1 (", pc12_raw[2], "%)")))
    )
    
    ################################################################################################
    # contribution of the variance
    ## raw 
    pvca_raw <- getPVCAresult(table_log2_clean_raw, doe)
    
    output$ctbv_raw <- renderPlotly(
      plot_ly(pvca_raw, x=~Variable, y=~Proportion.Variance, type="bar", color=~Variable) %>% 
        layout(yaxis=list(title="Weighted average proportion variance"), xaxis=list(title="Variables"))
    )
    
    ## lfq
    pvca_lfq <- getPVCAresult(table_log2_clean_lfq, doe)
    
    output$ctbv_lfq <- renderPlotly(
      plot_ly(pvca_lfq, x=~Variable, y=~Proportion.Variance, type="bar", color=~Variable) %>% 
        layout(yaxis=list(title="Weighted average proportion variance"), xaxis=list(title="Variables"))
    )
    
    updateProgressBar(session = session, id = "pb2", value =70, total = 100)
    ################################################################################################
    # quantro
    # check the number of sample.type is >= 2
    if(length(unique(doe$sample.type))<2){
      output$qt_notallow <- renderText(
        "The number of 'sample.type' is less than 2, Quantro test can not be applied!"
      )
    }else{
      output$qt_notallow <- renderText(
        ""
      )
      ## raw
      qt_raw <- doQuantro(table_log2_clean_raw, doe)
      output$qt_qt_raw <- renderTable(
        qt_raw
      )
      qt_raw_a <- doQuantroAnova(table_log2_clean_raw, doe)
      output$qt_anov_raw <- renderTable(
        qt_raw_a, striped = TRUE
      )
      
      ## lfq
      qt_lfq <- doQuantro(table_log2_clean_lfq, doe)
      output$qt_qt_lfq <- renderTable(
        qt_lfq
      )
      qt_lfq_a <- doQuantroAnova(table_log2_clean_lfq, doe)
      output$qt_anov_lfq <- renderTable(
        qt_lfq_a, striped = TRUE
      )
      
      # use pg_ints_raw_melt_doe or pg_ints_lfq_melt_doe
      # boxplot
      ## raw
      pg_ints_raw_melt_doe$value = as.numeric(pg_ints_raw_melt_doe$value)
      output$qt_box_raw <- renderPlotly(
        plot_ly(pg_ints_raw_melt_doe %>% arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
          layout(yaxis=list(title="Log2 Intensity "),
                 xaxis=list(title="Sample", 
                            categoryorder= "array",
                            categoryarray=~sample.type))
      )
      
      ## lfq
      pg_ints_lfq_melt_doe$value = as.numeric(pg_ints_lfq_melt_doe$value)
      output$qt_box_lfq <- renderPlotly(
        plot_ly(pg_ints_lfq_melt_doe %>% arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
          layout(yaxis=list(title="Log2 Intensity "),
                 xaxis=list(title="Sample", 
                            categoryorder= "array",
                            categoryarray=~sample.type))
      )
      
      # density
      ## raw
      output$qt_den_raw <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe, aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()
        )
      )
      
      output$qt_den_raw_split <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe, aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()+
            facet_wrap(~ sample.type)
        )
      )
      ## lfq
      output$qt_den_lfq <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe, aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()
        )
      )
      output$qt_den_lfq_split <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe, aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()+
            facet_wrap(~ sample.type)
        )
      )
    }
    
    
    # quantro --- No control
    nc_doe <- doe %>% filter(type=="sample")
    if(length(unique(nc_doe$sample.type))<2){
      output$nc_qt_notallow <- renderText(
        "The number of 'sample.type' after filtering out the control samples is less than 2, Quantro test can not be applied!"
      )
    }else{
      ## raw
      nc_qt_raw <- doQuantroNoControl(table_log2_clean_raw, doe)
      output$nc_qt_qt_raw <- renderTable(
        nc_qt_raw
      )
      nc_qt_raw_a <- doQuantroNoControlAnova(table_log2_clean_raw, doe)
      output$nc_qt_anov_raw <- renderTable(
        nc_qt_raw_a, striped = TRUE
      )
      
      ## lfq
      nc_qt_lfq <- doQuantroNoControl(table_log2_clean_lfq, doe)
      output$nc_qt_qt_lfq <- renderTable(
        nc_qt_lfq
      )
      nc_qt_lfq_a <- doQuantroNoControlAnova(table_log2_clean_lfq, doe)
      output$nc_qt_anov_lfq <- renderTable(
        nc_qt_lfq_a, striped = TRUE
      )
      
      # use pg_ints_raw_melt_doe or pg_ints_lfq_melt_doe
      # boxplot
      ## raw
      pg_ints_raw_melt_doe$value = as.numeric(pg_ints_raw_melt_doe$value)
      output$nc_qt_box_raw <- renderPlotly(
        plot_ly(pg_ints_raw_melt_doe %>% filter(type=="sample") %>%  arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
          layout(yaxis=list(title="Log2 Intensity "),
                 xaxis=list(title="Sample", 
                            categoryorder= "array",
                            categoryarray=~sample.type))
      )
      
      ## lfq
      pg_ints_lfq_melt_doe$value = as.numeric(pg_ints_lfq_melt_doe$value)
      output$nc_qt_box_lfq <- renderPlotly(
        plot_ly(pg_ints_lfq_melt_doe %>% filter(type=="sample") %>% arrange(sample.type), x=~sample, y=~value, type="box", color=~sample.type) %>% 
          layout(yaxis=list(title="Log2 Intensity "),
                 xaxis=list(title="Sample", 
                            categoryorder= "array",
                            categoryarray=~sample.type))
      )
      
      # density
      ## raw
      output$nc_qt_den_raw <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()
        )
      )
      
      output$nc_qt_den_raw_split <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()+
            facet_wrap(~ sample.type)
        )
      )
      ## lfq
      output$nc_qt_den_lfq <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()
        )
      )
      output$nc_qt_den_lfq_split <- renderPlotly(
        ggplotly(
          ggplot(pg_ints_raw_melt_doe %>% filter(type=="sample"), aes(x=value))+
            stat_density(aes(group=sample, color=sample.type), position="identity", geom = "line")+
            labs(x="Log2 Intensity", y="Density")+
            theme_minimal()+
            facet_wrap(~ sample.type)
        )
      )
    }
    
    ################################################################################################
    # end
    
  })
  
  
  
  ###################################################################################
  
  ### major demo process ###
  observeEvent(input$demo, {
    useSweetAlert()
    progressSweetAlert(
      session = session,
      id = "pb2",
      title = "In progress",
      striped = TRUE,
      value = 0
    )
    ### 1 proteinGroups.txt
    updateProgressBar(session = session, id = "pb2", value = 10)
    pg_ex_response()
    updateProgressBar(session = session, id = "pb2", value = 20)
    
    ### 2 peptides.txt
    pt_ex_response()
    updateProgressBar(session = session, id = "pb2", value = 40)
    
    ### 3 msms.txt
    ### 4 summary.txt
    ms_sum_ex_response()
    updateProgressBar(session = session, id = "pb2", value = 50)
    
    ### 5 doe.txt
    updateProgressBar(session = session, id = "pb2", value = 60)
    doe_ex_response()
    updateProgressBar(session = session, id = "pb2", value = 85)
    updateProgressBar(session = session, id = "pb2", value = 100)
    closeSweetAlert(session = session)
    shinyalert(title = "Demo run is finished", text = "Please check the results by clicking the tabs in the sidebar.", type = "success")
  
  
  })
  
  
  ###################################################################################
  ###################
  # generate report #
  ###################
  output$report <- downloadHandler(
    
    filename = function(){
      paste0("QC-MQ_report_", Sys.time(),".html")
    },
    content = function(file){
      shiny::withProgress(
        message = "Generating the report",
        value = 0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          tempReport <- file.path(tempdir(), "report.Rmd")
          file.copy("report.Rmd", tempReport, overwrite = TRUE)
          tempCode <- file.path(tempdir(), "functions.R")
          file.copy("functions.R", tempCode, overwrite = TRUE)
          shiny::incProgress(3/10)
          params <- list(
            pgfile = input$pg_file$datapath,
            ptfile = input$pt_file$datapath,
            msfile = input$ms_file$datapath,
            smfile = input$sum_file$datapath,
            doefile = input$doe_file$datapath,
            doe_sep = input$sep
          )
          shiny::incProgress(4/10)
          rmarkdown::render(tempReport, output_file = file,
                            params = params,
                            envir = new.env(parent = globalenv()))
          
        }
      )
      
    }
  )
  
  
  
}








