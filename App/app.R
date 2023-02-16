# DrugSeq shiny app
#
# THINGS TO FIX:
## find Vas' code for processing L1000 data (highest priority)
## add support for 0-2 metadata groups (t-test, or no test - just ranking)
## allow user to select metadata category (low priority)

# required packages

library(shiny)
library(dplyr)
library(purrrlyr) # necessary for summary statistics
library(ggpubr) # necessary for boxplots
library(data.table)
library(DT)
library(tidyverse) # necessary to bind together dataframes for final datatable
library(grid)
library(gridExtra)


# other requirements - allow larger upload file size, preload L1000 signature data, set not_sel argument
options(shiny.maxRequestSize = 100 * 1024 ^ 2)
not_sel <- "Not Selected"
Drug_MOA_L10002017 <-
  read.csv("Drug_MOA_L10002017.csv")
L1000_signature_data <-
  read.table(file = "Drug_signatures.txt",
             row.names = 1,
             header = TRUE)

# define function to calculate discordance
create_L1000_discordance_table <-
  function(Disease_signature, Metadata_input) {
    L1000_discordance <-
      Metadata_input #use metadata file to define row names
    showModal(
      modalDialog(
        'Calculating drug discordance (this may take a few minutes)',
        footer = NULL
      )
    )
    for (i in row.names(L1000_signature_data)) {
      target_compound <- i
      compound_sig <-
        as.data.frame(t(
          subset(
            L1000_signature_data,
            rownames(L1000_signature_data) == target_compound
          )
        ))
      compound_sig <-
        subset(compound_sig, compound_sig != 0) #subset compound signature
      cmpdGenes <- rownames(compound_sig)
      cmpd_overlap <-
        rownames(Disease_signature)[which(rownames(Disease_signature) %in% cmpdGenes)] #subset this disease data for the genes in target compound response signature
      overlap <- Disease_signature[cmpd_overlap,]
      cmpd_ordered <- as.numeric(as.vector(t(compound_sig)))
      names(cmpd_ordered) <- rownames(compound_sig)
      cmpd_ordered2 <-
        cmpd_ordered[cmpd_overlap] # Character list in brackets orders to fit that character list...
      overlap2 <-
        overlap #set up each tumor as a vector to be able to calculate SC for each
      mylist <- list()
      options(warn = -1)
      for (i in 1:ncol(overlap2)) {
        cori <- cor(cmpd_ordered2, overlap[i], method = "spearman")
        mylist[[i]] <- cori
      }
      options(warn = 0)
      overlap[i]
      df <- do.call("rbind", mylist)
      as.data.frame(df)
      row.names(df) <- colnames(overlap2)
      colnames(df)[1] = paste(target_compound)
      L1000_discordance <-
        merge(L1000_discordance, df, by = "row.names", fill = TRUE)
      row.names(L1000_discordance) <- L1000_discordance$Row.names
      L1000_discordance <- L1000_discordance[-c(1)]
    }
    L1000_discordance_rmzero <-
      L1000_discordance[, colSums(is.na(L1000_discordance)) == 0] # remove any columns with NA values
    L1000_discordance_table <- L1000_discordance_rmzero
    L1000_discordance_table <-
      as.data.frame(t(L1000_discordance_table)) # transpose
    Drug <- rownames(L1000_discordance_table) # save drug names
    L1000_discordance_table$Drug <-
      Drug # make drug name into column
    L1000_discordance_table <-
      as.data.frame(L1000_discordance_table) # convert to data frame to add index
    L1000_discordance_table$index <-
      1:nrow(L1000_discordance_table) # add index column to sort later
    L1000_discordance_table <-
      merge(Drug_MOA_L10002017,
            L1000_discordance_table,
            by = "Drug",
            all.y = TRUE) # add in corrected drug names
    L1000_discordance_table <-
      L1000_discordance_table[order(L1000_discordance_table$index),] # reorder by index
    L1000_discordance_table <-
      subset(L1000_discordance_table, select = -c(Drug, index, MOA)) # remove unneeded columns
    L1000_discordance_table <-
      as.data.frame(L1000_discordance_table)
    L1000_discordance_table[1, 1] <- "meta"
    L1000_discordance_table_display <-
      data.table(L1000_discordance_table)
    t_L1000_discordance_table <- t(L1000_discordance_table)
    colnames(t_L1000_discordance_table) <-
      t_L1000_discordance_table[1, ]
    t_L1000_discordance_table <-
      as.data.frame(t_L1000_discordance_table[-1, ])
    t_L1000_discordance_table_samp <- t_L1000_discordance_table[1:2]
    t_L1000_discordance_table_num <-
      t_L1000_discordance_table[, -c(1:2)]
    t_L1000_discordance_table_num <-
      as.data.frame(apply(t_L1000_discordance_table_num, 2, as.numeric))
    t_L1000_discordance_table_bind <-
      cbind(t_L1000_discordance_table_samp,
            t_L1000_discordance_table_num)
    L1000_discordance_table_boxplot <-
      t_L1000_discordance_table_bind
    L1000_discordance_table_boxplot <-
      data.table(L1000_discordance_table_boxplot, keep.rownames = TRUE)
    removeModal()
    # summary statistics
    # calculate mean
    showModal(modalDialog('Calculating ANOVA with TukeyHSD', footer = NULL))
    summary_mean <-
      L1000_discordance_rmzero %>% slice_rows("meta") %>% dmap(mean)
    summary_mean <- as.data.frame(t(summary_mean))
    colnames(summary_mean) <- summary_mean[1,]
    colnames(summary_mean) <-
      paste(colnames(summary_mean), "_mean_discordance", sep = "")
    summary_mean <- summary_mean[-1, ]
    summary_mean_rownames <- rownames(summary_mean)
    summary_mean <-
      as.data.frame(apply(summary_mean, 2, as.numeric))
    summary_mean$Drug <- summary_mean_rownames
    summary_mean <- summary_mean %>% relocate (Drug, .before = 1)
    # calculate sd
    summary_sd <-
      L1000_discordance_rmzero %>% slice_rows("meta") %>% dmap(sd)
    summary_sd <- as.data.frame(t(summary_sd))
    colnames(summary_sd) <- summary_sd[1,]
    colnames(summary_sd) <-
      paste(colnames(summary_sd), "_stdev_discordance", sep = "")
    summary_sd <- summary_sd[-1, ]
    summary_sd_rownames <- rownames(summary_sd)
    summary_sd <- as.data.frame(apply(summary_sd, 2, as.numeric))
    summary_sd$Drug <- summary_sd_rownames
    summary_sd <- summary_sd %>% relocate (Drug, .before = 1)
    summary_stats_table_list <- list(summary_mean, summary_sd)
    summary_stats_table <-
      as.data.table(summary_stats_table_list %>% reduce(full_join, by = "Drug"))
    summary_stats_table <- data.table(summary_stats_table)
    #for disease signatures of 3 or more metadata subgroups, calculate ANOVA with tukey HSD post hoc test
    #calculate ANOVA and save
    L1000_discordance_no_metatdata <-
      L1000_discordance_rmzero[-c(1)]
    drugs <- colnames(L1000_discordance_no_metatdata)
    formulae <-
      lapply(drugs, function(x)
        as.formula(paste0(x, " ~ meta")))
    ANOVA <-
      lapply(formulae, function(x)
        summary(aov(x, data = L1000_discordance_rmzero)))
    names(ANOVA) <- format(formulae)
    ANOVA_p_value <-
      unlist(lapply(ANOVA, function(x)
        x[[1]]$"Pr(>F)"[1]))
    ANOVA_p_value_df <-
      data.frame(Drug = sub(' ~ meta', '', names(ANOVA_p_value)),
                 ANOVA_pvalue = ANOVA_p_value)
    ANOVA_p_value_table <<- ANOVA_p_value_df
    data.table(ANOVA_p_value_table)
    #calculate TukeyHSD and save
    TukeyHSD <-
      lapply(formulae, function(x)
        TukeyHSD(aov(x, data = L1000_discordance_rmzero)))
    names(TukeyHSD) <- format(formulae)
    list_of_drug_names_with_category <-
      paste0(colnames(L1000_discordance_no_metatdata), " ~ meta")
    TukeyHSD_p_value <- NULL
    TukeyHSD_p_value_list = list()
    for (i in list_of_drug_names_with_category) {
      TukeyHSD_i <- Reduce(full_join, TukeyHSD[[i]])
      colnames(TukeyHSD_i) <- c('diff', 'lwr', 'upr', i)
      TukeyHSD_i_t <- t(TukeyHSD_i)
      TukeyHSD_p_value_list[[i]] <- TukeyHSD_i_t
    }
    TukeyHSD_p_value = do.call(rbind, TukeyHSD_p_value_list)
    row_names_df_to_remove <- c("diff", "lwr", "upr")
    TukeyHSD_p_value <-
      TukeyHSD_p_value[!(row.names(TukeyHSD_p_value) %in% row_names_df_to_remove),]
    TukeyHSD_p_value <- as.data.frame(TukeyHSD_p_value)
    TukeyHSD_p_value_rownames <-
      gsub("\\ ~ meta*", "", row.names(TukeyHSD_p_value))
    colnames(TukeyHSD_p_value) <-
      paste(colnames(TukeyHSD_p_value), "_TukeyHSD", sep = "")
    TukeyHSD_p_value$Drug <- TukeyHSD_p_value_rownames
    TukeyHSD_p_value <-
      TukeyHSD_p_value %>% relocate (Drug, .before = 1)
    TukeyHSD_p_value_table <<- TukeyHSD_p_value
    data.table(TukeyHSD_p_value_table)
    # bind together everything to make drug table
    ANOVA_p_value_to_bind <- as.data.frame(ANOVA_p_value_df)
    TukeyHSD_p_value_to_bind <- as.data.frame(TukeyHSD_p_value)
    Summary_stats_to_bind <- as.data.frame(summary_mean)
    drug_table_list <-
      list(ANOVA_p_value_to_bind,
           TukeyHSD_p_value_to_bind,
           Summary_stats_to_bind)
    Drug_table <-
      drug_table_list %>% reduce(full_join, by = "Drug")
    Drug_table <- merge(Drug_MOA_L10002017, Drug_table, by = "Drug")
    Drug_table <- subset(Drug_table, select = -c(Drug))
    data.table(Drug_table)
    removeModal()
    list(
      Drug_table = Drug_table,
      L1000_discordance = L1000_discordance_table_display,
      L1000_discordance_boxplot = L1000_discordance_table_boxplot
    )
  }

# define UI
ui <- fluidPage(navbarPage(
  title = "DrugSeq",
  tabPanel(title = "About",
           titlePanel("About DrugSeq platform"),
           hr(),
           img(src = 'DrugSeq_schematic.jpg', width = '90%')),
  tabPanel(
    title = "DrugSeq: Drug Discordance Analysis",
    titlePanel("Drug Discordance Analysis"),
    sidebarLayout(
      sidebarPanel(
        h3("Step 1: Upload files"),
        a(
          href = "DrugSeq_templates.csv",
          "Download template files",
          download = NA,
          target = "_blank"
        ),
        br(),
        br(),
        fileInput(
          "disease_signature",
          "Upload disease signature csv file",
          accept =
            ".csv"
        ),
        fileInput("metadata", "Upload metadata csv file", accept =
                    ".csv"),
        selectInput(
          "L1000_dataset",
          "Select L1000 dataset (coming soon)",
          choices = c("2017", "2021"),
          selected = FALSE,
          multiple = FALSE
        ),
        hr(),
        h3("Step 2a: Calculate drug discordance"),
        h3("Step 2b: Calculate subgroup statistics"),
        actionButton("run_button", "Run Analysis", icon =
                       icon("play")),
        hr(),
        h3("Step 3: Select drugs to plot"),
        selectizeInput(
          "drugs_to_plot",
          "Select L1000 drug to plot",
          choices = NULL,
          selected = FALSE,
          multiple = TRUE
        )
      ),
      mainPanel(tabsetPanel(
        tabPanel(
          title = "Drug Discordance",
          br(),
          downloadButton(
            'download_drug_discordance_table',
            "Download drug discordance table"
          ),
          hr(),
          fluidRow(column(
            width = 12, DT::dataTableOutput("drug_discordance_table")
          ))
        ),
        tabPanel(
          title = "Statistics",
          br(),
          downloadButton('download_stats_table', "Download combined statistics table"),
          hr(),
          fluidRow(column(
            width = 12, DT::dataTableOutput("drug_stats_table")
          ))
        ),
        tabPanel(
          title = "Box Plot",
          br(),
          downloadButton(
            'pdf',
            "Download box plots"
          ),
          hr(),
          uiOutput("plots")
        )
      ))
    ),
  )
))

# define server function
server <- function(input, output, session) {
  L1000_signature_data <-
    read.table(file = "Drug_signatures.txt",
               row.names = 1,
               header = TRUE)
  Drug_signatures <- L1000_signature_data
  Disease_signature_input <- reactive({
    if (is.null(input$disease_signature))
    {
      return()
    }
    else {
      input$disease_signature$datapath
    }
  })
  Disease_signatures <-
    reactive({
      read.csv(file = Disease_signature_input(),
               row.names = 1,
               header = TRUE)
    })
  Metadata_input <- reactive({
    if (is.null(input$metadata))
    {
      return()
    }
    else {
      input$metadata$datapath
    }
  })
  Metadata <-
    reactive({
      read.csv(file = Metadata_input(),
               row.names = 1,
               header = TRUE)
    })
  observeEvent(Metadata_input(), {
    choices <- c(not_sel, colnames(Metadata_input()))
    updateSelectInput(inputId = "metadata_category", choices = choices)
  })
  
  L1000_discordance_table <- eventReactive(input$run_button, {
    req(input$disease_signature, input$metadata)
    create_L1000_discordance_table(Disease_signatures(), Metadata())
  })
  output$drug_stats_table <- DT::renderDataTable({
    DT::datatable(
      L1000_discordance_table()$Drug_table,
      options = list(autoWidth = TRUE),
      filter = list(position = 'top', clear = FALSE)
    )
  })
  output$download_stats_table <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "DrugSeq_combined_statistics.csv", sep = "_")
    },
    content = function(file) {
      write.csv(L1000_discordance_table()$Drug_table, file)
    }
  )
  output$drug_discordance_table <- DT::renderDataTable({
    DT::datatable(L1000_discordance_table()$L1000_discordance,
                  options = list(autoWidth = TRUE),)
  })
  output$download_drug_discordance_table <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "DrugSeq_drug_discordance.csv", sep = "_")
    },
    content = function(file) {
      write.csv(L1000_discordance_table()$L1000_discordance, file)
    }
  )
  observe({
    updateSelectizeInput(
      session,
      "drugs_to_plot",
      choices = unique(L1000_discordance_table()$Drug_table$Drug_name),
      server = TRUE
    )
  })
  observe({
    req(input$drugs_to_plot)
    lapply(input$drugs_to_plot, function(par) {
      p <-
        ggboxplot(
          L1000_discordance_table()$L1000_discordance_boxplot,
          x = "meta",
          y = par,
          color = "gray30",
          fill = "meta",
          palette = "rainbow",
          group_by = "meta",
          ylab = paste(par, " discordance"),
          xlab = "metadata group",
          width = .6,
          outlier.shape = NA,
          bxp.errorbar = TRUE,
          bxp.errorbar.width = 0.2,
          legend = 0
        ) +
        ggtitle(par) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept =
                     0,
                   linetype = "dashed",
                   color = "gray8") +
        theme(
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14)
        )
      output[[paste("plot", par, sep = "_")]] <- renderPlot({
        p
      },
      width = 380,
      height = 350)
    })
  })
  plot_output_list <- reactiveValues()
  output$plots <- renderUI({
    req(input$drugs_to_plot)
    plot_output_list <- lapply(input$drugs_to_plot, function(par) {
      plotname <- paste("plot", par, sep = "_")
      plotOutput(plotname, height = '250px', inline = TRUE)
    })
    do.call(tagList, plot_output_list)
    
  })
  output$pdf <- downloadHandler(
    filename = function(){
      paste(Sys.Date(), "DrugSeq_box_plots.pdf", sep = "_")
    },
    content = function(file) {
      plot_list <- lapply(input$drugs_to_plot, function(par) {
        ggboxplot(
          L1000_discordance_table()$L1000_discordance_boxplot,
          x = "meta",
          y = par,
          color = "gray30",
          fill = "meta",
          palette = "rainbow",
          group_by = "meta",
          ylab = paste(par, " discordance"),
          xlab = "metadata group",
          width = .6,
          outlier.shape = NA,
          bxp.errorbar = TRUE,
          bxp.errorbar.width = 0.2,
          legend = 0
        ) +
          ggtitle(par) +
          geom_point(alpha = 0.5) +
          geom_hline(yintercept =
                       0,
                     linetype = "dashed",
                     color = "gray8") +
          theme(
            plot.title = element_text(size = 20),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14)
          )
      })
      pdf(file)
      arrangeGrob(for (i in 1:length(plot_list)) {
        print(plot_list[[i]])
      },
      ncol = length(plot_list))
      dev.off()
    }
  )
}
# run the application
shinyApp(ui = ui, server = server)
