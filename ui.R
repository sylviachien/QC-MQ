# app_10


###################################################################################
# load libraries
library(shiny)
library(shinyBS)
library(shinyjs)
library(shinydashboard)
library(shinyFiles)
library(plotly)
library(shinyalert)
library(shinyWidgets)
library(BiocManager)
#options(repos = BiocManager::repositories())

###################################################################################
useShinyalert()
useSweetAlert()

# add logo along with the title in the header
head <- tags$a(href="https://report.pri.bms.com/QC-MQ/",
               tags$img(src="logo1.png", height=80, width=80), 
               "Quality Control Dashboard for MaxQuant Proteomics Data", target="_blank")



ui <- dashboardPage(
  
  dashboardHeader(
    # setup the height of header
    tags$li(
      class="dropdown",
      tags$style(".main-header {max-height: 80px}"),
      tags$style(".main-header .logo {height: 80px}")
    ),
    title = head,
    titleWidth = 800
    
  ),
  ###################################################################################
  
  dashboardSidebar(
    # adjust the sidebar height for the header height
    tags$style(".left-side, .main-sidebar {padding-top: 100px}"),
    
    sidebarMenu(
      menuItem("Upload MaxQuant/DOE files", tabName = "upload", icon = icon("upload")),
      menuItem("Data Structure", tabName = "datastructure", icon = icon("datastructure")),
      menuItem("Experimental Checkpoints", tabName = "expCheck", icon = icon("expCheck")),
      menuItem("Contaminants", tabName = "contaminants", icon = icon("contaminants")),
      menuItem("Protein/Peptide Counts", tabName = "pp_count", icon = icon("pp_count")),
      menuItem("Protein/Peptide Intensity", tabName = "pp_intensity", icon = icon("pp_intensity")),
      menuItem("Coefficient of Variation", tabName = "cv", icon = icon("cv")),
      menuItem("Top 20 Proteins", tabName = "top20", icon = icon("top20")),
      menuItem("Principal Component Analysis", tabName = "pca", icon = icon("pca")),
      menuItem("Contribution of Variance", tabName = "contriVariance", icon = icon("contriVariance")),
      menuItem("Quantro Global Shift Test", tabName = "quantrotest", icon = icon("quantrotest")),
      menuItem("Download Datasets", tabName = "downloaddataset", icon = icon("downloaddataset"))
    )
    
  ),
  ###################################################################################
  
  dashboardBody(
    
    useShinyjs(),
    useShinyalert(),
    useSweetAlert(),
    
    # add reference to CSS file in the www folder from the app directory
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    
    # set up different pages
    tabItems(
      # upload page
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "Instructions",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            h4(strong(span("Run demo QC", style = "color:red"))),
            HTML(
              paste(
                h5(HTML('&emsp;'), strong("   Steps:")),
                h5(HTML('&emsp;'), strong("   1. Click 'Run Demo' button.")),
                h5(HTML('&emsp;'), strong("   2. After the success message poping up, browse the results by clicking the tabs in the side menu."))
              )
            ),
            br(),
            h4(strong(span("Run your own data", style = "color:red"))),
            HTML(
              paste(
                h5(HTML('&emsp;'), strong("   Steps:")),
                h5(HTML('&emsp;'), strong("   1. Have the output files of MaxQuant(MQ) ready. They are in the combined/txt folder. The minimum uploading requirement to run this QC is the 'proteinGroups.txt'.")),
                h5(HTML('&emsp;'), strong("   2. If applicable, please upload the other MQ output files: 'peptides.txt', 'msms.txt', and 'summary.txt'.")),
                h5(HTML('&emsp;'), strong("   3. Prepare the DOE file (ex: doe.txt) in the suggested format:"))
              )
            ),
            HTML(
              paste(
                h5(HTML('&emsp;'), HTML('&emsp;'), "- 'Sample.id', 'Run.order', 'Type', and 'Sample.type' columns are required in the DOE."),
                h5(HTML('&emsp;'), HTML('&emsp;'), "- Only two categories: 'control' and 'sample' in the 'Type' column."),
                h5(HTML('&emsp;'), HTML('&emsp;'), "- Avoid naming the sample (Sample.id) starting with a number."),
                h5(HTML('&emsp;'), strong("   4. Upload the DOE file, wrong format of the DOE would stop running QC process.")),
                h5(HTML('&emsp;'), strong("   5. Click 'Start QC' button.")),
                h5(HTML('&emsp;'), strong("   6. After the success message poping up, browse the results by clicking the tabs in the side menu."))
                
              )
            ),
            br(),
            h4(strong(span("Generate report", style = "color:red"))),
            HTML(
              paste(
                h5(HTML('&emsp;'), strong("   Steps:")),
                h5(HTML('&emsp;'), strong("   1. Upload the files as described in 'Run your own data'.")),
                h5(HTML('&emsp;'), strong("   2. Click 'Download the report' button.")),
                h5(HTML('&emsp;'), strong("   3. QC-MQ report is downloaded."))
                
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Upload file(s)",
            status = "primary",
            solidHeader = TRUE,
            fileInput(
              inputId = "pg_file",
              label = "Choose the 'proteinGroups.txt' to upload"
            ),
            bsAlert("alert1"),
            fileInput(
              inputId = "pt_file",
              label = "Choose the 'peptides.txt' to upload (Optional)"
            ),
            bsAlert("alert2"),
            fileInput(
              inputId = "ms_file",
              label = "Choose the 'msms.txt' to upload (Optional)"
            ),
            bsAlert("alert3"),
            fileInput(
              inputId = "sum_file",
              label = "Choose the 'summary.txt' to upload (Optional)"
            ),
            bsAlert("alert4"),
            radioButtons(inputId = "sep",
                         label = "Separator of DOE/doe.txt",
                         choices = c(Tab='\t',Comma=',',Semicolon=';',Space=' '), selected = '\t'),
            fileInput(inputId = "doe_file",
                      label = "Choose the 'doe.txt' to upload (Optional)",
                      accept = "txt/tab-separated-values"),
            bsAlert("alert5a"),
            bsAlert("alert5b"),
            bsAlert("alert5c"),
            bsAlert("alert5d"),
            bsAlert("alert10"),
            helpText("The doe.txt must have at least the following columns: 'Sample.id', 'Run.order', 'Type', and 'Sample.type'. The 'Type' column contains only 'control' and 'sample' categories."),
            br(),
            h5(""),
            actionButton("doe_demo_table", label = "Example of DOE"),
            tableOutput("doe_example")
            
          ),
          box(
            title = "Input Parameters",
            status = "primary",
            solidHeader = TRUE,
            radioButtons("LFQRun",
                         label = "Did you run LFQ in MaxQuant?",
                         choices = list("Yes"="Y", "No"="N"), selected = 'N'),
            radioButtons("spikeinRun",
                         label = "Did you include spike-in (DigestIF, P00000) in the experiment?",
                         choices = list("Yes"="Y", "No"="N"), selected = 'N'),
            #textInput("custom_spikeinRun",
            #          label = "Please specify the 'Protein.IDs' of the protein group you would like to check as a spike-in: (only one protein group is allowed each time) (Optional)"),
            #bsAlert("alert11"),
            textInput("prefix",
                      label = "Output prefix: (Required)"),
            textInput("exp_title",
                      label = "Title of experiment: (Optional)"),
            
            br(),
            useSweetAlert(),
            actionButton("start", label = "Start QC"),
            #progressBar(
            #  id = "pb1",
            #  status = "primary",
            #  striped = TRUE,
            #  value = 0,
            #  display_pct = FALSE
            #),
            bsAlert("alert1a"), # pg required
            #bsAlert("alert6"),  # prefix required
            bsAlert("alert2a"),
            bsAlert("alert3a"),
            bsAlert("alert4a"),
            bsAlert("alert5e"),
            bsAlert("alert7"),
            br(),
            br(),
            actionButton("demo", label = "Run Demo"),
            #shinyWidgets::progressBar(
            #  id = "pb2",
            #  status = "info",
            #  striped = TRUE,
            #  value = 0,
            #  display_pct = FALSE
            #),
            br(),
            br(),
            br(),
            downloadButton("report", label = "Download the report")
            #progressBar(
            #  id = "pb3",
            #  status = "info",
            #  striped = TRUE,
            #  value = 0,
            #  display_pct = FALSE
            #)
          )
        )
      ),
      # one page, one table, one plot
      tabItem(
        tabName = "datastructure",
        h4(strong("DOE is required for this analysis.")),
        h4(textOutput("doe_descrp")),
        br(),
        column(width = 3, offset = 1, h4(tableOutput("datastructure_table"))),
        column(width = 6, plotlyOutput("datastructure_plot"))
      ),
      # two tabs, two boxes in each tab
      tabItem(
        tabName = "expCheck",
        tabsetPanel(
          tabPanel(
            "Digestion Efficiency",
            br(),
            fluidRow(
              box(
                title = "Occurrence of Cysteine/Missed Cleavage", 
                status = "info",
                width =12,
                br(),
                tableOutput("expck_table"),
                h4(strong("peptides.txt is required for this analysis.")),
                br(),
                br(),
                p("The occurence of cysteine should be over 10%."),
                p("It could be consider exceeding the standard (~20-25%), if the rate of missed-cleavage is above 40%.")
              )
            ),
            fluidRow(
              box(
                title = "Missed cleavages per sample",
                status = "info",
                width = 12,
                h4(strong("msms.txt and summary.txt are required for this analysis.")),
                plotlyOutput("msclvg_persample")
              )
            )
          ),
          tabPanel(
            "Spike-in(P00000 DIGESTIF)",
            br(),
            h4(strong("Spike-in protein must be included for this analysis.")),
            
            fluidRow(
              box(
                title = "Intensity of Spikein",
                status = "info",
                width = 12,
                tabBox(
                  title = "",
                  width = 12,
                  tabPanel(
                    "Raw Intensity",
                    br(),
                    plotlyOutput("spike_plot_raw")
                  ),
                  tabPanel(
                    "LFQ Intensity",
                    h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                    plotlyOutput("spike_plot_lfq")
                  )
                )
                
                
              ),
              box(
                title = "Samples with spike-in intensity falling outside 2SD range",
                status = "info",
                width = 12,
                tabBox(
                  title = "",
                  width = 12,
                  tabPanel(
                    "Raw Intensity",
                    br(),
                    h5("Sample with spike-in intensity less than mean - 2SD"),
                    tableOutput("spikelow_raw"),
                    br(),
                    h5("Sample with spike-in intensity more than mean + 2SD"),
                    tableOutput("spikehigh_raw")
                    
                  ),
                  tabPanel(
                    "LFQ Intensity",
                    br(),
                    h5("Sample with spike-in intensity less than mean - 2SD"),
                    tableOutput("spikelow_lfq"),
                    br(),
                    h5("Sample with spike-in intensity more than mean + 2SD"),
                    tableOutput("spikehigh_lfq")
                  )
                  
                )
              )
            )
          )
        )
      ),
      # two tabs, two tabs in each tab, 
      tabItem(
        tabName = "contaminants",
        br(),
        h5("The proteinGroups are divided into different contaminant categories by their origins (Bovine) and the uniqueness while mapping to protein IDs in this analysis: "),
        HTML(
          paste(
            h5(HTML('&emsp;'), span("- Non-Contaminant", style = "color:orange")),
            h5(HTML('&emsp;'), span("- Contaminant.Bovine.Unique", style = "color:orange")),
            h5(HTML('&emsp;'), span("- Contaminant.Bovine.NonUnique", style = "color:orange")),
            h5(HTML('&emsp;'), span("- Contaminant.NonBovine.Unique", style = "color:orange")),
            h5(HTML('&emsp;'), span("- Contaminant.NonBovine.NonUnique", style = "color:orange")),
            h5(HTML('&emsp;'), span("- Trypsin", style = "color:orange")),
            h5(HTML('&emsp;'), span("- LGB(Bovine)", style = "color:orange"))
          )
        ),
        br(),
        h5("The 'Unique' mapping is determined when the proteinGroup only contains the protein.id(s) that mapped to the defined contaminant(s)."),
        h5("The 'Unique' contaminants will be removed from the data for the most of QC analyses and the clean dataset can be downloaded in the 'Download Datasets' page."),
        h5("Among the contaminants, we also particularly examine the levels of 'Trpsin' and 'Bovine LGB (Beta-Lactoglobulin)' in the Contaminant analysis."),
        br(),
        br(),
        tabsetPanel(
          tabPanel(
            "Overall by Individual Sample",
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('contam_raw'),
                  plotlyOutput('contam_raw_log2')
                ),
                tabPanel(
                  "LFQ Intensity",
                  br(),
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('contam_lfq'),
                  plotlyOutput('contam_lfq_log2')
                )
              )
            )
          ),
          tabPanel(
            "By Sample Type",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('contam_raw_st')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('contam_lfq_st')
                )
              )
            )
          )
        )
      ),
      # two tabs, three tabs in each tab, two tabs in each tab
      tabItem(
        tabName = "pp_count",
        tabsetPanel(
          tabPanel(
            "Protein Count",
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "By Individual Sample",
                  fluidRow(
                    tabBox(
                      title = "",
                      width = 12,
                      tabPanel(
                        "Raw Intensity",
                        br(),
                        plotlyOutput('pg_count_raw')
                      ),
                      tabPanel(
                        "LFQ Intensity",
                        br(),
                        h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                        plotlyOutput('pg_count_lfq')
                      )
                    )
                  )
                ),
                tabPanel(
                  "By Sample Type",
                  br(),
                  h4(strong("DOE is required for this analysis.")),
                  fluidRow(
                    tabBox(
                      title = "",
                      width = 12,
                      tabPanel(
                        "Raw Intensity",
                        br(),
                        plotlyOutput("pg_count_st_raw")
                      ),
                      tabPanel(
                        "LFQ Intensity",
                        br(),
                        h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                        plotlyOutput('pg_count_st_lfq')
                      )
                    )
                  )
                  
                ),
                tabPanel(
                  "By Type (Control/Experimental Sample)",
                  br(),
                  h4(strong("DOE is required for this analysis.")),
                  fluidRow(
                    tabBox(
                      title = "",
                      width = 12,
                      tabPanel(
                        "Raw Intensity",
                        br(),
                        plotOutput('pg_count_t_raw')
                      ),
                      tabPanel(
                        "LFQ Intensity",
                        br(),
                        h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                        plotOutput('pg_count_t_lfq')
                      )
                    )
                  )
                )
              )
            )
          ),
          tabPanel(
            "Peptide Count",
            h4(strong("peptides.txt is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "By Individual Sample",
                  br(),
                  plotlyOutput('pt_count')
                ),
                tabPanel(
                  "By Sample Type",
                  h4(strong("DOE is required for this analysis.")),
                  plotlyOutput('pt_count_st')
                ),
                tabPanel(
                  "By Type (Control/Experimental Sample)",
                  h4(strong("DOE is required for this analysis.")),
                  plotlyOutput('pt_count_t')
                )
              )
            )
          )
        )
      ),
      # five tabs,two-three tabs in each tab
      tabItem(
        tabName = "pp_intensity",
        tabsetPanel(
          tabPanel(
            "Overall by Individual Sample",
            fluidRow(
              tabBox(
                title = "",
                width = 12, 
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pg_int_raw')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pg_int_lfq')
                ),
                tabPanel(
                  "Peptide Intensity",
                  h4(strong("peptides.txt is required for this analysis.")),
                  plotlyOutput('pt_int')
                )
              )
            )
          ),
          tabPanel(
            "Distribution of protein group intensities per sample",
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pg_int_vio_raw')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pg_int_vio_lfq')
                )
              )
            ),
            fluidRow(
              h4(strong("DOE is required for this analysis.")),
              tabBox(
                title = "Violin Plot by Run order",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pg_int_vio_raw_doe')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pg_int_vio_lfq_doe')
                )
              )
            )
          ),  
          tabPanel(
            "Mean Proportion Plot by type",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pg_int_mpp_t_raw')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pg_int_mpp_t_lfq')
                )
              )
            )
          ),
          tabPanel(
            "Cumulative Intensity Distribution for fraction of sample by sample type",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                height = 1200,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotOutput('pg_int_cdf_st_raw', height = 1000)
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotOutput('pg_int_cdf_st_lfq', height = 1000)
                )
              )
            )
          ),
          tabPanel(
            "Raw vs LFQ",
            h4(strong("LFQ run in MaxQuant is required for this analysis.")),
            fluidRow(
              box(
                title = "Overall Sample", 
                status = "info",
                width =12,
                br(),
                plotlyOutput('pg_int_rawlfq', height = 400, width = 400)
              )
            ),
            fluidRow(
              box(
                title = "By Sample Type", 
                status = "info",
                width =12,
                h4(strong("DOE is required for this analysis.")),
                plotlyOutput('pg_int_st_rawlfq', height = 1000)
              )
            )
          )
        )
      ),
      # three tabs, two tabs in each tab
      tabItem(
        tabName = "cv",
        tabsetPanel(
          tabPanel(
            "Overall",
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pg_int_cv_raw')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pg_int_cv_lfq')
                )
              )
            )
          ),
          tabPanel(
            "By Sample Type",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pg_int_cv_st_raw'),
                  plotlyOutput('pg_int_cv_st_raw_split')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pg_int_cv_st_lfq'),
                  plotlyOutput('pg_int_cv_st_lfq_split')
                )
              )
            )
          ),
          tabPanel(
            "By Type (Control/Experimental Sample)",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pg_int_cv_t_raw'),
                  plotlyOutput('pg_int_cv_t_raw_split')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pg_int_cv_t_lfq'),
                  plotlyOutput('pg_int_cv_t_lfq_split')
                )
              )
            )
          )
        )
      ),
      # three tabs, two tabs in each tab
      tabItem(
        tabName = "top20",
        h4("Full tables can be downloaded in the 'Download Datasets' page."),
        tabsetPanel(
          tabPanel(
            "Overall",
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  tableOutput('top20_raw')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  tableOutput('top20_lfq')
                )
              )
            )
          ),
          tabPanel(
            "By Sample Type",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  uiOutput("sampletype_raw"),
                  br(),
                  br(),
                  tableOutput('top20_st_raw')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  br(),
                  uiOutput("sampletype_lfq"),
                  br(),
                  br(),
                  tableOutput('top20_st_lfq')
                )
              )
            )
          ),
          tabPanel(
            "By Type (Control/Experimental Sample)",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  uiOutput("type_raw"),
                  br(),
                  br(),
                  tableOutput('top20_t_raw')
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  br(),
                  uiOutput("type_lfq"),
                  br(),
                  br(),
                  tableOutput('top20_t_lfq')
                )
              )
            )
          )
        )
      ),
      # three tabs, two tabs in each tab
      tabItem(
        tabName = "pca",
        tabsetPanel(
          tabPanel(
            "By Individual Sample",
            fluidRow(
              tabBox(
                title = "",
                width = 12, 
                height = 1000,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pca_raw', height = 700)
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pca_lfq', height = 700)
                )
              )
            )
          ),
          tabPanel(
            "By Sample Type",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                height = 1000,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pca_st_raw', height = 700)
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pca_st_lfq', height = 700)
                )
              )
            )
          ),
          tabPanel(
            "By Run Order",
            h4(strong("DOE is required for this analysis.")),
            fluidRow(
              tabBox(
                title = "",
                width = 12,
                height = 1000,
                tabPanel(
                  "Raw Intensity",
                  br(),
                  plotlyOutput('pca_ro_raw', height = 700)
                ),
                tabPanel(
                  "LFQ Intensity",
                  h4(strong("LFQ run in MaxQuant is required for this analysis.")),
                  plotlyOutput('pca_ro_lfq', height = 700)
                )
              )
            )
          )
        )
      ),
      # two tabs
      tabItem(
        tabName = "contriVariance",
        br(),
        h4(strong("DOE is required for this analysis.")),
        tabsetPanel(
          tabPanel(
            "Raw Intensity",
            br(),
            plotlyOutput('ctbv_raw')
          ),
          tabPanel(
            "LFQ Intensity",
            h4(strong("LFQ run in MaxQuant is required for this analysis.")),
            plotlyOutput('ctbv_lfq')
          )
        )
      ),
      # two tabs, four boxes in each tab
      tabItem(
        tabName = "quantrotest",
        br(),
        h4(strong("DOE is required for this analysis.")),
        p("The global differences between groups is evaluated to assess whether the global normalization methods should be applied.
          Two tests: ANOVA and QuantroStat are included in this analysis."),
        p("The ANOVA shows if the medians of the distributions are different across groups."),
        p("The 'quantroStat' is the ratio of variance between groups to variance within groups.
          The global adjusment for normalization of the data may not be appropriate if the variability between groups is sufficiently larger than that with groups"),
        p("Here we perform the test on both raw and LFQ data and 1000 permutation test are applied to assess the statistical significance of the test statistic. 
          'quantroPvalPerm' is the p-value associated with the proportion of times the test statistics resulting from the permuted samples were larger than quantroStat"),
        a("Reference: Hicks SC,Irizarry RA, quantro: a data-driven approach to guide the choice of an appropriate, Genome Biol. 2015 Jun 4;16:117.", href="https://www.ncbi.nlm.nih.gov/pubmed/?term=quantro%2C+bioinformatics"),
        br(),
        h4(strong(span(textOutput('qt_notallow'), style="color:red"))),
        h4(strong(span(textOutput('nc_qt_notallow'), style="color:red"))),
        tabsetPanel(
          tabPanel(
            "Raw Intensity",
            br(),
            fluidRow(
              box(
                title = "ANOVA Test",
                status = "info",
                tableOutput('qt_anov_raw')
              ),
              box(
                title = "Quantro Test",
                status = "info",
                tableOutput('qt_qt_raw')
              )
            ),
            fluidRow(
              box(
                title = "Intensity Plot",
                status = "info",
                width = 12,
                plotlyOutput('qt_box_raw')
              )
            ),
            fluidRow(
              box(
                title = "Density Plot",
                status = "info",
                width = 12,
                plotlyOutput('qt_den_raw'),
                plotlyOutput('qt_den_raw_split')
              )
            )
          ),
          tabPanel(
            "LFQ Intensity",
            h4(strong("LFQ run in MaxQuant is required for this analysis.")),
            fluidRow(
              box(
                title = "ANOVA Test",
                status = "info",
                tableOutput('qt_anov_lfq')
              ),
              box(
                title = "Quantro Test",
                status = "info",
                tableOutput('qt_qt_lfq')
              )
            ),
            fluidRow(
              box(
                title = "Intensity Plot",
                status = "info",
                width = 12,
                plotlyOutput('qt_box_lfq')
              )
            ),
            fluidRow(
              box(
                title = "Density Plot",
                status = "info",
                width = 12,
                plotlyOutput('qt_den_lfq'),
                plotlyOutput('qt_den_lfq_split')
              )
            )
          ),
          tabPanel(
            "No 'Control'- Raw Intensity",
            br(),
            fluidRow(
              box(
                title = "ANOVA Test",
                status = "info",
                tableOutput('nc_qt_anov_raw')
              ),
              box(
                title = "Quantro Test",
                status = "info",
                tableOutput('nc_qt_qt_raw')
              )
            ),
            fluidRow(
              box(
                title = "Intensity Plot",
                status = "info",
                width = 12,
                plotlyOutput('nc_qt_box_raw')
              )
            ),
            fluidRow(
              box(
                title = "Density Plot",
                status = "info",
                width = 12,
                plotlyOutput('nc_qt_den_raw'),
                plotlyOutput('nc_qt_den_raw_split')
              )
            )
          ),
          tabPanel(
            "No 'Control'- LFQ Intensity",
            h4(strong("LFQ run in MaxQuant is required for this analysis.")),
            fluidRow(
              box(
                title = "ANOVA Test",
                status = "info",
                tableOutput('nc_qt_anov_lfq')
              ),
              box(
                title = "Quantro Test",
                status = "info",
                tableOutput('nc_qt_qt_lfq')
              )
            ),
            fluidRow(
              box(
                title = "Intensity Plot",
                status = "info",
                width = 12,
                plotlyOutput('nc_qt_box_lfq')
              )
            ),
            fluidRow(
              box(
                title = "Density Plot",
                status = "info",
                width = 12,
                plotlyOutput('nc_qt_den_lfq'),
                plotlyOutput('nc_qt_den_lfq_split')
              )
            )
          )
          
        )
      ),
      # one page
      tabItem(
        tabName = "downloaddataset",
        fluidRow(
          box(
            title = "Descriptions of Datasets",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            HTML(
              paste(
                h5(strong(span("Log2_Raw/LFQ_data", style = "color:purple"))),
                p(HTML('&emsp;'), ("The log2 value of raw/lfq intensities of proteinGroups for each sample. The 'Reverse' and 'Only.identified.by.sites' results are filtered in this dataset.")),
                h5(strong(span("Log2_Raw/LFQ_cleaned_data", style = "color:purple"))),
                p(HTML('&emsp;'), ("The log2 value of raw/lfq intensities of proteinGroups for each sample. Besides the 'Reverse' and 'Only.identified.by.sites', the proteinGroups that uniquely mapped to the defined 'Potential Contaminant/Contaminant' are removed in this dataset.")),
                h5(strong(span("Log2_Raw/LFQ_Contaminant_table", style = "color:purple"))),
                p(HTML('&emsp;'), ("The log2 value of raw/lfq intensities of proteinGroups that are defined as 'Potential Contaminant/Contaminant' for each sample.")),
                h5(strong(span("Mean_Log2_Raw/LFQ_data", style = "color:purple"))),
                p(HTML('&emsp;'), ("The log2 mean value of raw/lfq intensities for each proteinGroups across all the samples.")),
                
                h5(strong(span("Mean_Log2_Raw/LFQ_bySampleType_data", style = "color:purple"))),
                p(HTML('&emsp;'), span("Only available when the DOE is provided.", style = "color:blue"), ("The log2 mean value of raw/lfq intensities for each proteinGroups across the samples within the same sample type groups.")),
                h5(strong(span("Mean_Log2_Raw/LFQ_byType_data", style = "color:purple"))),
                p(HTML('&emsp;'), span("Only available when the DOE is provided.", style = "color:blue"), ("The log2 mean value of raw/lfq intensities for each proteinGroups across the samples within the 'control' and 'sample' groups."))
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Raw Data",
            status = "primary",
            solidHeader = TRUE,
            h4(),
            downloadButton("dtable_log2_raw", "Log2_Raw_data"),
            br(),
            br(),
            downloadButton("dtable_log2_clean_raw", "Log2_Raw_cleaned_data"),
            br(),
            br(),
            downloadButton("dtable_contam_log2_raw", "Log2_Raw_Contaminant_table"),
            br(),
            br(),
            downloadButton("dannotated_mean_table_raw", "Mean_Log2_Raw_data"),
            br(),
            br(),
            downloadButton("dwide_sampletype_sort_table_raw", "Mean_Log2_Raw_bySampleTye_data"),
            br(),
            br(),
            downloadButton("dwide_type_sort_table_raw", "Mean_Log2_Raw_byType_data")
          ),
          box(
            title = "LFQ Data",
            status = "primary",
            solidHeader = TRUE,
            h4(strong("LFQ run in MaxQuant is required for this analysis.")),
            downloadButton("dtable_log2_lfq", "Log2_LFQ_data"),
            br(),
            br(),
            downloadButton("dtable_log2_clean_lfq", "Log2_LFQ_cleaned_data"),
            br(),
            br(),
            downloadButton("dtable_contam_log2_lfq", "Log2_LFQ_Contaminant_table"),
            br(),
            br(),
            downloadButton("dannotated_mean_table_lfq", "Mean_Log2_LFQ_data"),
            br(),
            br(),
            downloadButton("dwide_sampletype_sort_table_lfq", "Mean_Log2_LFQ_bySampleTye_data"),
            br(),
            br(),
            downloadButton("dwide_type_sort_table_lfq", "Mean_Log2_LFQ_byType_data")
          )
        )
      )
      #
    )
  )
  ###################################################################################
  
)




