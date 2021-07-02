library(shiny)
library(rgl)
library(periscope)
library("magrittr")
library(plotly)
library(Cairo)

# library(rgl)
# library(crosstalk)
# library(manipulateWidget)
# library(webshot)



# # Define UI for application of Lifespan Study Data


fluidPage(
  
  tags$head(tags$style("
                  #container * {  
   display: inline;inline-block;
   vertical-align: top;
                     }")),
  div(id="container",img(src = "mitochondrial-signaling-lab-logo.png", height = 36, width = 84), h2("Cellular Lifespan Study")),
  h6("All data was generated at the Mitochondrial Psychobiology Lab, Columbia University Medical Center"),
  br(),
  
  # allow for unlimited time to keep webpage active
  tags$head(
    HTML(
      "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 15000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
    )
  ),
  textOutput("keepAlive"),
  
  
  sidebarLayout(
    
    sidebarPanel(
      
      helpText("Type/Select one or more Cell Lines:"),
      uiOutput("CellLineSelector"),
      
      helpText("Type/Select one or more Treatments:"),
      uiOutput("TreatmentSelector"),
      
      
      # checkboxGroupInput("replicates", "Replicate: ",
      #   choices = c(1,2,3),
      #   selected = 1
      # ),
      checkboxInput('legend', 'Show Legend', TRUE),
      
      helpText("Select one or more Experimental Replicates:"),
      checkboxInput('replicates_1', 'Replicate 1', TRUE),
      checkboxInput('replicates_2', 'Replicate 2 (Betahydroxybutyrate, 2-Deoxyglucose, Galactose)', FALSE),
      checkboxInput('replicates_3', 'Replicate 3 (3% Oxygen)', FALSE),
      
      helpText("Normalize to:"),
      checkboxInput('Norm_Cell_Volume', 'Cell Volume', FALSE),
      checkboxInput('Division_Norm', 'Divisions', FALSE),
      checkboxInput('Age_Norm', 'Age of Donor', FALSE),
      
      radioButtons("scale", "Select Y-axis scale:",
                   c("Linear" = "linear",
                     "Log" = "log")),
      radioButtons("outliers", "Remove Outliers:",
                   c("No" = "no",
                     "Yes" = "yes")),
      radioButtons("x_axis", "Select X axis:",
                   c("Days Grown" = "days_grown",
                     "MiAge Corrected Days Grown" = "miage_days_grown",
                     "Population Doublings" = "doublings_axis",
                     "MiAge Population Doublings" = "miage_doublings_axis",
                     "Culture Passages" = "passage_axis")),
      checkboxInput('second_axis', 'Overlay Growth Curve', FALSE),
      checkboxInput('fit', 'Show Fit Line', TRUE),
      checkboxInput('fit_text', 'Show Fit Text', TRUE),
      radioButtons("annotation", "Annotation Text:",
                   c("None" = "none",
                     "Cell Line Group" = "cell_line",
                     "All stats" = "all")),
      #checkboxInput('cell_line_text', 'Show Cell Line', TRUE),
      #checkboxInput('fit_text', 'Show Fit Line Equation', FALSE),
      checkboxInput('se', 'Show Confidence Intervals', FALSE),
      
      sliderInput("poly", "Degrees of Freedom in Polynomial Fit:",
                  min = 1, max = 10, step = 1, value = 3
      ),
      
      
      br(),
      helpText("Type/Select which dataset to download:"),
      uiOutput("DatasetSelector"),
      downloadButton("downloadData", "Download selected Data"),
      #br(), br(),
      #actionButton("saveGraphs", "Save Plots")
      #downloadButton("downloadGraphs", "Download all Graphs")
      
      # radioButtons("fileType", "Download Plot as:",
      #              c("png" = "png",
      #                "pdf" = "pdf",
      #                "jpeg" = "jpeg",
      #                "tiff" = "tiff"))
      
      # # Input: Simple integer interval ----
      # helpText("Select Time Grown for Energy Phenotypes"),
      # sliderInput("days_grown", "Days Grown:",
      #             min = 0, max = 300,
      #             value = 300),
      
      # animationOptions(interval = 10, loop = FALSE, playButton = TRUE,
      #                  pauseButton = TRUE),
      

    ),
    
    mainPanel(
      tabsetPanel(id = "tabs", 
        #tabPanel("Seahorse Data", plotOutput("Seahorse_Plots", click = "time_plot_click"), verbatimTextOutput("DataPoint")),
        tabPanel("Growth Curves", id = "Growth_Curves", 
                 plotlyOutput("Growth_Plots"),
                 br(),
                 checkboxInput('adj_population_doublings', 'Correct population doublings for pre-study divisions based off MiAge Clock', FALSE),
                 # column(11, align="right",
                 # downloadButton("downloadPlot_Growth_Curve","Download Plot")),
                 # br(),
                 # br(),
                 # verbatimTextOutput("Growth_DataPoint"),
                 br(),
                 plotlyOutput("Growth_Rate_Plots"),
                 br(),
                 plotlyOutput("Cell_Size_Plots"),
                 checkboxInput('diameter_conversion', 'Convert to diameter of cells', FALSE),
                 # column(11, align="right",
                 #        downloadButton("downloadPlot_Cell_Size","Download Plot")),
                 # br(),
                 # br(),
                 # verbatimTextOutput("Size_DataPoint"),
                 br(),
                 plotlyOutput("Percent_Dead_Plots")
                 # column(11, align="right",
                 #       downloadButton("downloadPlot_Percent_Dead","Download Plot")),
                 # br(),
                 # br(),
                 # verbatimTextOutput("Dead_DataPoint")
        ),
        tabPanel("DNA Methylation Age", id = "DNAmAge", 
                 helpText("Type/Select a DNAmAge Clock:"),
                 uiOutput("ClockSelector"),
                 br(), br(),
                 plotlyOutput("DNAmAge_Plots"), 
                 br()
                 #  column(11, align="right",
                 #      downloadButton("downloadPlot_DNAm","Download Plot")),
                 #      
                 #             #plotOutput("DNAmAge_BarPlot") 
                 # 
                 # br(),
                 # br(),
                 # verbatimTextOutput("DNAmAge_DataPoint")
        ),
        # tabPanel("DNA Methylation CpGs", id = "CpG",
        #          helpText("Type/Select a CpG:"),
        #          #uiOutput("CpGSelector"),
        #          selectInput("CpG", "CpG", multiple = FALSE, choices = character(0)),
        #          br(), br(),
        #          plotlyOutput("CpG_Plot"),
        #          br(), br(),
        #          actionButton("PCA_button", "Run Principal Component Analysis"),
        #          br(), br(),
        #          plotlyOutput("PCA_Plot"),
        #          br()
        # 
        #          #verbatimTextOutput("Correlation_DataPoint"),
        #          #br(), br(), br(), br(), br(), br(), br()
        #          #helpText("3D PCA"),
        #          #rglwidgetOutput("PCA_3D_Plot", width = 800, height = 600)
        # ),
        
        tabPanel("RNAseq Genes", id = "RNAseq",
                 helpText("Type/Select a Gene:"),
                 selectInput("Gene", "Gene", multiple = FALSE, choices = character(0)),
                 checkboxInput('log2_genes', 'Show Log2 expression values', F),
                 br(), br(),
                 plotlyOutput("RNA_Gene_Plot"),
                 br(), br(),
                 #helpText("Type/Select a Transcript:"),
                 #selectInput("Transcript", "Transcript", multiple = FALSE, choices = character(0)),
                 #br(), br(),
                 #plotlyOutput("RNA_Transcript_Plot"),
                 br(), br(),
                 
                 actionButton("PCA_button_2D_RNA", "Run Principal Component Analysis 2D"),
                 br(), br(),
                 plotOutput("PCA_Plot_2D_RNA"),
                 br(),
                 actionButton("PCA_button_3D_RNA", "Run Principal Component Analysis 3D"),
                 br(),
                 plotlyOutput("PCA_Plot_3D_RNA"),
                 br()
                 
                 #verbatimTextOutput("Correlation_DataPoint"),
                 #br(), br(), br(), br(), br(), br(), br()
                 #helpText("3D PCA"),
                 #rglwidgetOutput("PCA_3D_Plot", width = 800, height = 600)
        ),
        
        tabPanel("Telomere Length", id = "TL", 
                 plotlyOutput("Telomere_Plots"), 
                 br()
                 # column(11, align="right",
                 #        downloadButton("downloadPlot_TL","Download Plot")),
                 # br(),
                 # br(),
                 # verbatimTextOutput("TL_DataPoint")
        ),
        tabPanel("Copy Number", id = "copy_number", 
                 plotlyOutput("CN_Plots"), 
                 br()
                 # column(11, align="right",
                 #        downloadButton("downloadPlot_Copy_Number","Download Plot")),
                 # br(),
                 # br(),
                 # verbatimTextOutput("CN_DataPoint")
        ),
        tabPanel("mtDNA Deletions", id = "Deletions", 
                 helpText("Select mtDNA Deletion Parameter:"),
                 uiOutput("Deletions_ParameterSelector"),
                 plotlyOutput("Deletions_Plots"),
                 br(), br()
                 # checkboxInput('Cytokine_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 # br()
        ),
        tabPanel("cell free DNA", id = "cf-DNA", 
                 radioButtons("cf-DNA", "Select DNA:",
                              c("mitochondrial DNA" = "mtDNA",
                                "nuclear DNA" = "nDNA",
                                "mtDNA/nDNA" = "both",
                                "mtDNA/copy number" = "copy_number")),
                 plotlyOutput("cf_Plots"), 
                 br(), br(),
                 checkboxInput('cf_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 br()
                 # column(11, align="right",
                 #        downloadButton("downloadPlot_cfDNA","Download Plot")),
                 # br(),
                 # br(),
                 # verbatimTextOutput("cfDNA_DataPoint")
        ),
        tabPanel("cell free IL-6", id = "IL6", 
                 plotlyOutput("IL6_Plots"), 
                 br(), br(),
                 checkboxInput('IL6_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 br()
                 # column(11, align="right",
                 #        downloadButton("downloadPlot_IL6","Download Plot")),
                 # br(),
                 # br(),
                 # verbatimTextOutput("IL6_DataPoint")
        ),
        tabPanel("Cytokines", id = "cytokines", 
                 helpText("Select Cytokine:"),
                 uiOutput("Cytokines_ParameterSelector"),
                 plotlyOutput("Cytokine_Plots"),
                 br(), br(),
                 checkboxInput('Cytokine_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 br()
        ),
        tabPanel("Seahorse Bioenergetics", id = "Seahorse", 
                helpText("Type/Select one or more Bioenergetic Parameters:"),
                uiOutput("ParameterSelector"),
                plotlyOutput("Seahorse_Plots"),
                br(), br(),
                checkboxInput('Run_norm', 'Normalize out Plate-to-Plate Variation', FALSE),
                br()
                # column(11, align="right",
                #        downloadButton("downloadPlot_Seahorse","Download Plot")),
                # br(),
                # br(),
                # verbatimTextOutput("Seahorse_DataPoint")
        ),
        tabPanel("Correlations", id = "correlation", 
                 helpText("2D Correlation Plot"),
                 column(3, 
                        helpText("Type/Select X axis parameter:"),
                        uiOutput("ParameterSelector2D_x")
                  ),
                 column(3, 
                        helpText("Type/Select Y axis parameter:"),
                        uiOutput("ParameterSelector2D_y")
                  ),
                 br(),br(), br(), br(), br(), 
                 plotOutput("Correlation_2D_Plots"),
                 br(), br(), br(), br(), br(),
                 checkboxInput('Fun_Regression', 'Perform Functional Regression', TRUE),
                 #checkboxInput('Optomize_Fit', 'Optomize the Degrees of Freedom for best Fit', FALSE),
                 sliderInput("poly_x", "Degrees of Freedom in Polynomial Fit for Parameter X:",
                             min = 1, max = 10, step = 1, value = 3
                 ),
                 sliderInput("poly_y", "Degrees of Freedom in Polynomial Fit for Parameter Y:",
                             min = 1, max = 10, step = 1, value = 3
                 ),
                 br(),br(),
                 #verbatimTextOutput("Correlation_DataPoint"),
                 br(), br(), br(), br(), br(), br(), br(),
                 helpText("3D Correlation Plot"),
                 column(3, 
                        helpText("Type/Select X axis parameter:"),
                        uiOutput("ParameterSelector3D_x")),
                 column(3, 
                        helpText("Type/Select Y axis parameter:"),
                        uiOutput("ParameterSelector3D_y")),
                 column(3, 
                       helpText("Type/Select Z axis parameter:"),
                       uiOutput("ParameterSelector3D_z")),
                 checkboxInput('ellipsoid', 'Show Ellipoids', TRUE),
                 rglwidgetOutput("Correlation_3D_Plots", width = 800, height = 600)
        ),
       
        tabPanel("Data Table", 
                 br(),
                 helpText("**Use sidebar to select dataset of interest**"),
                 br(),
                 tableOutput("DataMatrix")
        ),
         tabPanel("About", 
                  br(),
                  h4("Shiny App Info"),
                  div(
                    p("This app was built using the R (v4.0.2) library 'shiny' (v1.6.0) in Rstudio (v1.3.10) . These data were generated as part of the Cellular Lifespan study which
                    contains genome-wide, multimodal data collected at high temporal frequency across the replicative lifespan 
                    in cultured primary human fibroblasts. \n
                       For a more detailed explanation of the study design see: Sturm et al., Scientific Data 2021 "),
                    p("Primary findings related to mitochondrial disease were published here: Sturm et al., Nature Metabolism 2021"),
                    p("Primary findings related to metabolic aging were published here: Sturm et al., Nature Aging 2021"),
                    #br(),
                    p("The full RNAseq dataset can be found here: GEO #XXXXX"),
                    p("The full DNAmethylation dataset can be found here: GEO #XXXXX"),
                    #br(),
                    hr(),
                    h6("This app was built by Gabriel (Gav) Sturm at the Mitochondrial Pyschobiology Lab, Columbia University Medical Center."),
                    h6("The app's source code is available at: https://github.com/gav-sturm/Cellular_Lifespan_Study"),
                    h6("All app-related questions can be pointed towards; gs2934@cumc.columbia.edu"),
                    br(),
                    h6("We are additionally looking to expand this database to include additional cellular lifespan projects."),
                    h6("For potential collaborations please contact; martin.picard@columbia.edu"),
                    #br(),
                    #h6("To learn more about the lab visit:"), 
                    hr(),
                    uiOutput("labLink",style = "font-size:14px")
                    #br()
                  )
        )
      )
    )
  )
)