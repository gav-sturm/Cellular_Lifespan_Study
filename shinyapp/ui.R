#####################################
##
## Title: Cellular Lifespan Study ShinyApp UI code
## Author: Gabriel (Gav) Sturm
## Contact Info: gs2934@cumc.columbia.edu
## Date: 2021-01-01
## 
##
##

library(shiny)
library(periscope)
library("magrittr")
library(plotly)
library(Cairo)
library(shinyBS)

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
  h6("Mitochondrial Psychobiology Lab, Columbia University Medical Center"),
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
      
      checkboxInput('legend', 'Show Legend', TRUE),
      
      helpText("Select one or more phases of the study (for matching ctrls):"),
      checkboxInput('replicates', 'Phase I (mitoQ,NAC,a-keto) ', TRUE),
      checkboxInput('replicates_1', 'Phase II (DEX,mitoNUIts,Oligo) ', TRUE),
      checkboxInput('replicates_2', 'Phase III (contact inhibition,hypoxia)', FALSE),
      checkboxInput('replicates_3', 'Phase IV (galactose,hydroxybuturate,2DG)', FALSE),
      checkboxInput('replicates_4', 'Phase V (SURF1 hypoxia)', FALSE),
      
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
                     "Days after treatment" = "days_treatment",
                     #"MiAge Corrected Days Grown" = "miage_days_grown",
                     "Population Doublings" = "doublings_axis",
                     #"MiAge Population Doublings" = "miage_doublings_axis",
                     "Culture Passages" = "passage_axis")),
      sliderInput("time_slider", "Percent lifespan (x-axis) shown:",
                  min = 0, max = 100, step = 5, value = c(0,100)
      ),
      
      checkboxInput('second_axis', 'Overlay Growth Curve', FALSE),
      checkboxInput('fit', 'Show Fit Line', TRUE),
      checkboxInput('fit_text', 'Show Fit Text', TRUE),
      radioButtons("annotation", "Annotation Text:",
                   c("None" = "none",
                     "Cell Line Group" = "cell_line",
                     "All stats" = "all")),
      #checkboxInput('cell_line_text', 'Show Cell Line', TRUE),
      #checkboxInput('fit_text', 'Show Fit Line Equation', FALSE),
      checkboxInput('se', 'Show 95% Confidence Intervals', FALSE),
      
      sliderInput("poly", "Degrees of Freedom in Polynomial Fit:",
                  min = 1, max = 10, step = 1, value = 3
      ),
      
      
      br(),

    ),
    
    mainPanel(
      tabsetPanel(id = "tabs", 
        #tabPanel("Seahorse Data", plotOutput("Seahorse_Plots", click = "time_plot_click"), verbatimTextOutput("DataPoint")),
        tabPanel("Growth Curves", id = "Growth_Curves", icon = icon("chart-bar"),
                 h4("Population doublings"),
                 #splitLayout(cellWidths = c("75%", "25%"),plotlyOutput("Growth_Plots"), plotlyOutput("Growth_BarPlot")),
                 plotlyOutput("Growth_Plots"),
                 checkboxInput('adj_population_doublings', 'Correct population doublings for pre-study divisions based off MiAge Clock', FALSE),
                 #actionButton("btn", "see method details"),
                 tipify(el = bsButton("methods_btn_PopDoublings", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Cells were passaged approximately every 6 days (+/- 1 day), with decreasing passaging frequency as cells enter quiescence, for up to 270 days. Population doublings were determined by cell counting at each passage using a Countess II Automated Cell Counter (Thermofisher).", 
                        placement="right", trigger = "hover"), 
                 br(),
                 h4("Growth rates"),
                 plotlyOutput("Growth_Rate_Plots"),
                 #splitLayout(cellWidths = c("75%", "25%"),plotlyOutput("Growth_Rate_Plots"), plotlyOutput("Growth_Rate_BarPlot")),
                 tipify(el = bsButton("methods_btn_GR", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Growth rates were determined at each passage using a Countess II Automated Cell Counter (Thermofisher).", 
                        placement="right", trigger = "hover"), 
                 br(),
                 h4("Cell size"),
                 #splitLayout(cellWidths = c("75%", "25%"),plotlyOutput("Cell_Size_Plots"), plotlyOutput("Cell_Size_BarPlot")),
                 plotlyOutput("Cell_Size_Plots"),
                 checkboxInput('diameter_conversion', 'Convert to diameter of cells', FALSE),
                 tipify(el = bsButton("methods_btn_Size", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Cell volume was determined at each passage using a Countess II Automated Cell Counter (Thermofisher).", 
                        placement="right", trigger = "hover"), 
                 br(),
                 #br(),
                 h4("Dead cells"),
                 plotlyOutput("Percent_Dead_Plots"),
                 #splitLayout(cellWidths = c("75%", "25%"),plotlyOutput("Percent_Dead_Plots"), plotlyOutput("Percent_Dead_BarPlot")),
                 tipify(el = bsButton("methods_btn_Dead", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "The proporation of dead cells at each passage was determined using trypan blue staining under a Countess II Automated Cell Counter (Thermofisher).", 
                        placement="right", trigger = "hover"), 
                 br(),
                 br(), br()
        ),
        tabPanel("DNAmethylation", id = "DNAmethylation", icon = icon("chart-bar"),
                 h4("DNAmethylation Age"),
                 helpText("Type/Select a Epigenetic Clock:"),
                 uiOutput("ClockSelector"),
                 plotlyOutput("DNAmAge_Plots"),
                 tipify(el = bsButton("methods_btn_Clocks", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Epigenetic clocks were computed for four clocks of chronological age: Horvath1 (i.e. PanTissue clock), Horvath2 (i.e. Skin&Blood clock), Hannum, and PedBE clocks; two clocks designed to predict mortality, the PhenoAge and GrimAge clocks; a clock to measure telomere length, DNAmTL; a clock designed to measure mitotic age, MiAge; a clock trained to predict cellular senescence, DNAmSen, and two DNA methylation measure of the rate of deterioration in physiological integrity, DunedinPoAm, and DundedinPACE.",
                        placement="right", trigger = "hover"), 
                 br(),
                 
                 h4("DNAmethylation CpGs"),
                 helpText("Type/Select a CpG: (Note: only age-related CpGs shown, ~40k. See 'About' page for full dataset)"),
                 selectInput("CpG", "CpG", multiple = FALSE, choices = character(0)),
                 plotlyOutput("CpG_Plot"),
                 tipify(el = bsButton("methods_btn_DNAm", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Global DNA methylation was measured using Illumina EPIC microarray after bisulfite conversion and hybridization using the Infinium Methylation EPIC BeadChip kit. Values are expressed as normalized beta values after RCP and Combat adjustments.",
                        placement="right", trigger = "hover"), 
                 br(),
                 actionButton("PCA_button", "Run Principal Component Analysis 2D"),
                 br(),
                 plotlyOutput("PCA_Plot"),
                 br(),
                 actionButton("PCA_button_3D_DNAm", "Run Principal Component Analysis 3D"),
                 helpText("Hover and click to explore 3D space:"),
                 plotlyOutput("PCA_Plot_3D_DNAm"),
                 br(),br(),br(),br(),br(),br(),br(),br()
        ),
        tabPanel("Gene expression", id = "RNAseq",icon = icon("chart-bar"),
                 h4("Gene Expression (RNAseq)"),
                 helpText("Type/Select a Gene:"),
                 selectInput("Gene", "Gene", multiple = FALSE, choices = character(0)),
                 checkboxInput('log2_genes', 'Show Log2 expression values', F),
                 #br(), br(),
                 plotlyOutput("RNA_Gene_Plot"),
                 tipify(el = bsButton("methods_btn_RNAseq", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Gene expression was measured using RNA sequencing on RbiZero gold extracted mRNA converted to cDNA with NEBNext® Ultra™ II RNA Library Prep Kit and run on a Illumina HiSeq 4000 with single index paired-end 150bp chemistry.",
                        placement="right", trigger = "hover"), 
                 br(),
                 actionButton("PCA_button_2D_RNA", "Run Principal Component Analysis 2D"),
                 #br(), br(),
                 plotOutput("PCA_Plot_2D_RNA"),
                 actionButton("PCA_button_3D_RNA", "Run Principal Component Analysis 3D"),
                 helpText("Hover and click to explore 3D space:"),
                 plotlyOutput("PCA_Plot_3D_RNA"),
                 br(),br(),br(),br(),br(),br(),br(),br()
        ),
        
        tabPanel("Telomere Length", id = "TL", icon = icon("chart-bar"),
                 h4("Telomere Length"),
                 plotlyOutput("Telomere_Plots"), 
                 br(),
                 tipify(el = bsButton("methods_btn_TL", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Relative telomere length was measured by qPCR and expressed as the ratio of telomere to the single-copy gene, human beta-globin, (T/S ratio), as previously described Cawthon et al. Nucleic Acids Res. 2002 and Lin et al. J. Immunol. Methods, 2010.",
                        placement="right", trigger = "hover"), 
                 br(), br()
        ),
        tabPanel("mtDNA", id = "mtDNA", icon = icon("chart-bar"),
                 h4("mtDNA Copy Number"),
                 plotlyOutput("CN_Plots"), 
                 tipify(el = bsButton("methods_btn_cn", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Cellular mtDNA content was quantified by qPCR on the same genomic material used for other DNA-based measurements. Duplex qPCR reactions with Taqman chemistry were used to simultaneously quantify mitochondrial (mtDNA, ND1) and nuclear (nDNA, B2M) amplicons, as described previously Picard et al. PNAS 2014.",
                        placement="right", trigger = "hover"), 
                 br(),br(),
                 h4("mtDNA sequencing"),
                 helpText("Select mtDNA Deletion Parameter:"),
                 uiOutput("Deletion_ParameterSelector"),
                 plotlyOutput("Deletions_Plots"),
                 tipify(el = bsButton("methods_btn_deletions", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "mtDNA next-generation sequencing was perfomed using an Ion Torrent S5XL platform followed by eKLIPse analysis as described in Goudenège et al. Nature Genetics 2018",
                        placement="right", trigger = "hover"), 
                 br(), br()
        ),
        tabPanel("cell-free DNA", id = "cf-DNA", icon = icon("chart-bar"),
                 h4("cell-free DNA"),
                 radioButtons("cf-DNA", "Select DNA:",
                              c("mitochondrial DNA" = "mtDNA",
                                "nuclear DNA" = "nDNA",
                                "mtDNA/nDNA" = "both",
                                "mtDNA/copy number" = "copy_number")),
                 plotlyOutput("cf_Plots"), 
                 checkboxInput('cf_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 tipify(el = bsButton("methods_btn_cfDNA", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Cell free mtDNA and nDNA levels were measured simultaneously by taqman-based duplex qPCR reactions targeted mitochondrial-encoded ND1 and nuclear-encoded B2M sequences as described by Ware et al., Journal of Biological Chemistry, 2020.",
                        placement="right", trigger = "hover"), 
                 br(), br()
        ),
        tabPanel("Cytokines", id = "IL6", icon = icon("chart-bar"),
                 h4("Cytokine Array"),
                 helpText("Select Cytokine:"),
                 uiOutput("Cytokines_ParameterSelector"),
                 plotlyOutput("Cytokine_Plots"),
                 checkboxInput('Cytokine_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 tipify(el = bsButton("methods_btn_Cytokines", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Cell-free cytokines were quantified on a Luminex 200 using two multiplex fluorescence-based arrays with custom-designed with selected cytokines and chemokines based on human age-related proteins. Analytes were selected based on reported correlations of their levels in human plasma with chronological age idenified in Tanaka et al., Aging Cell, 2018.",
                        placement="right", trigger = "hover"),
                 br(),
                 h4("cell-free IL6 ELISA"),
                 #helpText("cell-free IL6 ELISA:"),
                  plotlyOutput("IL6_Plots"), 
                 checkboxInput('IL6_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 tipify(el = bsButton("methods_btn_IL6", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Cell-free IL6 was quantified using enzyme-linked immunosorbent assay (ELISA,  Abcam #ab229434).",
                        placement="right", trigger = "hover"),
                 br(),
                 h4("cell-free GDF15 ELISA"),
                 #helpText("cell-free GDF15 ELISA:"),
                 plotlyOutput("GDF15_Plots"), 
                 checkboxInput('GDF15_norm', 'Normalize to Cells Counted at Passage', TRUE),
                 tipify(el = bsButton("methods_btn_GDF15", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                        "Cell-free GDF15 was quantified using enzyme-linked immunosorbent assay (ELISA, R&D #DGD150).",
                        placement="right", trigger = "hover"),
                 br(),
                 br(), br()
        ),
        tabPanel("Seahorse Bioenergetics", id = "Seahorse", icon = icon("chart-bar"),
                h4("Seahorse Bioenergetics"),
                helpText("Type/Select one or more Bioenergetic Parameters:"),
                uiOutput("ParameterSelector"),
                plotlyOutput("Seahorse_Plots"),
                checkboxInput('Run_norm', 'Normalize out Plate-to-Plate Variation', FALSE),
                tipify(el = bsButton("methods_btn_Seahorse", "hover to see method details", icon = icon("info-circle"), style = "inverse", size = "small"), 
                       "Bioenergetic parameters were measured using the XFe96 Seahorse extracellular flux analyzer (Agilent) using the MitoStress Test. ATP production rates (JATP) were estimated as described by Mookerjee et al., Journal of Biological Chemistry, 2017",
                       placement="right", trigger = "hover"),
                br(),
                br(),br()
        ),
        tabPanel("Correlations", id = "correlation", icon = icon("chart-bar"),
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
                 helpText("Hover and click to explore 3D space:"),
                 plotlyOutput("Correlation_3D_Plots"),
                 br(),br(),br(),br(),br(),br(),br(),br()
        ),
       
        tabPanel("Download Data", icon = icon("table"),
                 br(),
                  helpText("Type/Select which parameter to download:"),
                 uiOutput("DatasetSelector"),
                
                 downloadButton("downloadData", "Download selected Data"),
                 br(),
                 helpText("**Use control panel on left to select cell line(s) & treatments(s) of interest, then use the above dropdown menu to select the parameter of interest.**"),
                 checkboxInput('full_download', 'Click to download the full sample-set (~2000k timepoints)', F),
                 h5("Raw Dataset Links:"),
                 uiOutput("geoSuperSeriesLink",style = "font-size:14px"),
                 uiOutput("geoDNAmLink2",style = "font-size:14px"),
                 uiOutput("geoRNAseqLink2",style = "font-size:14px"),
                 uiOutput("figshareImageLink2",style = "font-size:14px"),
                 uiOutput("figshareDataLink2",style = "font-size:14px"),
                 br(),br(),
                 tableOutput("DataMatrix"),
                 br(), br()
        ),
         tabPanel("About", icon = icon("list-alt"),
                  br(),
                  h4("Shiny App Info"),
                  div(
                    p("This app was built using the R (v4.0.2) library 'shiny' (v1.6.0) in Rstudio (v1.3.10) . These data were generated as part of the Cellular Lifespan study which
                    contains genome-wide, multimodal data collected at high temporal frequency across the replicative lifespan 
                    in cultured primary human fibroblasts. \n
                        "),
                    h5("Manuscript Links:"),
                    uiOutput("scidataLink",style = "font-size:14px"), # 
                    uiOutput("surf1Link",style = "font-size:14px"), # https://www.biorxiv.org/content/10.1101/2022.02.22.481548v1
                    uiOutput("dexLink",style = "font-size:14px"), #
                    p("Primary findings related to metabolic aging were published here: Sturm et al., XXXXXX 2021"),
                    h6("See 'Download Data' page for portal to raw and processed datasets."),
                    #br(),
                    #h5("Raw Dataset Links:"),
                    #uiOutput("geoSuperSeriesLink",style = "font-size:14px"),
                    #uiOutput("geoDNAmLink",style = "font-size:14px"),
                    #uiOutput("geoRNAseqLink",style = "font-size:14px"),
                    #uiOutput("figshareImageLink",style = "font-size:14px"),
                    #uiOutput("figshareDataLink",style = "font-size:14px"),
                    #p("The full DNA methylation dataset can be found here: GSE179847"),
                    #p("The full RNAseq dataset can be found here: GSE179848"),
                    #br(),
                    hr(),
                    h6("This app was built by Gabriel (Gav) Sturm at the Mitochondrial Pyschobiology Lab, Columbia University Medical Center."),
                    uiOutput("githubLink",style = "font-size:14px"),
                    #h6("The app's source code is available at: https://github.com/gav-sturm/Cellular_Lifespan_Study"),
                    h6("All app-related questions can be pointed towards: gs2934@cumc.columbia.edu"),
                    br(),
                    h6("We are additionally looking to expand this database to include additional cellular lifespan projects."),
                    h6("For potential collaborations please contact: martin.picard@columbia.edu"),
                    #br(),
                    #h6("To learn more about the lab visit:"), 
                    hr(),
                    uiOutput("labLink",style = "font-size:14px"),
                    br(), br()
                  )
        )
      )
    )
  )
)