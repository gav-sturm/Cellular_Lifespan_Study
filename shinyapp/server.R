#####################################
##
## Title: Cellular Lifespan Study ShinyApp Server code
## Author: Gabriel (Gav) Sturm 
## Date: 2021-01-01
## 
##
##


library(shiny)
library(ggplot2)
library(ggpmisc)
library(rgl)
library(car)
library(gridExtra)
library(plotly)
library(Cairo)
library(polynom)
library(stringr)
library(mgcv)
library(processx)
library(orca)

## Load Dataset
# setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp")
# dir()
LS_Data <- read.csv("Lifespan_Study_Data.csv")


# Parameter Lists
Seahorse_Parameters <- c("ATPtotal","ATPglyc", "ATPox", "ATPtotal_max", "ATPglyc_max","ATPox_max"	, "ATPtotal_spare",	"ATPox_spare","ATPglyc_spare","Resting_Metabolic_Rate","Max_Metabolic_Rate", 
                         "PicoWatts_per_cell",	"Max_PicoWatts_per_cell",	"Watts_per_kilogram",	"Max_Watts_per_kilogram",
                          "Basal_Respiration", "Max_Respiration", "Oxidative_Spare_Capacity", "Proton_Leak", "Non_Mitochondrial_Respiration", 
                        "Baseline_ECAR", "Max_ECAR", "Glycolytic_Spare_Capacity", "BHI", "OCRcoupled", "PPRtotal", "PPRrespiration", 
                        "PPRglycolysis", "PPRtotal_max",	"PPRgly_max",	"PPRresp_max",	"ATP_Linked_Respiration", "Coupling_Efficiency", "Percent_Spare_Capacity")

Clock_Parameters <- c("Horvath1",	"Horvath2",	"Hannum",	"PhenoAge",	"DNAmTL",	"GrimAge",	"DNAmADM",	"DNAmB2M",	
  "DNAmCystatinC",	"DNAmGDF15",	"DNAmLeptin",	"DNAmPAI1",	"DNAmTIMP1",	"DNAmPACKYRS",	
  "PCHorvath1",	"PCHorvath2",	"PCHannum",	"PCPhenoAge",	"PCDNAmTL",	"PCGrimAge",	
  "Mitotic_Age",	"DunedinPoAm_Age",	"Mol_Skin_Clock",	"DNAmSen")
scol <- which(colnames(LS_Data)=="Lumican_BR13_18_plex" )
ncol <- which(colnames(LS_Data)=="IL_18_IL_1F4_BR78_35_plex" )
Cytokine_Parameters <- as.character(colnames(LS_Data)[scol:ncol])
Deletions_Parameters <- c("mtDNA_deletion_Haplogroupe",	"mtDNA_deletion_Nb_de_délétions",	
          "mtDNA_deletion_Mean_lenght",	"mtDNA_deletion_Min_lenght",	"mtDNA_deletion_Max_lenght",	"mtDNA_deletion_Mean_HT",	
          "mtDNA_deletion_Min_HT",	"mtDNA_deletion_Max_HT",	"mtDNA_deletion_Mean_Repeat_lenght",	"mtDNA_deletion_Min_Reapeat_lenght",
          	"mtDNA_deletion_Max_repeat_Lenght")


# calculator for the x axis limits and steps
select_x_max <- function (data) {
  x_max <- 300
  step <- 50
  if('Contact_Inhibition_Regrowth' %in% data$Treatments) {
    x_max <- 225
    step <- 25
  }
  else if(is.element('4', data$Study_Part)) {
    x_max <- 175
    step <- 25
  }
  
  else if('SURF1_Mutation' %in% data$Clinical_Condition) {
    x_max <- 150
    step <- 25
  }
  
  else if('hFB1' %in% data$Cell_Line | 'hFB2' %in% data$Cell_Line) {
    x_max <- 150
    step <- 25
  }
  
  else if('Contact_Inhibition' %in% data$Treatments) {
    x_max <- 125
    step <- 25
  }
  else if('hFB11' %in% data$Cell_Line | '3' %in% data$Percent_Oxygen | '3' %in% data$Replicate_Line) {
    x_max <- 80
    step <- 20
  }
  
  if('HEK293' %in% data$Cell_Line) {
    x_max <- 600
    step <- 100
  }

  else if (is.element('2', data$Study_Part)) {
    if('hFB12' %in% data$Cell_Line | 'hFB13' %in% data$Cell_Line | 'hFB14' %in% data$Cell_Line) {
      x_max <- 300
      step <- 50
    }
    
  }
  
  x_data <- c(x_max, step)
  return(x_data)
}

# pick treatment colors
getColors <- function(treatments) {
  colors <- c()
  for(i in 1:length(treatments)) {
    treatment <- treatments[i]
    if(!is.vector(treatments)) {
      treatment <- treatments
    }
    if(treatment =="Control_21") {
      colors <- append(colors, "dimgray")
    }
    else if(treatment =="Control_3") {
      colors <- append(colors, "cadetblue4")
    }
    else if(treatment =="DEX_21") {
      colors <- append(colors, "red")
    }
    else if(treatment =="DEX_3") {
      colors <- append(colors, "chocolate4")
    }
    else if(treatment =="Pulsated_DEX_21") {
      colors <- append(colors, "orange")
    }
    else if(treatment =="MitoQ_21") {
      colors <- append(colors, "lightblue")
    }
    else if(treatment =="MitoQ+DEX_21") {
      colors <- append(colors, "blue")
    }
    else if(treatment =="NAC_21") {
      colors <- append(colors, "lightgreen")
    }
    else if(treatment =="NAC+DEX_21") {
      colors <- append(colors, "green")
    }
    else if(treatment =="a-ketogluturate_21") {
      colors <- append(colors, "yellow")
    }
    else if(treatment =="Modulators_21") {
      colors <- append(colors, "green")
    }
    else if(treatment =="Modulators+DEX_21") {
      colors <- append(colors, "darkgreen")
    }
    else if(treatment =="Oligomycin_21") {
      colors <- append(colors, "blue")
    }
    else if(treatment =="Oligomycin+DEX_21") {
      colors <- append(colors, "darkblue")
    }
    else if(treatment =="5-azacytidine_21") {
      colors <- append(colors, "pink")
    }
    else if(treatment =="5-azacytidine+Modulators_21") {
      colors <- append(colors, "purple")
    }
    else if(treatment =="Contact_Inhibition_21") {
      colors <- append(colors, "cyan2")
    }
    else if(treatment =="Contact_Inhibition_Regrowth_21") {
      colors <- append(colors, "cyan3")
    }
    else if(treatment =="Contact_Inhibition_3") {
      colors <- append(colors, "darkcyan")
    }
    else if(treatment =="betahydroxybutyrate_21") {
      colors <- append(colors, "darkorchid")
    }
    else if(treatment =="2-Deoxyglucose_21") {
      colors <- append(colors, "gold1")
    }
    else if(treatment =="Galactose_21") {
      colors <- append(colors, "darkorange1")
    }
  }
  return(colors)
}

# pick cell line colors
getColors_2 <- function(cell_lines) {
  colors <- c()
  for(i in 1:length(cell_lines)) {
    cell_line <- cell_lines[i]
    if (cell_line =="hFB12" | cell_line =="hFB13" | cell_line =="hFB14" | cell_line =="hFB1" | cell_line =="hFB2" | cell_line =="hFB11") {
      colors <- append(colors, "lightgray")
    }
    else if (cell_line =="hFB6" | cell_line =="hFB7" | cell_line =="hFB8") {
      colors <- append(colors, "blue")
    }
    else if (cell_line == "HEK293") {
      colors <- append(colors, "red")
    }
  }
  return(colors)
}

# pick shapes of cell line symbols
getShapes <- function(cell_lines) {
  shapes <- c()
  for(i in 1:length(cell_lines)) {
    cell_line <- cell_lines[i]
    if(cell_line =="hFB12") {
      shapes <- append(shapes, 21)
    }
    else if(cell_line =="hFB13") {
      shapes <- append(shapes, 22)
    }
    else if(cell_line =="hFB14") {
      shapes <- append(shapes, 23)
    }
    else if(cell_line =="hFB11") {
      shapes <- append(shapes, 24)
    }
    else if(cell_line =="hFB1") {
      shapes <- append(shapes, 21)
    }
    else if(cell_line =="hFB2") {
      shapes <- append(shapes, 22)
    }
    else if(cell_line =="hFB6") {
      shapes <- append(shapes, 21)
    }
    else if(cell_line =="hFB7") {
      shapes <- append(shapes, 22)
    }
    else if(cell_line =="hFB8") {
      shapes <- append(shapes, 23)
    }
    else if(cell_line =="HEK293") {
      shapes <- append(shapes, 21)
    }
  }
  return(shapes)
}


##### Define server logic #####
shinyServer(function(input, output, session) {
  
  # hyperlink
  url <- a("http://www.picardlab.org/", href="http://www.picardlab.org/")
  output$labLink <- renderUI({
    tagList("Learn more about our lab here:", url)
  })
  
  #Cell Line list
  output$CellLineSelector<-renderUI({
    selectInput('Cell_Line', 'Cell Lines',
                c("HC1 (hFB12): primary human fibroblast, breast, healthy, male, 18yo" = "hFB12",
                  "HC2 (hFB13): primary human fibroblast, breast, healthy, female, 18yo" = "hFB13",
                  "HC3 (hFB14): primary human fibroblast, foreskin, healthy, male, 0yo" = "hFB14",
                  "HC4 (hFB11): primary human fibroblast, breast, healthy, male, 18yo" = "hFB11",
                  "HC5 (hFB1): primary human fibroblast, upper arm, healthy, male, 29yo"= "hFB1",
                  "HC6 (hFB2): primary human fibroblast, upper arm, healthy, female, 36yo" = "hFB2",
                  "SURF1_1 (hFB6): primary human fibroblast, upper arm, SURF1 Mutation, male, 1yo" = "hFB6",
                  "SURF1_2 (hFB7): primary human fibroblast, upper arm, SURF1 Mutation, male, 11yo" = "hFB7",
                  "SURF1_3 (hFB8): primary human fibroblast, upper arm, SURF1 Mutation, female, 9yo" = "hFB8",
                  "HEK293 Immortalized, female" = "HEK293"), 
                multiple=TRUE, 
                selectize=TRUE, 
                selected="hFB12") #default value
  })
  
  
  #Treatment list
  output$TreatmentSelector<-renderUI({
    selectInput('Treatment', 'Treatments',
                c("Untreated Control (HC1,2,3,4,5,6, SURF1_1,2,3)" = "Control_21",
                  "3% Oxygen (HC1,2,4)" = "Control_3",
                  "Contact Inhibition (HC1,2,4,5)" = "Contact_Inhibition_21",
                  "3% Oxygen + Contact_Inhibition_21 (HC1,2,4)" = "Contact_Inhibition_3",
                  "Contact Inhibition Regrowth (HC5)" = "Contact_Inhibition_Regrowth_21",
                  "100nM Dexamethasone (HC1,2,3,4,5,6, SURF1_1,2,3)" = "DEX_21",
                  "3% Oxygen + 100nM Dexamethasone (HC1,2,4)" = "DEX_3",
                  "mitoNUITS (HC1,2,3)" = "Modulators_21",
                  "mitoNUITS + 100nM Dexamethasone (HC1,2,3)" = "Modulators+DEX_21",
                  "1nM Oligomycin (HC1,2,3)" = "Oligomycin_21",
                  "1nM Oligomycin + 100nm Dexamethasone (HC1,2,3)" = "Oligomycin+DEX_21",
                  "1ug/mL 5-azacytidine (HC1,2)" = "5-azacytidine_21",
                  "1ug/mL 5-azacytidine + mitoNUITS (HC1,2)" = "5-azacytidine+Modulators_21",
                  "Pulsated 100nM Dexamethasone (30min/passage) (HC5,6)" = "Pulsated_DEX_21",
                  "10nM MitoQ (HC5,6)" = "MitoQ_21",
                  "10nM MitoQ + 100nM Dexamethasone (HC5,6)" = "MitoQ+DEX_21",
                  "5uM n-acetylcysteine (HC5,6)" = "NAC_21",
                  "5uM n-acetylcysteine + 100nM Dexamethasone (HC5,6)" = "NAC+DEX_21",
                  "1mM a-ketogluturate (HC5,6)" = "a-ketogluturate_21",
                  "10mM Betahydroxybutyrate (HC1,2)" = "betahydroxybutyrate_21",
                  "1mM 2-Deoxyglucose (HC1,2)" = "2-Deoxyglucose_21",
                  "5.5mM Galactose (HC1,2)" = "Galactose_21"
                ),
                multiple=TRUE, 
                selectize=TRUE, 
                selected="Control_21") #default value
  })
  
  #get the selected Cell Lines
  Selected_Cell_Lines <-reactive({
    if(is.null(input$Cell_Line) || length(input$Cell_Line)==0)
      return()
    as.vector(input$Cell_Line)
    
  })
  
  #get the selected Treatments
  Selected_Treatments <-reactive({
    if(is.null(input$Treatment) || length(input$Treatment)==0)
      return()
    as.vector(input$Treatment)
    
  })
  
  #filter the data according to the selected cell lines, treatments, and parameters
  Selected_Data <- reactive({
    Data <- LS_Data[LS_Data$Treatment %in% Selected_Treatments(),]
    Data <- Data[Data$Cell_Line %in% Selected_Cell_Lines(),]
    if(input$replicates_1 == FALSE) {
      Data <- Data[!(Data$Replicate_Line %in% 1),]
    }
    if (input$replicates_2 == FALSE) {
      Data <- Data[!(Data$Replicate_Line %in% 2),]
    }
    if (input$replicates_3 == FALSE) {
      Data <- Data[!(Data$Replicate_Line %in% 3),]
    }
    else {
      Data <- Data[Data$Replicate_Line %in% c(0,1,2,3),]
    }

  })
  

  # Main Function for Plotting a parameter over time.
  plotTime <- function(Data, parameter, Y_title) {
    scale <- input$scale
    outlier <- input$outliers
    SE <- input$se
    DoF <- input$poly
    second_axis <- input$second_axis
    
    x_data <- select_x_max(Data)
    x_max <- x_data[1]
    x_step <- x_data[2]
    
    if(input$Norm_Cell_Volume == TRUE) {
      cell_volume <- Data$Cell_Volume
      parameter <- parameter / cell_volume
      Y_title <- paste(Y_title, "/cell volume", sep = "")
    }
    
    if(input$Division_Norm == TRUE & input$tabs != "Seahorse Bioenergetics" & input$tabs != "Correlations") {
      div_per_passage <- Data$Divisions_Per_Passage
      div_per_passage[round(div_per_passage,1) == 0.0] <- NA 
      parameter <- parameter / div_per_passage
      #parameter[!is.nan(parameter)] <- NA
      Y_title <- paste(Y_title, "/division")
    }
    
    if(input$Age_Norm == TRUE) {
      parameter <- parameter / Data$Donor_Age_with_gestation
      Y_title <- paste(Y_title, " normalized to the donor's age")
    }
    
    if(scale == "log") {
      parameter <- log10(parameter)
      #parameter[round(parameter,1) == 0.0] <- NA 
      parameter[is.infinite(parameter)] <- NA 
      Y_title <- paste("Log10 ", Y_title, sep="")
    }
    
    Y_title <- paste(Y_title, ")", sep = "")
    
    outliers <- NULL
    if(outlier == "yes") {
      outliers <- boxplot(parameter, plot=FALSE)$out
      outliers_data <- Data[(parameter %in% outliers),]
      Data <- Data[!(parameter %in% outliers),]
      parameter <- parameter[!(parameter %in% outliers)]
    }
    
    
    annotation <- ""
    if(!is.null(outliers)) {
      n <- length(outliers_data$Days_Grown)
      tag <- rep("days", each=n)
      outliers_days_grown <- paste(round(outliers_data$Days_Grown,0), tag, sep = " ")
      concate_cols <- paste(outliers_data$Cell_line_group_new, outliers_days_grown, sep = ": ")
      outlier_string <- paste(concate_cols , collapse = ", ")
      annotation <- paste("Extreme Outlier(s) Removed: ", outlier_string)
    }
    
    formula <- y ~ poly(x, DoF, raw = TRUE)
    
    x_axis <- "Days_Grown"
    x_title <- "Time Grown (Days)"
    if(input$x_axis == "miage_days_grown") {
      x_axis <- "MiAge_Days_Grown"
      x_title <- "MiAge adjusted Days Grown"
    }
    if(input$x_axis == "doublings_axis") {
      x_axis <- "Population_Doublings_DI"
      x_title <- "Population Doublings"
      x_max <- 80
      x_step <- 10
      
    }
    if(input$x_axis == "miage_doublings_axis") {
      x_axis <- "MiAge_Population_Doublings"
      x_title <- "MiAge adjusted Population Doublings"
      x_max <- 80
      x_step <- 10
      
    }
    if(input$x_axis == "passage_axis") {
      x_axis <- "Passage"
      x_title <- "Passage"
      x_max <- 50
      x_step <- 5
      
    }
    
    treatments <- sort(Selected_Treatments())
    if(is.null(treatments)) {
      colors <- 'dimgray'
      colors_2 <- 'dimgray'
      shapes <- 21
    }
    else {
      colors <- getColors(treatments)
      colors_2 <- getColors_2(sort(Selected_Cell_Lines()))
      shapes <- getShapes(sort(Selected_Cell_Lines()))
    }

   
    
    p <- ggplot(Data, aes_string(x=x_axis, y = parameter, color ="Treatment")) +
      geom_point(size = 2 , stroke = 0.5, aes(shape = Cell_line_new, color = Treatment, fill = Cell_line_new)) +
      #geom_text(aes(label=Percent_Change),hjust=-0.2, vjust=1.2, size = 3) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors_2) +
      scale_shape_manual(values=shapes) +
      scale_y_continuous(name = Y_title) +
      scale_x_continuous(name =  x_title) +
      theme_classic() +
      annotate("text", label = annotation, x = max(as.numeric(Data$Days_Grown),na.rm=T)/2, y = 0, size = 1, hjust = 1, vjust = 1) +
      coord_cartesian(clip = "off")
    groups <- unique(Data$Cell_line_group_new)
    
    if(input$legend == TRUE) {
        p = p + theme(legend.title = element_blank())
    }
    else {
      p = p + theme(legend.position='none')
    }
    
    
    if(input$fit == TRUE & input$tabs != "Correlations") { 
      p = p + geom_smooth(method = "lm", formula = formula, size =0.4, aes(group = Cell_line_group_new),show.legend = FALSE, alpha = 0.05, se = SE)
  
    }

    if(second_axis == TRUE) {
      SF <- mean(Data$Population_Doublings_DI, na.rm = TRUE) / mean(parameter, na.rm = TRUE)
      Popdata <- Data$Population_Doublings_DI / SF
      p = p + geom_point(size = 2 , stroke = 0.5, aes(y = Popdata, shape = Cell_line_new, color = Treatment, fill = Cell_line_new))
      p = p + scale_y_continuous(name = Y_title, sec.axis = sec_axis(~.*SF, name = "Population Doublings"))
    }
    
    if(input$tabs != "Correlations") {
      p <-  ggplotly(p, tooltip = c("Cell_line_new", "Treatment", "y", "Days_Grown")) %>%
        config(displaylogo = FALSE)
      if(input$annotation != "none") {
        startCount <- 1
        count <- length(unique(Data$Cell_line_group_new))
        for (i in 1:count) {
          group <- unique(Data$Cell_line_group_new)[i]
          endCount <-  startCount + nrow(Data[Data$Cell_line_group_new %in% group,]) - 1
          
          time_points <- Data[startCount:endCount,]$Days_Grown
          if(input$x_axis == "miage_time_grown") {
            time_points <- Data[startCount:endCount,]$MiAge_Days_Grown
          }
          if(input$x_axis == "doublings_axis") {
            time_points <- Data[startCount:endCount,]$Population_Doublings_DI
          }
          if(input$x_axis == "miage_doublings_axis") {
            time_points <- Data[startCount:endCount,]$MiAge_Population_Doublings
          }
          if(input$x_axis == "passage_axis") {
            time_points <- Data[startCount:endCount,]$Passage
          }
          group_data <- parameter[startCount:endCount]
          
          Eq_Text <- as.character(group)
          
          if(input$annotation == "all") {
            model <- lm(group_data ~ poly(time_points, DoF, raw = TRUE))
            eq <- as.character(signif(as.polynomial(coef(model)), 2))
            r2 <- summary(model)$r.squared
            Eq_Text <- paste(group, ":\n y = ", eq, "\n r^2: ",round(r2,2), sep = "")
            
          }
         
          
          xpos <- tail(time_points[!is.na(group_data)], n=1)
          #xpos <- xpos + (0.15*xpos)
          ypos <- tail(group_data[!is.na(group_data)], n=1)
          #ypos <- ypos - (0.05*ypos)
          
          treatment <- Data[Data$Cell_line_group_new == group,]$Treatment[1]
          color <- 'dimgray'
          if(is.null(treatments)) {
            color <- 'dimgray'
          }
          else {
            color <- getColors(treatment)
          }
          
          
          p = p %>% add_annotations(text = Eq_Text, x=xpos, y=ypos, font = list(size = 8, color = color), xanchor = "left", yanchor = "top", showarrow = F) 
          startCount <- endCount
        }
      }
      }
    
    return(p)
  }
  
  #Cell Line list
  output$DatasetSelector<-renderUI({
    selectInput('Dataset', 'Dataset',
                c("Cytology" = "cytology",
                  "DNAmAges" = "DNAm",
                  "Seahorse Bioenergetics" = "bioenergetics",
                  "Telomere Length" = "telomere_length",
                  "mtDNA copy number" = "mtDNA_copy_number",
                  "mtDNA sequencing"= "mtDNA_seq",
                  "RNAseq metadata" = "RNAseq",
                  "Cytokines" = "cytokines",
                  "IL6 & GDF15 ELISAs" = "ELISAs",
                  "cell-free DNA" = "cfDNA",
                  "All data" = "all"), 
                multiple=FALSE, 
                selectize=TRUE, 
                selected="all") #default value
  })
  
  #get the selected Cell Lines
  Selected_Dataset<-reactive({
    if(is.null(input$Dataset) || length(input$Dataset)==0)
      return()
    as.vector(input$Dataset)
    
  })
  
  # Download selected Data
  output$downloadData <- downloadHandler(
    filename = paste0("Lifespan_Study_selected",Selected_Dataset(),"_data.csv"),
    content = function(file) {
      metadata <- read.csv(paste0("downloadable_data/Cellular_lifespan_study_",Selected_Dataset(),".csv"))
      output_data <- Selected_Data()
      output_data <- metadata[metadata$Sample %in% output_data$Sample,]
      if(Selected_Dataset() != "all") {
        output_data <- output_data[!is.na(output_data[,length(output_data)]),]
      }
      write.csv(output_data, file, na = "", row.names = FALSE)
    }
  )
  
  # save graphs
  observeEvent(input$saveGraphs, {
   filename =  "Lifespan_Study_Graphs.pdf"
   plist <- list(Growth_Plot(), Cell_Size_Plot(), Percent_Dead_Plot(),TL_Plot(), CN_Plot(), cfDNA_Plot(), IL6_Plot())
   plots <- subplot(plist, nrows = length(plist), titleX = TRUE, titleY = TRUE)
   orca(plots,filename, height = (250)*length(plist))
  
  })
  
  

  #### Growth Curves ###
  ############  
  
  
  # Insert Growth Curve Plots
  output$Growth_Plots <- renderPlotly({
    
    Data <- Selected_Data()
    parameter <- Data$Population_Doublings_DI
    Y_title <- "Population Doublings (divisions"
    if(input$adj_population_doublings == TRUE) {
      parameter <- Data$MiAge_Population_Doublings
      Y_title <- "MiAge adjusted Population Doublings (divisions"
    }

    p <- plotTime(Data, parameter,Y_title)
   p
    
  })
  
  
  Growth_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$Population_Doublings_DI
    Y_title <- "Population Doublings (divisions"
    p <- plotTime(Data, parameter,Y_title)
    p
  })
  
  
  
  # Insert Growth Curve Plots
  output$Growth_Rate_Plots <- renderPlotly({
    
    
    Data <- Selected_Data()
    parameter <- Data$Doubling_Rate
    Y_title <- "Doubling Rate (divisions / day"
    
    
    p <- plotTime(Data, parameter,Y_title)
    p
    
  })
  
  
  Growth_Rate_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$Doubling_Rate
    Y_title <- "Doubling Rate (divisions / day"
    
    p <- plotTime(Data, parameter,Y_title)
    p
  })
  
  # observe({
  #   if(input$tabs == "Growth Curves") {
  #     updateSliderInput(session, "poly", value = 3)
  #   }
  # })
  
  
  # Insert Cell Size Plots
  output$Cell_Size_Plots <- renderPlotly({
    Data <- Selected_Data()
    parameter <- Data$Cell_Volume
    Y_title <- "Cell Volume (um^3"
    
    if(input$diameter_conversion == TRUE) {
      parameter <- Data$Cell_Size
      Y_title <- "Cell Diameter (um"
    }
    
    p <- plotTime(Data, parameter,Y_title)
    p
    
  })
  
  Cell_Size_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$Cell_Volume
    Y_title <- "Cell Volume (um^3"
    
    if(input$diameter_conversion == TRUE) {
      parameter <- Data$Cell_Size
      Y_title <- "Cell Diameter (um"
    }
    
    p <- plotTime(Data, parameter,Y_title)
    p
  })
  # 
  # output$downloadPlot_Cell_Size <- downloadHandler(
  #   filename = paste('Cell_Size_Plot.', "jpeg", sep = ""),
  #   content = function(file) {
  #     ggsave(file, plot = Cell_Size_Plot(), width = 20, height = 12, units = "cm", device = "jpeg")
  # })
  
  # Insert Percent Dead Plots
  output$Percent_Dead_Plots <- renderPlotly({
    
    Data <- Selected_Data()
    parameter <- Data$Percent_Dead
    Y_title <- "Dead Cells (%"
  
    
    p <- plotTime(Data, parameter,Y_title)
    p
    
  })
  

  Percent_Dead_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$Percent_Dead
    Y_title <- "Dead Cells (%"
    p <- plotTime(Data, parameter,Y_title)
    p
  })
  #############
  
  #### DNAmAge ###
  ############  
  
  
  #Parameter list
  output$ClockSelector<-renderUI({
    selectInput('Clock_Parameter', 'DNA Methylation Age Clocks',
                Clock_Parameters, 
                multiple=FALSE, 
                selectize=TRUE, 
                selected="Skin_Blood_Clock") #default value
  })
  
  #get the selected Parameters
  Selected_Clock <-reactive({
    if(is.null(input$Clock_Parameter) || length(input$Clock_Parameter)==0)
      return()
    as.vector(input$Clock_Parameter)
  })
  
  

  output$DNAmAge_Plots <- renderPlotly({     
        Data <- Selected_Data()
        parameter <- as.character(Selected_Clock())
        Y_title <- "Epigenetic Age ("
        parameter_col <- as.numeric(which(names(Data)==parameter))
        parameter <- as.numeric(Data[,parameter_col])
        Y_title <- paste(Y_title, as.character(Selected_Clock()), sep = "") 
        
        
        p <- plotTime(Data, parameter,Y_title)
        p
    })


  # 
  # DNAmAge_Plot <- reactive({
  #   Data <- Selected_Data()
  #   parameter <- as.character(Selected_Clock())
  #   Y_title <- "Epigenetic Age ("
  #   parameter_col <- as.numeric(which(names(Data)==parameter))
  #   parameter <- as.numeric(Data[,parameter_col])
  #   Y_title <- paste(Y_title, as.character(Selected_Clock()), sep = " ") 
  #   Y_title <- paste(Y_title, ")", sep = "")
  #   
  #   
  #   p <- plotTime(Data, parameter,Y_title)
  #   p
  # })
  # 
  # output$downloadPlot_DNAm <- downloadHandler(
  #   filename = paste('Epigenetic_Age_Plot.', "jpeg", sep = ""),
  #   content = function(file) {
  #     ggsave(file, plot = DNAmAge_Plot(), width = 20, height = 12, units = "cm", device = "jpeg")
  #   })
 
  ############# 
  
  ### DNA Methylation CpGs ###
  ############
  
  observe({
     if(input$tabs == "DNA Methylation CpGs") {
       withProgress(message = 'Loading Methylation Data (~1 min)', value = 0.25, {
           load("betas.rcp.FunNorm.RData") # 1 min load time
       })
        CpGs <- as.vector(rownames(beta.rcp))
        updateSelectizeInput(session, "CpG", choices = CpGs, server = TRUE, selected = "cg14817997")
        
        # Insert CpG Time Course Plot
        output$CpG_Plot <- renderPlotly({
          
          Data <- Selected_Data()
          parameter_name <-  as.character(input$CpG)
          Y_title <- paste("CpG ",parameter_name, " (",sep="")
          
          #parameter_name <- "cg14817997"
          #Data <- LS_Data
          Data$basename <- as.character(Data$basename)
          
          parameter_row <- as.numeric(which(CpGs==parameter_name))
          parameter <- as.data.frame(beta.rcp[parameter_row,])
          colnames(parameter) <- parameter_name
          CpG <- as.numeric(matrix(ncol = 1, nrow = nrow(Data)))
          Data <- cbind(Data, CpG)
          
          
          for(i in 1:nrow(Data)) {
            if(Data$basename[i] %in% rownames(parameter)) {
              Data$CpG[i] <- parameter[which(rownames(parameter)==Data$basename[i]),1]  
            }
          }
          parameter <- as.numeric(Data$CpG)
          print(parameter)
          p <- plotTime(Data, parameter,Y_title)
          
          #p <- p %>% layout(yaxis = list(range = c(0,1.0), dtick = 0.20))
          p 

          
        })
        
        # Insert PCA Plot
        output$PCA_Plot <- renderPlotly({
          input$PCA_button
          isolate({
            Data <- Selected_Data()
            Data$basename <- as.character(Data$basename)
            betas <- beta.rcp[,colnames(beta.rcp) %in% Data$basename]
            basenames<-Data$basename[Data$basename != ""]
            betas <- betas[,match(basenames,colnames(betas))]
            
            withProgress(message = 'Loading PCA (~2 min)', value = 0.25, {
              pca = prcomp(t(betas), center = TRUE, scale = TRUE)
            })
            
            pca_x_data <- as.numeric(pca$x[,1])
            pca_y_data <- as.numeric(pca$x[,2])
            
            
            parameter_x_data <- rep(NA, nrow(Data))
            parameter_y_data <- rep(NA, nrow(Data))
            count <- 1
            for(i in 1:nrow(Data)) {
              if(Data$basename[i] == "") {
                parameter_x_data[i] <- NA
                parameter_y_data[i] <- NA
              }
              else {
                parameter_x_data[i] <- pca_x_data[count]
                parameter_y_data[i] <- pca_y_data[count]
                count <- count + 1
              }
              
            }
            
            sum <- summary(pca)
            PoV1 <- round(sum$importance[2]*100,2) # Proportion of Variance Explained for PC1
            PoV2 <- round(sum$importance[5]*100,2) # Proportion of Variance Explained for PC2
            
            X_title <- paste("Principle Component 1 (", PoV1, "% variance)", sep="")
            Y_title <- paste("Principle Component 2 (", PoV2, "% variance)", sep="")
            
            
            outliers_data <- NULL
            if(input$outliers == "yes") {
              outliers_x <- boxplot(parameter_x_data, plot=FALSE)$out
              #Data[(parameter_x_data %in% outliers_x),] <- NA
              parameter_x_data[(parameter_x_data %in% outliers_x)] <- NA
              outliers_data_x <- Data[(parameter_x_data %in% outliers_x),]
              
              outliers_y <- boxplot(parameter_y_data, plot=FALSE)$out
              #Data[(parameter_y_data %in% outliers_y),] <- NA
              parameter_y_data[(parameter_y_data %in% outliers_y)] <- NA
              outliers_data_y <- Data[(parameter_y_data %in% outliers_y),]
              
              outliers_data <- merge(outliers_data_x, outliers_data_y, by = "Unique_Variable_Name")
              
            }
            
            annotation <- ""
            if(!is.null(outliers_data)) {
              n <- length(outliers_data$Days_Grown)
              tag <- rep("days", each=n)
              outliers_days_grown <- paste(round(as.numeric(outliers_data$Days_Grown),0), tag, sep = " ")
              concate_cols <- paste(outliers_data$Cell_Line_Group, outliers_days_grown, sep = ": ")
              outlier_string <- paste(concate_cols , collapse = ", ")
              annotation <- paste("Extreme Outlier(s) Removed: ", outlier_string)
            }
            
            treatments <- sort(Selected_Treatments())
            if(is.null(treatments)) {
              colors <- 'dimgray'
              colors_2 <- 'dimgray'
              shapes <- 21
              
            }
            else {
              colors <- getColors(treatments)
              colors_2 <- getColors_2(sort(Selected_Cell_Lines()))
              shapes <- getShapes(sort(Selected_Cell_Lines()))
            }
            scale <- input$scale
            SE <- input$se
            DoF <- input$poly
            formula <- y ~ poly(x, DoF, raw = TRUE)
            
            p <- ggplot(Data, aes_string(x=parameter_x_data, y=parameter_y_data, color ="Treatment")) +
              scale_color_manual(values = colors) +
              scale_fill_manual(values = colors_2) +
              scale_shape_manual(values=shapes) +
              scale_y_continuous(name = Y_title) +
              scale_x_continuous(name =  X_title) +
              theme_classic() +
              annotate("text", label = annotation, x = 2, y = 0, size = 3, hjust = 0) +
              coord_cartesian(clip = "off")+
              theme(text = element_text(size = 14),
                    legend.position="none",
                    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                    plot.margin = margin(t = 6, r = 6, b = 6, l = 6, unit = "pt"),
                    plot.title = element_text(size = 32, hjust = 0.05, vjust = -.1))
            
            groups <- unique(Data$Cell_Line_Group)
            
            
            p = p + geom_point(size = 2 , stroke = 0.5, aes(shape = Cell_Line, color = Treatment, fill = Cell_Line))
            p = p + geom_text(aes(label=paste(round(Data$Days_Grown, digits = 0), "days", sep = " ")),hjust=-1, vjust=-1, size = 2, alpha = 0.7)
            p = p + geom_hline(yintercept=0)
            p = p + geom_vline(xintercept=0)
            
            
            # if(input$fit == TRUE) { 
            #   p = p + geom_smooth(method = "lm", formula = formula, size = 1.5, aes(group = Cell_Line_Group), alpha = 0.05, se = SE)
            # }
            
            # Convert 2 Graphs to Plotly and arrange in grid
            p <- ggplotly(p)
            
            
            # Add equation info to first and second graph
            # if (input$fit_text == TRUE) {
            #   count <- length(unique(Data$Cell_Line_Group))
            #   p = p + stat_poly_eq(formula = formula, aes(group = Cell_Line_Group, label =  paste(groups, stat(eq.label), stat(rr.label), sep = "~~")),
            #                        rr.digits = 2, coef.digits = 2,parse = TRUE,  label.x = "left", label.y = "bottom", show.legend = TRUE)
            # }
            
            
            p
          })
          
        })
        
        
        
     }
  })
  
  #Parameter list
  # output$CpGSelector<-renderUI({
  #   selectInput('CpG_Parameter', 'CpG Site',
  #               CpGs, 
  #               multiple=FALSE, 
  #               selectize=TRUE, 
  #               selected="cg03278726") #default value
  # })
  # 
  # #get the selected Parameters
  # Selected_CpG <-reactive({
  #   if(is.null(input$CpG_Parameter) || length(input$CpG_Parameter)==0)
  #     return()
  #   as.vector(input$CpG_Parameter)
  # })
  
  
  
  
  
  ############  
  
  ### RNAseq Genes ###
  ############
  
  observe({
    if(input$tabs == "RNAseq Genes") {
      withProgress(message = 'Loading RNAseq Data (~1 min)', value = 0.25, {
        #RNAseq_Data <- read.csv("allRNA_Transcript_RNAseq_data_no_cutoff.csv", header=TRUE, row.names = 1) # 1 min load time
        #RNAseq_Gene_Data <- read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", header=TRUE, row.names = 1)
        #setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp")
        RNAseq_Gene_Data <- read.csv("RNAseq_allGenes_vst_normalized.csv",header=T,row.names=1)
        # Rename Samples
        corrected_sample_name <- substring(colnames(RNAseq_Gene_Data),str_locate(colnames(RNAseq_Gene_Data), "X")[,1]+1)
        colnames(RNAseq_Gene_Data) <- corrected_sample_name
      })
     
      Genes <- as.vector(rownames(RNAseq_Gene_Data))
      updateSelectizeInput(session, "Gene", choices = Genes, server = TRUE, selected = "UCP2")

      
      # Insert Gene Time Course Plot
      output$RNA_Gene_Plot <- renderPlotly({
        Data <- Selected_Data()
        parameter_name <-  as.character(input$Gene)
        Y_title <- paste("Gene ", parameter_name, " (Log2 normalized exp",sep="")
        
        if(input$log2_genes == F) {
          RNAseq_Gene_Data <- 2^RNAseq_Gene_Data
          Y_title <- paste("Gene ", parameter_name, " (normalized exp",sep="")
        }
        

       
        #parameter_name <- "lcl|NC_000001.10;NM_001160234.1;ATP1A1;mRNA"
        #Data <- LS_Data
        Data$RNAseq_ID <- as.character(Data$RNAseq_ID)
        Data <- Data[!is.na(Data$RNAseq_ID),]
        Data <- Data[Data$RNAseq_ID %in% colnames(RNAseq_Gene_Data),]
        RNAseq_Data_filtered <- RNAseq_Gene_Data[,colnames(RNAseq_Gene_Data) %in% Data$RNAseq_ID]
        RNAseq_Data_filtered <- RNAseq_Data_filtered[,match(Data$RNAseq_ID, colnames(RNAseq_Data_filtered))]
        
        
        parameter_row <- as.numeric(which(Genes==parameter_name))
        parameter <- RNAseq_Data_filtered[parameter_row,]
        

        
        colnames(parameter) <- parameter_name
        # Gene <- as.numeric(matrix(ncol = 1, nrow = nrow(Data)))
        parameter <- as.numeric(parameter)
        Data <- cbind(Data, parameter)
        
        
        #print(parameter)
        p <- plotTime(Data, parameter, Y_title)
        
        #p <- p %>% layout(yaxis = list(range = c(0,1.0), dtick = 0.20))
        p 
        
        
      })
      
      # Transcripts <- as.vector(rownames(RNAseq_Data))
      # updateSelectizeInput(session, "Transcript", choices = Transcripts, server = TRUE, selected = "lcl|NC_000011.9;NM_003355.2;UCP2;mRNA")
      # 
      
      # Insert Transcript Time Course Plot
      # output$RNA_Transcript_Plot <- renderPlotly({
      #   
      #   Data <- Selected_Data()
      #   parameter_name <-  as.character(input$Transcript)
      #   title_gene_name <- substring(parameter_name, str_locate(parameter_name, ";")[,1]+1)
      #   title_gene_name <- substring(title_gene_name, str_locate(title_gene_name, ";")[,1]+1)
      #   Y_title <- paste("Transcript ", title_gene_name, " (TPM",sep="")
      # 
      #   #parameter_name <- "lcl|NC_000001.10;NM_001160234.1;ATP1A1;mRNA"
      #   #Data <- LS_Data
      #   Data$RNAseq_ID <- as.character(Data$RNAseq_ID)
      #   Data <- Data[!is.na(Data$RNAseq_ID),]
      #   Data <- Data[Data$RNAseq_ID %in% colnames(RNAseq_Data),]
      #   RNAseq_Data_filtered <- RNAseq_Data[,colnames(RNAseq_Data) %in% Data$RNAseq_ID]
      #   RNAseq_Data_filtered <- RNAseq_Data_filtered[,match(Data$RNAseq_ID, colnames(RNAseq_Data_filtered))]
      # 
      # 
      #   parameter_row <- as.numeric(which(Transcripts==parameter_name))
      #   parameter <- RNAseq_Data_filtered[parameter_row,]
      #   
      #   colnames(parameter) <- parameter_name
      #   # Gene <- as.numeric(matrix(ncol = 1, nrow = nrow(Data)))
      #   parameter <- as.numeric(parameter)
      #   Data <- cbind(Data, parameter)
      #   
      #   #print(parameter)
      #   p <- plotTime(Data, parameter, Y_title)
      #   
      #   #p <- p %>% layout(yaxis = list(range = c(0,1.0), dtick = 0.20))
      #   p 
      #   
      #   
      # })
      # 
      # Insert PCA Plot
      output$PCA_Plot_2D_RNA <- renderPlot({
        input$PCA_button_2D_RNA
        isolate({
          
          Data <- Selected_Data()
          #Data <- LS_Data
          Data$RNAseq_ID <- as.character(Data$RNAseq_ID)
          Data <- Data[!is.na(Data$RNAseq_ID),]
          Data <- Data[Data$RNAseq_ID %in% colnames(RNAseq_Gene_Data),]
          RNAseq_Data_filtered <- RNAseq_Gene_Data[,colnames(RNAseq_Gene_Data) %in% Data$RNAseq_ID]
          RNAseq_Data_filtered <- RNAseq_Data_filtered[,match(Data$RNAseq_ID, colnames(RNAseq_Data_filtered))]
          rna <- as.matrix(t(RNAseq_Data_filtered))
          rna <- rna[,which(apply(rna, 2, var)!=0)]
          withProgress(message = 'Loading PCA (~2 min)', value = 0.25, {
            pca = prcomp(rna, center = TRUE, scale = TRUE)
          })
          
          pca_x_data <- as.numeric(pca$x[,1])
          pca_y_data <- as.numeric(pca$x[,2])
          pca_z_data <- as.numeric(pca$x[,3])
          
          parameter_x_data <- rep(NA, nrow(Data))
          parameter_y_data <- rep(NA, nrow(Data))
          parameter_z_data <- rep(NA, nrow(Data))
          count <- 1
          for(i in 1:nrow(Data)) {
            if(Data$RNAseq_ID[i] == "") {
              parameter_x_data[i] <- NA
              parameter_y_data[i] <- NA
              parameter_z_data[i] <- NA
            }
            else {
              parameter_x_data[i] <- pca_x_data[count]
              parameter_y_data[i] <- pca_y_data[count]
              parameter_z_data[i] <- pca_z_data[count]
              count <- count + 1
            }
            
          }
          
          sum <- summary(pca)
          PoV1 <- round(sum$importance[2]*100,2) # Proportion of Variance Explained for PC1
          PoV2 <- round(sum$importance[5]*100,2) # Proportion of Variance Explained for PC2
          PoV3 <- round(sum$importance[8]*100,2) # Proportion of Variance Explained for PC3
          
          X_title <- paste("Principle Component 1 (", PoV1, "% variance)", sep="")
          Y_title <- paste("Principle Component 2 (", PoV2, "% variance)", sep="")
          Z_title <- paste("Principle Component 3 (", PoV2, "% variance)", sep="")
          
          
          outliers_data <- NULL
          if(input$outliers == "yes") {
            outliers_x <- boxplot(parameter_x_data, plot=FALSE)$out
            #Data[(parameter_x_data %in% outliers_x),] <- NA
            parameter_x_data[(parameter_x_data %in% outliers_x)] <- NA
            outliers_data_x <- Data[(parameter_x_data %in% outliers_x),]
            
            outliers_y <- boxplot(parameter_y_data, plot=FALSE)$out
            #Data[(parameter_y_data %in% outliers_y),] <- NA
            parameter_y_data[(parameter_y_data %in% outliers_y)] <- NA
            outliers_data_y <- Data[(parameter_y_data %in% outliers_y),]
            
            outliers_z <- boxplot(parameter_z_data, plot=FALSE)$out
            #Data[(parameter_y_data %in% outliers_z),] <- NA
            parameter_z_data[(parameter_z_data %in% outliers_z)] <- NA
            outliers_data_z <- Data[(parameter_z_data %in% outliers_z),]
            
            outliers_data <- merge(outliers_data_x, outliers_data_y, by = "Unique_Variable_Name")
            
          }
          
          annotation <- ""
          if(!is.null(outliers_data)) {
            n <- length(outliers_data$Days_Grown)
            tag <- rep("days", each=n)
            outliers_days_grown <- paste(round(as.numeric(outliers_data$Days_Grown),0), tag, sep = " ")
            concate_cols <- paste(outliers_data$Cell_Line_Group, outliers_days_grown, sep = ": ")
            outlier_string <- paste(concate_cols , collapse = ", ")
            annotation <- paste("Extreme Outlier(s) Removed: ", outlier_string)
          }
          
          treatments <- sort(Selected_Treatments())
          if(is.null(treatments)) {
            colors <- 'dimgray'
            colors_2 <- 'dimgray'
            shapes <- 21
            
          }
          else {
            colors <- getColors(treatments)
            colors_2 <- getColors_2(sort(Selected_Cell_Lines()))
            shapes <- getShapes(sort(Selected_Cell_Lines()))
          }
          scale <- input$scale
          SE <- input$se
          DoF <- input$poly
          formula <- y ~ poly(x, DoF, raw = TRUE)
          
          p <- ggplot(Data, aes_string(x=parameter_x_data, y=parameter_y_data, size = "Days_Grown", color ="Treatment")) +
            scale_color_manual(values = colors) +
            scale_fill_manual(values = colors_2) +
            scale_shape_manual(values = shapes) +
            scale_y_continuous(name = Y_title) +
            scale_x_continuous(name =  X_title) +
            theme_classic() +
            annotate("text", label = annotation, x = 2, y = 0, size = 3, hjust = 0) +
            coord_cartesian(clip = "off")+
            theme(text = element_text(size = 14),
                  #legend.position="none",
                  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  plot.margin = margin(t = 6, r = 6, b = 6, l = 6, unit = "pt"),
                  plot.title = element_text(size = 32, hjust = 0.05, vjust = -.1))

          groups <- unique(Data$Cell_line_group_new)
          
          
          p = p + geom_point(stroke = 2, aes(size = Days_Grown, shape = Cell_line_new, color = Treatment, fill = Cell_line_new))
          p = p + scale_size(range = c(1, 10), name="Days Grown")
          p = p + geom_text(aes(label=paste(round(Days_Grown, digits = 0), "days", sep = " ")),hjust=-0.5, vjust=-0.5, size = 3, alpha = 0.9)
          p = p + geom_hline(yintercept=0)
          p = p + geom_vline(xintercept=0)


          # if(input$fit == TRUE) {
          #    p = p + geom_smooth(method = "lm", formula = formula, size = 1.5, aes(group = Cell_Line_Group), alpha = 0.05, se = SE)
          # }

          # Convert 2 Graphs to Plotly and arrange in grid
          #p <- ggplotly(p)
          # p <- p %>% add_markers(size = ~Data$Days_Grown, sizes= c(10,50), opacity=1)
          
            
          X_title <- paste("PC1 ", PoV1, "%", sep="")
          Y_title <- paste("PC2 ", PoV2, "%", sep="")
          Z_title <- paste("PC3 ", PoV3, "%", sep="")
      
    

          # Add equation info to first and second graph
          # if (input$fit_text == TRUE) {
          #   count <- length(unique(Data$Cell_Line_Group))
          #   p = p + stat_poly_eq(formula = formula, aes(group = Cell_Line_Group, label =  paste(groups, stat(eq.label), stat(rr.label), sep = "~~")),
          #                        rr.digits = 2, coef.digits = 2,parse = TRUE,  label.x = "left", label.y = "bottom", show.legend = TRUE)
          # }
          
          
          p
        })
        
      })
      
      # Insert PCA Plot
      output$PCA_Plot_3D_RNA <- renderPlotly({
        input$PCA_button_3D_RNA
        isolate({
          
          Data <- Selected_Data()
          #Data <- LS_Data
          Data$RNAseq_ID <- as.character(Data$RNAseq_ID)
          Data <- Data[!is.na(Data$RNAseq_ID),]
          Data <- Data[Data$RNAseq_ID %in% colnames(RNAseq_Gene_Data),]
          RNAseq_Data_filtered <- RNAseq_Gene_Data[,colnames(RNAseq_Gene_Data) %in% Data$RNAseq_ID]
          RNAseq_Data_filtered <- RNAseq_Data_filtered[,match(Data$RNAseq_ID, colnames(RNAseq_Data_filtered))]
          rna <- as.matrix(t(RNAseq_Data_filtered))
          rna <- rna[,which(apply(rna, 2, var)!=0)]
          withProgress(message = 'Loading PCA (~2 min)', value = 0.25, {
            pca = prcomp(rna, center = TRUE, scale = TRUE)
          })
          
          pca_x_data <- as.numeric(pca$x[,1])
          pca_y_data <- as.numeric(pca$x[,2])
          pca_z_data <- as.numeric(pca$x[,3])
          
          parameter_x_data <- rep(NA, nrow(Data))
          parameter_y_data <- rep(NA, nrow(Data))
          parameter_z_data <- rep(NA, nrow(Data))
          count <- 1
          for(i in 1:nrow(Data)) {
            if(Data$RNAseq_ID[i] == "") {
              parameter_x_data[i] <- NA
              parameter_y_data[i] <- NA
              parameter_z_data[i] <- NA
            }
            else {
              parameter_x_data[i] <- pca_x_data[count]
              parameter_y_data[i] <- pca_y_data[count]
              parameter_z_data[i] <- pca_z_data[count]
              count <- count + 1
            }
            
          }
          
          sum <- summary(pca)
          PoV1 <- round(sum$importance[2]*100,2) # Proportion of Variance Explained for PC1
          PoV2 <- round(sum$importance[5]*100,2) # Proportion of Variance Explained for PC2
          PoV3 <- round(sum$importance[8]*100,2) # Proportion of Variance Explained for PC3
          
          X_title <- paste("Principle Component 1 (", PoV1, "% variance)", sep="")
          Y_title <- paste("Principle Component 2 (", PoV2, "% variance)", sep="")
          Z_title <- paste("Principle Component 3 (", PoV2, "% variance)", sep="")
          
          
          outliers_data <- NULL
          if(input$outliers == "yes") {
            outliers_x <- boxplot(parameter_x_data, plot=FALSE)$out
            #Data[(parameter_x_data %in% outliers_x),] <- NA
            parameter_x_data[(parameter_x_data %in% outliers_x)] <- NA
            outliers_data_x <- Data[(parameter_x_data %in% outliers_x),]
            
            outliers_y <- boxplot(parameter_y_data, plot=FALSE)$out
            #Data[(parameter_y_data %in% outliers_y),] <- NA
            parameter_y_data[(parameter_y_data %in% outliers_y)] <- NA
            outliers_data_y <- Data[(parameter_y_data %in% outliers_y),]
            
            outliers_z <- boxplot(parameter_z_data, plot=FALSE)$out
            #Data[(parameter_y_data %in% outliers_z),] <- NA
            parameter_z_data[(parameter_z_data %in% outliers_z)] <- NA
            outliers_data_z <- Data[(parameter_z_data %in% outliers_z),]
            
            outliers_data <- merge(outliers_data_x, outliers_data_y, by = "Unique_Variable_Name")
            
          }
          
          annotation <- ""
          if(!is.null(outliers_data)) {
            n <- length(outliers_data$Days_Grown)
            tag <- rep("days", each=n)
            outliers_days_grown <- paste(round(as.numeric(outliers_data$Days_Grown),0), tag, sep = " ")
            concate_cols <- paste(outliers_data$Cell_Line_Group, outliers_days_grown, sep = ": ")
            outlier_string <- paste(concate_cols , collapse = ", ")
            annotation <- paste("Extreme Outlier(s) Removed: ", outlier_string)
          }
          
          treatments <- sort(Selected_Treatments())
          if(is.null(treatments)) {
            colors <- 'dimgray'
            colors_2 <- 'dimgray'
            shapes <- 21
            
          }
          else {
            colors <- getColors(treatments)
            colors_2 <- getColors_2(sort(Selected_Cell_Lines()))
            shapes <- getShapes(sort(Selected_Cell_Lines()))
          }
          scale <- input$scale
          SE <- input$se
          DoF <- input$poly
          formula <- y ~ poly(x, DoF, raw = TRUE)
          
          X_title <- paste("PC1 ", PoV1, "%", sep="")
          Y_title <- paste("PC2 ", PoV2, "%", sep="")
          Z_title <- paste("PC3 ", PoV3, "%", sep="")
          
          plot_3d <- plot_ly(type="scatter3d", data = Data, x = parameter_x_data, y = parameter_y_data, z = parameter_z_data,
                             color = Data$Treatment, colors = colors, 
                             size = ~Data$Days_Grown,
                             sizes= c(5,30),
                             symbol = Data$Cell_line_new,
                             symbols = shapes,
                             text = ~paste("Group: ", Data$Cell_line_group_new, '$<br>Cut:', Data$Days_Grown),
                             marker = list(opacity = 1, sizemode = 'diameter')
                      )
          # plot_3d <- plot_3d %>% add_markers(size = ~Days_Grown, opacity=1)
          plot_3d <- plot_3d %>% layout(scene = list(xaxis = list(title = X_title),
                                                     yaxis = list(title = Y_title),
                                                     zaxis = list(title = Z_title)))
          
          plot_3d
        })
        
      })
    }
  })
  
  #Parameter list
  # output$CpGSelector<-renderUI({
  #   selectInput('CpG_Parameter', 'CpG Site',
  #               CpGs, 
  #               multiple=FALSE, 
  #               selectize=TRUE, 
  #               selected="cg03278726") #default value
  # })
  # 
  # #get the selected Parameters
  # Selected_CpG <-reactive({
  #   if(is.null(input$CpG_Parameter) || length(input$CpG_Parameter)==0)
  #     return()
  #   as.vector(input$CpG_Parameter)
  # })
  
  
  
  
  
  ############  
  
  
  #### Telomere Length ###
  ############  
    
    # Insert Telomere Length Plots
    output$Telomere_Plots <- renderPlotly({
      Data <- Selected_Data()
      parameter <- Data$Telomere_Length
      Y_title <- "Telomere Length (T/S"
      
      
      p <- plotTime(Data, parameter, Y_title)
      p
      
      
    })
    
  # # Allow for hovering/clicking of data points
  # output$TL_DataPoint <- renderText({
  #   parameter <- Selected_Data()$Telomere_Length
  #   xy_str <- function(e) {
  #     if(is.null(e)) return("Select Data Point for more information...\n")
  #     data_points <- Selected_Data()[which(abs(Selected_Data()$Days_Grown-e$x)%in% sort(abs(Selected_Data()$Days_Grown-e$x))[1:5]),]
  #     data_point <- data_points[which(abs(data_points$Telomere_Length-e$y) %in% sort(abs(data_points$Telomere_Length-e$y))[1]),]
  #     group <- data_point$Cell_Line_Group
  #     CV <- data_point$Telomere_Length_CV *100
  #     paste0("Cell Line Group: ", group, "\n" , "Time Grown: ", round(e$x, 1), " days","\n", "Telomer Length: ", round(e$y, 1), " T/S", "\n", "CV: ", CV, " %")
  #   }
  #   
  #   paste0(
  #     xy_str(input$time_plot_click)[1]
  #   )
  # })
  # 
  TL_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$Telomere_Length
    Y_title <- "Telomere Length (T/S"
    
    p <- plotTime(Data, parameter,Y_title)
    p
  })
  # 
  # output$downloadPlot_TL <- downloadHandler(
  #   filename = paste('Telomere_Length_Plot.', "jpeg", sep = ""),
  #   content = function(file) {
  #     ggsave(file, plot = TL_Plot(), width = 20, height = 12, units = "cm", device = "jpeg")
  #   })
  
  
  #############  
  
  ### Copy Number ###
  ############   
  
  # Insert cf-DNA Plots
  output$CN_Plots <- renderPlotly({
    Data <- Selected_Data()
    parameter <- Data$Copy_Number
    Y_title <- "Copy Number (copies/cell"
    
    
    p <- plotTime(Data, parameter,Y_title)
    p
    
  })
  # 
  # # Allow for hovering/clicking of data points
  # output$CN_DataPoint <- renderText({
  #   
  #   xy_str <- function(e) {
  #     if(is.null(e)) return("Select Data Point for more information...\n")
  #     data_points <- Selected_Data()[which(abs(Selected_Data()$Days_Grown-e$x)%in% sort(abs(Selected_Data()$Days_Grown-e$x))[1:5]),]
  #     data_point <- data_points[which(abs(data_points$Copy_Number-e$y) %in% sort(abs(data_points$Copy_Number-e$y))[1]),]
  #     group <- data_point$Cell_Line_Group
  #     paste0("Cell Line Group: ", group, "\n" , "Time Grown: ", round(e$x, 1), " days","\n", "Copy Number: ", round(e$y, 1), " copies/cell", "\n")
  #   }
  #   
  #   paste0(
  #     xy_str(input$time_plot_click)[1]
  #   )
  # })
  # 
  CN_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$Copy_Number
    Y_title <- "Copy Number (copies/cell"


    p <- plotTime(Data, parameter,Y_title)
    p
  })

  # output$downloadPlot_Copy_Number <- downloadHandler(
  #   filename = paste('Copy_Number_Plot.', "jpeg", sep = ""),
  #   content = function(file) {
  #     ggsave(file, plot = CN_Plot(), width = 20, height = 12, units = "cm", device = "jpeg")
  #   })
  
  ############ 
  
  ### mtDNA Deletions ###
  ############ 
  
  #Parameter list
  output$Deletions_ParameterSelector<-renderUI({
    selectInput('Deletions_Parameter', 'Deletions_Parameters',
                Deletions_Parameters, 
                multiple=FALSE, 
                selectize=TRUE, 
                selected="mtDNA_deletion_Nb_de_délétions") #default value
  })
  
  #get the selected Parameters
  Deletions_Selected_Parameters <-reactive({
    if(is.null(input$Deletions_Parameter) || length(input$Deletions_Parameter)==0)
      return()
    as.vector(input$Deletions_Parameter)
  })
  
  # Insert the right number of plot output objects into the web page
  output$Deletions_Plots <- renderPlotly({
    Data <- Selected_Data()
    parameter <- as.character(sort(Deletions_Selected_Parameters()))
    Y_title <- parameter
    parameter_col <- as.numeric(which(names(Data)==parameter))
    parameter <- as.numeric(Data[,parameter_col])
    #Y_title <- paste0(Y_title, " ")
    # if(input$Cytokine_norm == TRUE) {
    #   parameter <-  parameter / Data$Cells_Counted
    #   Y_title <- paste(Y_title, "/1M cells", sep = "")
    # }
    #print(parameter)
    
    p <- plotTime(Data, parameter,Y_title)
    p
  })
  
  ##############
  
  
  ### cf-DNA ###
  ############   
    
    # Insert cf-DNA Plots
    output$cf_Plots <- renderPlotly({
      Data <- Selected_Data()
      parameter <- Data$cf_mtDNA
      Y_title <- "cf-mtDNA (copies/ml"
      
      if(input$`cf-DNA` == "nDNA") {
        parameter <- Data$cf_nDNA
        Y_title <- "cf-nDNA (copies/ml"
      }
      else if(input$`cf-DNA` == "both") {
        parameter <- Data$cf_mtDNA / Data$cf_nDNA
        Y_title <- "cf-mtDNA/cf-nDNA (copies"
      }
      else if(input$`cf-DNA` == "copy_number") {
        parameter <- Data$cf_mtDNA / Data$Copy_Number
        Y_title <- "cf-mtDNA/copy number (copies"
      }
      
      if(input$cf_norm == TRUE) {
        parameter <-  parameter / Data$Cells_Counted
        Y_title <- paste(Y_title, "/1M cells", sep = "")
      }
      #Y_title <- paste(Y_title, ")", sep = "")
      
      p <- plotTime(Data, parameter,Y_title)
      p
      
    })
    
  # # Allow for hovering/clicking of data points
  # output$cfDNA_DataPoint <- renderText({
  #  
  #   xy_str <- function(e) {
  #     if(is.null(e)) return("Select Data Point for more information...\n")
  #     data_points <- Selected_Data()[which(abs(Selected_Data()$Days_Grown-e$x)%in% sort(abs(Selected_Data()$Days_Grown-e$x))[1:5]),]
  #     data_point <- data_points[which(abs(data_points$cf_mtDNA-e$y) %in% sort(abs(data_points$cf_mtDNA-e$y))[1]),]
  #     group <- data_point$Cell_Line_Group
  #     paste0("Cell Line Group: ", group, "\n" , "Time Grown: ", round(e$x, 1), " days","\n", "cell free mtDNA: ", round(e$y, 1), " copies/ml", "\n")
  #   }
  #   
  #   paste0(
  #     xy_str(input$time_plot_click)[1]
  #   )
  # })
  # 
  cfDNA_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$cf_mtDNA
    Y_title <- "cf-mtDNA (copies/ml"

    if(input$`cf-DNA` == "nDNA") {
      parameter <- Data$cf_nDNA
      Y_title <- "cf-nDNA (copies/ml"
    }
    else if(input$`cf-DNA` == "both") {
      parameter <- Data$cf_mtDNA / Data$cf_nDNA
      Y_title <- "cf-mtDNA/cf-nDNA (copies/ml"
    }

    if(input$cf_norm == TRUE) {
      parameter <-  parameter / Data$Cells_Counted
      Y_title <- paste(Y_title, "/1M cells", sep = "")
    }
    #Y_title <- paste(Y_title, ")", sep = "")

    p <- plotTime(Data, parameter,Y_title)
    p
  })
  # 
  # output$downloadPlot_cfDNA <- downloadHandler(
  #   filename = paste('cell_free_DNA_Plot.', "jpeg", sep = ""),
  #   content = function(file) {
  #     ggsave(file, plot = cfDNA_Plot(), width = 20, height = 12, units = "cm", device = "jpeg")
  # })
  
  ############   
  
  #### IL-6 ###
  #############  
  
  # Insert IL6 Plots
  output$IL6_Plots <- renderPlotly({
    Data <- Selected_Data()
    parameter <- Data$IL6
    Y_title <- "IL6 (pg/ml"
    
    if(input$IL6_norm == TRUE) {
      parameter <-  Data$IL6/Data$Cells_Counted
      Y_title <- paste(Y_title, "/1M cells", sep = "")
    }
    
    #Y_title <- paste(Y_title, ")", sep = "")
    
    
    p <- plotTime(Data, parameter,Y_title)
    p
    
    
  })
  
  # # Allow for hovering/clicking of data points
  # output$IL6_DataPoint <- renderText({
  #   parameter <- Selected_Data()$IL6
  #   xy_str <- function(e) {
  #     if(is.null(e)) return("Select Data Point for more information...\n")
  #     data_points <- Selected_Data()[which(abs(Selected_Data()$Days_Grown-e$x)%in% sort(abs(Selected_Data()$Days_Grown-e$x))[1:5]),]
  #     data_point <- data_points[which(abs(data_points$IL6-e$y) %in% sort(abs(data_points$IL6-e$y))[1]),]
  #     group <- data_point$Cell_Line_Group
  #     paste0("Cell Line Group: ", group, "\n" , "Time Grown: ", round(e$x, 1), " days","\n", "IL6: ", round(e$y, 1), " pg/ml", "\n")
  #   }
  #   
  #   paste0(
  #     xy_str(input$time_plot_click)[1]
  #   )
  # })
  # 
  IL6_Plot <- reactive({
    Data <- Selected_Data()
    parameter <- Data$IL6
    Y_title <- "IL6 (pg/ml"

    if(input$IL6_norm == TRUE) {
      parameter <-  Data$IL6/Data$Cells_Counted
      Y_title <- paste(Y_title, "/1M cells", sep = "")
    }
    #Y_title <- paste(Y_title, ")", sep = "")

    p <- plotTime(Data, parameter,Y_title)
    p
  })
  # 
  # output$downloadPlot_IL6 <- downloadHandler(
  #   filename = paste('cell_free_IL6_Plot.', "jpeg", sep = ""),
  #   content = function(file) {
  #     ggsave(file, plot = IL6_Plot(), width = 20, height = 12, units = "cm", device = "jpeg")
  #   })
  
  ############ 
  
  ### Cytokines ###
  ############ 
  
  #Parameter list
  output$Cytokines_ParameterSelector<-renderUI({
    selectInput('Cytokine_Parameter', 'Cytokine_Parameters',
                Cytokine_Parameters, 
                multiple=FALSE, 
                selectize=TRUE, 
                selected="Lumican_BR13_18_plex") #default value
  })
  
  #get the selected Parameters
  Cytokine_Selected_Parameters <-reactive({
    if(is.null(input$Cytokine_Parameter) || length(input$Cytokine_Parameter)==0)
      return()
    as.vector(input$Cytokine_Parameter)
  })
  
  # Insert the right number of plot output objects into the web page
  output$Cytokine_Plots <- renderPlotly({
    
    Data <- Selected_Data()
    parameter <- as.character(sort(Cytokine_Selected_Parameters()))
    Y_title <- parameter
    parameter_col <- as.numeric(which(names(Data)==parameter))
    parameter <- as.numeric(Data[,parameter_col])
    Y_title <- paste0(Y_title, " pg")
    if(input$Cytokine_norm == TRUE) {
      parameter <-  parameter / Data$Cells_Counted
      Y_title <- paste(Y_title, "/1M cells", sep = "")
    }
    #print(parameter)
    
    
    
    
    p <- plotTime(Data, parameter,Y_title)
    p
  })
  
  
  
  # # Allow for hovering/clicking of data points
  # output$Seahorse_DataPoint <- renderText({
  #   parameter_name <- as.character(sort(Selected_Parameters()))
  #   parameter_col <- as.numeric(which(names(Selected_Data())==parameter_name))
  #   parameter <- as.numeric(Selected_Data()[,parameter_col])
  #   xy_str <- function(e) {
  #     if(is.null(e)) return("Select Data Point for more information...\n")
  #     data_points <- Selected_Data()[which(abs(Selected_Data()$Days_Grown-e$x)%in% sort(abs(Selected_Data()$Days_Grown-e$x))[1:5]),]
  #     data_point <- data_points[which(abs(data_points[,parameter_col]-e$y) %in% sort(abs(data_points[,parameter_col]-e$y))[1]),]
  #     group <- data_point$Cell_Line_Group
  #     paste0("Cell Line Group: ", group, "\n" , "Time Grown: ", round(e$x, 1), " days","\n", parameter_name, ": ", round(e$y, 1), "\n")
  #   }
  #   
  #   paste0(
  #     xy_str(input$time_plot_click)[1]
  #   )
  # })
  # 
  # Seahorse_Plot <- reactive ({
  #   Data <- Selected_Data()
  #   parameter <- as.character(sort(Selected_Parameters()))
  #   Y_title <- parameter
  #   parameter_col <- as.numeric(which(names(Data)==parameter))
  #   parameter <- as.numeric(Data[,parameter_col])
  #   
  #   if(Y_title == "Baseline_ECAR" | Y_title == "Max_ECAR") {
  #     Y_title <- paste(Y_title, "(mpH/")
  #   }
  #   else if(Y_title == "ATPglyc" | Y_title == "ATPox" | Y_title == "ATPproduction" ) {
  #     Y_title <- paste(Y_title, "(pmol ATP/")
  #   }
  #   else {
  #     Y_title <- paste(Y_title, "(pmol O2/")
  #   }
  #   
  #   
  #   if(input$Division_Norm == TRUE) {
  #     parameter <- parameter * 1440 * Data$Days_Per_Passage / Data$Divisions_Per_Passage
  #     Y_title <- paste(Y_title, "division/AUcells)")
  #   }
  #   
  #   else {
  #     Y_title <- paste(Y_title, "min/AU cells)", sep = "")
  #   }
  #   
  #   
  #   p <- plotTime(Data, parameter,Y_title)
  #   p
  # })
  # 
  # output$downloadPlot_Seahorse <- downloadHandler(
  #   filename = paste('Seahorse_Plot.', "jpeg", sep = ""),
  #   content = function(file) {
  #     ggsave(file, plot = Seahorse_Plot(), width = 20, height = 12, units = "cm", device = "jpeg")
  #   })
  
  ##############

  
  ### Seahorse ###
  ############ 
    
    #Parameter list
    output$ParameterSelector<-renderUI({
      selectInput('Parameter', 'Parameters',
                  Seahorse_Parameters, 
                  multiple=FALSE, 
                  selectize=TRUE, 
                  selected="ATPproduction") #default value
    })
    
    #get the selected Parameters
    Selected_Parameters <-reactive({
      if(is.null(input$Parameter) || length(input$Parameter)==0)
        return()
      as.vector(input$Parameter)
    })
    
    # Insert the right number of plot output objects into the web page
    output$Seahorse_Plots <- renderPlotly({

          Data <- Selected_Data()
          parameter <- as.character(sort(Selected_Parameters()))
          Y_title <- paste0(parameter, " (")
          if(grepl("ATP", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "pmolATP/min/20kcells")
          }
          else if(grepl("respiration", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "pmolO2/min/20kcells")
          }
          else if(grepl("Oxidative", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "pmolO2/min/20kcells")
          }
          else if(grepl("Proton", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "pmolO2/min/20kcells")
          }
          else if(grepl("ECAR", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "mpH/min/20kcells")
          }
          else if(grepl("PPR", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "pmolH/min/20kcells")
          }
          else if(grepl("Coupling", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "%")
          }
          else if(grepl("Metabolic", parameter, fixed = F)) {
            Y_title <- paste0(Y_title, "mlO2/min/kg")
          }
          parameter_col <- as.numeric(which(names(Data)==parameter))
          parameter <- as.numeric(Data[,parameter_col])
          
          
          if(input$Run_norm == TRUE) {
            #LS_Data$Seahorse_Run_Number <- as.numeric(LS_Data$Seahorse_Run_Number)
            parameter_name <- as.character(sort(Selected_Parameters()))
            parameter_col <- as.numeric(which(names(Selected_Data())==parameter_name))
            
            runs <- unique(LS_Data$Seahorse_Run_Number[!is.na(LS_Data$Seahorse_Run_Number)])
            mean_run <- matrix(ncol=2, nrow = length(runs))
            for(i in 1:length(runs)) {
              #i <- 1
              run <- as.numeric(runs[i])
              mean_run[i,1] <- run
              
              run_data <- LS_Data[LS_Data$Seahorse_Run_Number %in% run,]
              mean_run[i,2] <- mean(run_data[,parameter_col], na.rm = TRUE)
            }
            #percent_mean_run <- mean_run
            correction_factor <- mean_run[,2] / mean(mean_run[,2], na.rm = TRUE)
            mean_run <- cbind(mean_run, correction_factor)
            #percent_mean_run <- cbind(mean[,1], percent_mean_run)
            correction_col <- matrix(ncol=1, nrow = nrow(Data))
            for(i in 1:nrow(Data)) {
              run <- Data[i,]$Seahorse_Run_Number
              if(is.na(run)) {
                correction_col[i] <- NA
              }
              else {
                correction_col[i] <- mean_run[mean_run[,1] %in% run,3]
              }
              
            }
            
            parameter <- parameter / correction_col
            
                
          }
          
          if(input$Division_Norm == TRUE) {
            div_per_passage <- Data$Divisions_Per_Passage
            div_per_passage[round(div_per_passage,1) == 0.0] <- NA
            parameter <- parameter * 1440 * Data$Days_Per_Passage / div_per_passage
            Y_title <- paste(Y_title, "/division", sep = "")
          }

          p <- plotTime(Data, parameter,Y_title)
          p
    })
  ##############
  
  ### Correlation Graphs ###
  ############ 
  
  #Parameter list 2D Correlation Graph
  output$ParameterSelector2D_x <-renderUI({
    selectInput('Cor_Parameter_x', 'X Axis Parameter',
                colnames(LS_Data), 
                multiple=FALSE, 
                selectize=TRUE, 
                selected="ATPglyc") #default value
  })
  
  output$ParameterSelector2D_y <-renderUI({
    selectInput('Cor_Parameter_y', 'Y Axis Parameter',
                colnames(LS_Data), 
                multiple=FALSE, 
                selectize=TRUE, 
                selected="ATPox") #default value
  })
  
  #get the selected 2D X parameter
  Selected_Correlation_Parameter_x <-reactive({
    if(is.null(input$Cor_Parameter_x) || length(input$Cor_Parameter_x)==0)
      return()
    as.vector(input$Cor_Parameter_x)
  })
  
  #get the selected 2D y parameter
  Selected_Correlation_Parameter_y <-reactive({
    if(is.null(input$Cor_Parameter_y) || length(input$Cor_Parameter_y)==0)
      return()
    as.vector(input$Cor_Parameter_y)
  })
  
  observe({

    if(input$tabs == "Correlations") {
      updateSliderInput(session, "poly", value = 1)
      updateCheckboxInput(session, "fit", value = FALSE)  
      #updateCheckboxInput(session, "fit",  value = FALSE)

    }
    else {
      updateSliderInput(session, "poly", value = 3)
      updateCheckboxInput(session, "fit", value = TRUE)  
      #updateCheckboxInput(session, "fit",  value = TRUE)
    }
  })
  
    # Geneate 2D Correlation Graphs
    output$Correlation_2D_Plots <- renderPlot({

        Data <- Selected_Data()
        
        parameter_x <- Selected_Correlation_Parameter_x()
        parameter_y <- Selected_Correlation_Parameter_y()
        
        X_title <- parameter_x
        Y_title <- parameter_y
        
        parameter_x_col <- as.numeric(which(names(Data)==parameter_x))
        parameter_y_col <- as.numeric(which(names(Data)==parameter_y))
        
        parameter_x_data <- as.numeric(Data[,parameter_x_col])
        parameter_y_data <- as.numeric(Data[,parameter_y_col])
        
        if(parameter_x == "IL6" | parameter_x == "cf_mtDNA" | parameter_x == "cf_nDNA") {
          parameter_x_data <- parameter_x_data / Data$Cells_Counted
        }
        
        if(parameter_y == "IL6" | parameter_y == "cf_mtDNA" | parameter_x == "cf_nDNA") {
          parameter_y_data <- parameter_y_data / Data$Cells_Counted
        }
        
        if(input$Division_Norm == TRUE) {
          if(str_detect(Y_title, "clock")) {
            parameter_y_data <- parameter_y_data * 1440 * Data$Days_Per_Passage
          }
          
          if(str_detect(X_title, "clock")) {
            parameter_x_data <- parameter_x_data * 1440 * Data$Days_Per_Passage 
          }
          
          div_per_passage <- Data$Divisions_Per_Passage
          div_per_passage[round(div_per_passage,1) == 0.0] <- NA 
          
          parameter_x_data <- parameter_x_data / div_per_passage
          X_title <- paste(X_title, "/division")
          
          parameter_y_data <- parameter_y_data / div_per_passage
          Y_title <- paste(Y_title, "/division")
        }
        
        if(input$Norm_Cell_Volume == TRUE) {
          Cell_Volume <- (Data$Cell_Size /2)^3 * pi * (4/3)
          parameter_x_data <- parameter_x_data / Cell_Volume
          X_title <- paste(X_title, "/cell volume", sep ="")
          
          parameter_y_data <- parameter_y_data / Cell_Volume
          Y_title <- paste(Y_title, "/cell volume", sep = "")
        }
        
        
        if(input$scale == "log") {
          parameter_x_data <- log10(parameter_x_data)
          parameter_x_data[is.infinite(parameter_x_data)] <- NA 
          #parameter_x_data[is.nan(parameter_x_data)] <- NA
          
          parameter_y_data <- log10(parameter_y_data)
          parameter_y_data[is.infinite(parameter_y_data)] <- NA 
          #parameter_y_data[is.nan(parameter_y_data)] <- NA
          
          X_title <- paste("Log10 ", X_title, sep = "")
          Y_title <- paste("Log10 ", Y_title, sep = "")
        }
        
        outliers_data <- NULL
        if(input$outliers == "yes") {
          outliers_x <- boxplot(parameter_x_data, plot=FALSE)$out
          #Data[(parameter_x_data %in% outliers_x),] <- NA
          parameter_x_data[(parameter_x_data %in% outliers_x)] <- NA
          outliers_data_x <- Data[(parameter_x_data %in% outliers_x),]
          
          outliers_y <- boxplot(parameter_y_data, plot=FALSE)$out
          #Data[(parameter_y_data %in% outliers_y),] <- NA
          parameter_y_data[(parameter_y_data %in% outliers_y)] <- NA
          outliers_data_y <- Data[(parameter_y_data %in% outliers_y),]
          
          outliers_data <- merge(outliers_data_x, outliers_data_y, by = "Unique_Variable_Name")
          
        }
        
        annotation <- ""
        if(!is.null(outliers_data)) {
          n <- length(outliers_data$Days_Grown)
          tag <- rep("days", each=n)
          outliers_days_grown <- paste(round(as.numeric(outliers_data$Days_Grown),0), tag, sep = " ")
          concate_cols <- paste(outliers_data$Cell_Line_Group, outliers_days_grown, sep = ": ")
          outlier_string <- paste(concate_cols , collapse = ", ")
          annotation <- paste("Extreme Outlier(s) Removed: ", outlier_string)
        }
        
        treatments <- sort(Selected_Treatments())
        if(is.null(treatments)) {
          colors <- 'dimgray'
          colors_2 <- 'dimgray'
          shapes <- 21
        }
        else {
          colors <- getColors(treatments)
          colors_2 <- getColors_2(sort(Selected_Cell_Lines()))
          shapes <- getShapes(sort(Selected_Cell_Lines()))
        }
        scale <- input$scale
        SE <- input$se
        DoF <- input$poly
        formula <- y ~ poly(x, DoF, raw = TRUE)
        
        x_data <- parameter_x_data
        y_data <- parameter_y_data
        X_title2 <- X_title
        Y_title2 <- Y_title
        
        time_points <- Data$Days_Grown 
        if(input$Fun_Regression == TRUE) {
          poly_x <- as.numeric(input$poly_x)
          poly_y <- as.numeric(input$poly_y)
          
         #  if(input$Optomize_Fit == TRUE) {
         #    Optomize_Polynomials <- function(data, time_points, data2, poly2) {
         #      best_DoF <- 1
         #      best_r2 <- 0
         #      
         #      # data <- data.frame(data = data)
         #      #tp <- time_points #data.frame(time_points = time_points)
         #      # data2 <- data.frame(data2 = data2)
         #      # poly2 <- poly2
         #      
         #      model2 <- lm(data2 ~ poly(time_points, poly2, raw = TRUE))
         #      tp <- data.frame(tp = time_points)
         #      model2_data <- as.numeric(predict(model2, newdata = tp))
         #      
         #      for(i in 1:9) {
         #        DoF <- 11 - i
         #        model <- lm(data ~ poly(time_points, DoF, raw = TRUE))
         #        model_data <- as.numeric(predict(model, newdata = tp))
         #        
         #        final_model <- lm(model2_data ~ poly(model_data, 3, raw = TRUE))
         #        r2 <- as.numeric(summary(final_model)$r.squared)
         #        
         #        if(r2 > best_r2) {
         #          best_r2 <- r2
         #          best_DoF <- DoF
         #        }
         #      }
         #      return(best_DoF)
         #   }
         #    
         #    poly_x <- Optomize_Polynomials(parameter_x_data, time_points, parameter_y_data, poly_y)
         #    updateSliderInput(session, "poly_x", value = poly_x)
         #    poly_y <- Optomize_Polynomials(parameter_y_data, time_points, parameter_x_data, poly_x)
         #    updateSliderInput(session, "poly_y", value = poly_y)
         # }
          
          startCount <- 1
          final_data_x <- NULL
          final_data_y <- NULL
          final_error_x <- NULL
          final_error_y <- NULL
          x_r2 <- c()
          y_r2 <- c()
          for (i in 1:length(unique(Data$Cell_line_group_new))) {
            
            group <- unique(Data$Cell_line_group_new)[i]
            
            endCount <-  startCount + nrow(Data[Data$Cell_line_group_new %in% group,]) - 1
            
            time_points <- Data[startCount:endCount,]$Days_Grown
            group_data_x <- parameter_x_data[startCount:endCount]
            group_data_y <- parameter_y_data[startCount:endCount]
              
            x_model <- gam(group_data_x ~ s(time_points,bs="cr", k = poly_x))
            y_model <- gam(group_data_y ~ s(time_points,bs="cr", k = poly_y))
            
            
            x_r2 <- append(x_r2, summary(x_model)$r.sq)
            y_r2 <- append(y_r2, summary(y_model)$r.sq)
            
            
            #time_points <- seq(min(Data$Days_Grown), max(Data$Days_Grown), length = nrow(Data))
            time_points <- data.frame(time_points = time_points)
            
            predicted_x <- predict(x_model, newdata = time_points, se.fit=TRUE)
            predicted_y <- predict(y_model, newdata = time_points, se.fit=TRUE)
            
            predicted_x <- as.numeric(unlist(predicted_x[[1]]))
            predicted_y <- as.numeric(unlist(predicted_y[[1]]))
            final_data_x <- c(final_data_x, predicted_x)
            final_data_y <- c(final_data_y, predicted_y)
            
            predicted_error_x <- as.numeric(unlist(predicted_x[[2]]))
            predicted_error_y <- as.numeric(unlist(predicted_y[[2]]))
            final_error_x <- c(final_error_x, predicted_error_x)
            final_error_y <- c(final_error_y, predicted_error_y)
            
            startCount <- endCount
          }
          
          X_title <- paste("Fun(", X_title, ", ", input$poly_x, ")", sep = "")
          Y_title <- paste("Fun(", Y_title, ", ", input$poly_y, ")", sep = "")
          
          parameter_x_data <- final_data_x
          parameter_y_data <- final_data_y
          
          xMin <- final_data_x-final_error_x
          xMax <-  final_data_x +final_error_x
          
          yMin <-  final_data_y-final_error_y
          yMax <-  final_data_y+final_error_y
      
        }
        
        p <- ggplot(Data, aes_string(x=parameter_x_data, y = parameter_y_data, color ="Treatment")) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors_2) +
          scale_shape_manual(values=shapes) +
          scale_y_continuous(name = Y_title) +
          scale_x_continuous(name =  X_title) +
          theme_classic() +
          annotate("text", label = annotation, x = 2, y = 0, size = 3, hjust = 0) +
          coord_cartesian(clip = "off")+
          theme(text = element_text(size = 14),
                legend.position="none",
                axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                plot.margin = margin(t = 6, r = 6, b = 6, l = 6, unit = "pt"),
                plot.title = element_text(size = 32, hjust = 0.05, vjust = -.1))
        
        groups <- unique(Data$Cell_line_group_new)
         
        if (input$Fun_Regression == TRUE) {
          count <- length(groups)
          #groups <- paste(unique(Data$Cell_Line_Group), rep("ParameterX_r2:", count), round(x_r2,2), rep("ParameterY_r2:", count), round(y_r2, 2),  sep = "~~")
          #ypos <- max(parameter_y_data, na.rm=T)
          #xpos <-  max(parameter_x_data, na.rm=T)
          #p = p + annotate("text", label = paste("Parameter X r^2: ", toString(round(x_r2,2)), "\n","Parameter Y r^2: ", toString(round(y_r2, 2)), sep = ""), size = 4, x = xpos, y= ypos)
          p = p + geom_line(size = 0.5 , aes(shape = Cell_line_new, color = Treatment, fill = Cell_line_new))
          p = p + geom_text(aes(label=paste(round(Days_Grown, digits = 0), "days", sep = " ")),hjust=-0.2, vjust=-0.2, size = 2, alpha = 0.7)
          if(input$se == TRUE) {
            p = p + geom_errorbar(aes(ymin=yMin, ymax=yMax), width = 10)
            p = p + geom_errorbarh(aes(xmin=xMin, xmax=xMax), height=10)
          }
        }
        else {
          p = p + geom_point(size = 2 , stroke = 1, aes(shape = Cell_line_new, color = Treatment, fill = Cell_line_new))
          p = p + geom_text(aes(label=paste(round(Days_Grown, digits = 0), "days", sep = " ")),hjust=-0.2, vjust=-0.2, size = 2, alpha = 0.7)
          
        }
        
        if(input$fit == TRUE) { 
          p = p + geom_smooth(method = "lm", formula = formula, size = 1.5, aes(group = Cell_line_group_new), alpha = 0.05, se = SE)
        }
        
        
        # Plot 2nd graph overlaying the Time Plots of the Two Variables
        p2 <- plotTime(Data, y_data, Y_title2)
      
        # Transform Right Y axis data to allow it to be comparable
        SF <- mean(x_data, na.rm = TRUE) / mean(y_data, na.rm = TRUE)
        adjusted_data <- x_data / SF
        p2 = p2 + geom_point(size = 2, stroke = 1, aes(y = adjusted_data, shape = Cell_line_new, color = Treatment, fill = Cell_line_new))
        p2 = p2 + scale_y_continuous(name = Y_title2, sec.axis = sec_axis(~.*SF, name = X_title2))
        
        #p2 = p2 + theme(axis.title.y = element_text(colour = "red"))
        
        # Correct Smooth Function of left Y axis
        formulaY <- y ~ poly(x, input$poly_y, raw = TRUE)
        p2 = p2 + geom_smooth(method = "lm", formula = formulaY, size =0.8, aes(group = Cell_line_group_new), alpha = 0.05, se = SE)
        
        
        # Add smooth function to right Y axis
        formulaX <- y ~ poly(x, input$poly_x, raw = TRUE)
        p2 = p2 + geom_smooth(method = "lm", formula = formulaX, size =0.8, aes(y = adjusted_data, group = Cell_line_group_new), alpha = 0.05, se = SE)
        
        # Convert 2 Graphs to Plotly and arrange in grid
        #p1 <- ggplotly(p)
        #p2 <- ggplotly(p2)
        
        # Add second Y axis to Plotly second graph
        #p2 = p2 %>% layout(yaxis2 = list(overlaying = "y", side = "right"))
        
        # Add equation info to first and second graph
        if (input$fit_text == TRUE) {
          count <- length(unique(Data$Cell_line_group_new))
          p = p + stat_poly_eq(formula = formula, aes(group = Cell_line_group_new, label =  paste(groups, stat(eq.label), stat(rr.label), sep = "~~")),
                                rr.digits = 2, coef.digits = 2,parse = TRUE,  label.x = "left", label.y = "bottom", show.legend = TRUE)
          p2 = p2 + stat_poly_eq(formula = formulaY, aes(group = Cell_line_group_new, label =  paste(groups, stat(rr.label), sep = "~~")),
                               rr.digits = 2, coef.digits = 2,parse = TRUE,  label.x = "left", label.y = "top", show.legend = TRUE)
          
          p2 = p2 + stat_poly_eq(formula = formulaX, aes(y = adjusted_data, group = Cell_line_group_new, label =  paste(groups, stat(rr.label), sep = "~~")),
                                 rr.digits = 2, coef.digits = 2,parse = TRUE,  label.x = "right", label.y = "bottom", show.legend = TRUE)
          
        }
        
        p_list <- list(p,p2)
        grid.arrange(grobs = p_list, widths = c(10,10), ncol=2)
        #subplot(p1, p2, titleX = TRUE, titleY = TRUE)

  })
    
    
    
    #Parameter list 3D Correlation Graph
    output$ParameterSelector3D_x <-renderUI({
      selectInput('Cor_Parameter3D_x', 'X Axis Parameter',
                  colnames(LS_Data), 
                  multiple=FALSE, 
                  selectize=TRUE, 
                  selected="Population_Doublings_DI") #default value
    })
    
    output$ParameterSelector3D_y <-renderUI({
      selectInput('Cor_Parameter3D_y', 'Y Axis Parameter',
                  colnames(LS_Data), 
                  multiple=FALSE, 
                  selectize=TRUE, 
                  selected="cf_mtDNA") #default value
    })
    
    output$ParameterSelector3D_z <-renderUI({
      selectInput('Cor_Parameter3D_z', 'Z Axis Parameter',
                  colnames(LS_Data), 
                  multiple=FALSE, 
                  selectize=TRUE, 
                  selected="IL6") #default value
    })
    
    #get the selected 3D X parameter
    Selected_Correlation_Parameter3D_x <-reactive({
      if(is.null(input$Cor_Parameter3D_x) || length(input$Cor_Parameter3D_x)==0)
        return()
      as.vector(input$Cor_Parameter3D_x)
    })
    
    #get the selected 3D y parameter
    Selected_Correlation_Parameter3D_y <-reactive({
      if(is.null(input$Cor_Parameter3D_y) || length(input$Cor_Parameter3D_y)==0)
        return()
      as.vector(input$Cor_Parameter3D_y)
    })
    
    #get the selected 3D z parameter
    Selected_Correlation_Parameter3D_z <-reactive({
      if(is.null(input$Cor_Parameter3D_z) || length(input$Cor_Parameter3D_z)==0)
        return()
      as.vector(input$Cor_Parameter3D_z)
    })
    
    # Generate 3D Correlation Graph
    output$Correlation_3D_Plots <- renderRglwidget({
      
      Data <- Selected_Data()
      
      parameter_x <- Selected_Correlation_Parameter3D_x()
      parameter_y <- Selected_Correlation_Parameter3D_y()
      parameter_z <- Selected_Correlation_Parameter3D_z()
      
      
      X_title <- parameter_x
      Y_title <- parameter_y
      Z_title <- parameter_z
      
      parameter_x_col <- as.numeric(which(names(Data)==parameter_x))
      parameter_y_col <- as.numeric(which(names(Data)==parameter_y))
      parameter_z_col <- as.numeric(which(names(Data)==parameter_z))
      
      parameter_x_data <- as.numeric(Data[,parameter_x_col])
      parameter_y_data <- as.numeric(Data[,parameter_y_col])
      parameter_z_data <- as.numeric(Data[,parameter_z_col])
      
      if(parameter_x == "IL6" | parameter_x == "cf_mtDNA" | parameter_x == "cf_nDNA") {
        parameter_x_data <- parameter_x_data / Data$Cells_Counted
      }
      
      if(parameter_y == "IL6" | parameter_y == "cf_mtDNA" | parameter_x == "cf_nDNA") {
        parameter_y_data <- parameter_y_data / Data$Cells_Counted
      }
      
      if(parameter_z == "IL6" | parameter_z == "cf_mtDNA" | parameter_z == "cf_nDNA") {
        parameter_z_data <- parameter_z_data / Data$Cells_Counted
      }
      
      
      if(input$scale == "log") {
        parameter_x_data <- log10(parameter_x_data)
        parameter_y_data <- log10(parameter_y_data)
        parameter_z_data <- log10(parameter_z_data)
        
        X_title <- paste("Log10 ", X_title, sep = "")
        Y_title <- paste("Log10 ", Y_title, sep = "")
        Z_title <- paste("Log10 ", Z_title, sep = "")
      }
      
      outliers_x <- NULL
      outliers_y <- NULL
      outliers_z <- NULL
      if(input$outliers == "yes") {
        outliers_x <- boxplot(parameter_x_data, plot=FALSE)$out
        Data <- Data[!(parameter_x_data %in% outliers_x),]
        
        outliers_y <- boxplot(parameter_y_data, plot=FALSE)$out
        Data <- Data[!(parameter_y_data %in% outliers_y),]
        
        outliers_z <- boxplot(parameter_z_data, plot=FALSE)$out
        Data <- Data[!(parameter_z_data %in% outliers_z),]
      }
      
      annotation <- ""
      if(!is.null(outliers_x) | !is.null(outliers_y)) {
        outlier_string_x <- paste(outliers_x, collapse=", ")
        outlier_string_y <- paste(outliers_y, collapse=", ")
        outlier_string_z <- paste(outliers_z, collapse=", ")
        annotation <- paste("Extreme Outlier(s) Removed: ", outlier_string_x, outlier_string_y, outlier_string_z)
      }
      
      treatments <- sort(Selected_Treatments())
      if(is.null(treatments)) {
        colors <- 'dimgray'
        colors_2 <- 'dimgray'
        shapes <- 21
        
      }
      else {
        colors <- getColors(treatments)
        colors_2 <- getColors_2(sort(Selected_Cell_Lines()))
        shapes <- getShapes(sort(Selected_Cell_Lines()))
      }
      
      
      rgl.open(useNULL=T)
      scatter3d(x = parameter_x_data, y = parameter_y_data, z = parameter_z_data,
                data=Data, 
                id= list(labels = paste(round(Data$Days_Grown, 0), " days"), n=nrow(Data)),
                id.n=nrow(Data),
                groups = as.factor(Data$Cell_Line_Group),
                #fit = "smooth",
                ellipsoid = input$ellipsoid,
                ellipsoid.alpha=0.01,
                #labels = as.character(Data$Days_Grown),
                #point.col = colors,
                #grid = TRUE,
                #surface.col = colors,
                xlab = X_title,
                ylab = Y_title,
                zlab = Z_title,
                axis.ticks=TRUE,
                #axis.col = c("blue", "red", "purple"),
                surface = FALSE
                #fill = TRUE,
                #panel.cols = "yellow"
                #col.grid = "darkblue"
                )
      
      # Identify3d(x = parameter_x, y = parameter_y, z = parameter_z, groups = as.factor(Data$Cell_Line_Group), labels = as.character(Data$Days_Grown))
      rglwidget()
      
      
    })
  
  ############
  
  
  # Generate an HTML table view of the data
  output$DataMatrix <- renderTable({
    metadata <- read.csv(paste0("downloadable_data/Cellular_lifespan_study_",Selected_Dataset(),".csv"))
    output_data <- Selected_Data()
    output_data <- metadata[metadata$Sample %in% output_data$Sample,]
    if(Selected_Dataset() != "all") {
      output_data <- output_data[!is.na(output_data[,length(output_data)]),]
    }
    output_data
  })
  


})
