# Load package and the start function
# library(shiny)
library(bslib)
library(shinyjs)  # Enable dynamic JS control (like changing button color)

## Our App

### Loading Data Into The App
# library(here)
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(patchwork)
library(cowplot)
library(magick)
library(pdftools)
library(grid)
library(scales)
# Automatically set working directory to script location (for RStudio only)
# if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
#   this_path <- dirname(rstudioapi::getSourceEditorContext()$path)
#   setwd(this_path)
#   print(paste("Working directory set to:", this_path))
# }

# Now load files relative to this directory
source(file = "helpers.R")


### Setup `ui.R` #####################################################################

ui <- page_navbar(
  title = tags$b("Mobile Element Insertion (MEI) Visualization"),
  id = "main_navbar", 
  
  #### add header
  header = tagList(
    # Activate shinyjs
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .indented-checkbox .form-check {
          margin-left: 10px;
        }
        .indented-checkbox .form-check-label {
        padding-left: 3px;
        }
        hr {
        margin-top: 5px !important;
        margin-bottom: 5px !important;
        }
      "))
    )
  ),
  
  #### set theme
  theme = bs_theme(
    preset = "bootstrap",
    "navbar-padding-bottom" = "20px"
  ) %>% bs_add_rules(
    list(
      ".navbar-brand { margin-left: 5px; margin-right: 100px; }",
      ".navbar-nav { position: relative; left: 99px; }", # move right
      ".nav-link[data-value='overall_dnam_panel'] { margin-left: 18px !important; }"
    )
  ),
  
  #### add sidebar
  sidebar = sidebar(
    selectizeInput(
      inputId = "INS_id",
      # label = "Select an INS (MEI):",  # ■ ● » 🧬 
      label = strong("📍 Select an INS (MEI):"), 
      choices = INS_IDs,
      selected = INS_IDs[1],
      multiple = FALSE
    ),
    
    # br(),
    tags$p("Use the checkbox below to display DNA methylation (DNAm) levels 
           for the selected MEI across sample groups:",
           style = "font-size: 16px; font-weight: normal; color: black; 
           background-color: #f9f9f9; padding: 10px; border-left: 4px solid #64748b; 
           margin-top: 10px; border-radius: 4px;margin-top: 10px;"),
    
    div(
    checkboxInput(
      inputId = "one_emi_boxplot",
      label = "Show DNAm boxplot for the selected MEI",
      value = FALSE),
    
    conditionalPanel(
      condition = "input.one_emi_boxplot == true",
      div(
        class = "indented-checkbox", # #f9f9f9 #f1f5f9
        style = "margin-left: 5px; padding: 10px; background-color: #f9f9f9; border-radius: 5px;",
      
        # Compare group selection
        checkboxGroupInput(
          inputId = "compare_group",
          label = "Compare groups:",
          choices = c("B", "PMM"), inline = TRUE,
          selected = c("B", "PMM")),
        
        # Compare region selection
        checkboxGroupInput(
          inputId = "compare_location",
          label = "Compare regions:",
          choices = c("flank_up", "INS", "flank_down"), inline = TRUE,
          selected = c("flank_up", "INS", "flank_down"))
    ))),
    
    hr(), # Add a horizontal rule #■ ●
    tags$p("📍 Overview of all MEIs:", 
           style = "font-weight: bold; color: black;"),
    tags$p("Click the button below to switch between 
           individual and overall DNAm views (loads in ~10 seconds):", 
           style = "font-size: 16px; font-weight: normal; color: black; 
           background-color: #f9f9f9; padding: 10px; border-left: 4px solid #64748b; 
           margin-top: 10px; border-radius: 4px;margin-top: 10px;"),
    actionButton(
      inputId = "all_emi_boxplot",
      label = "View overall DNAm summary", 
      # class = "btn-warning",
      class = "btn-primary"
    )
  ),
  nav_panel("Selected MEI", value = "selected_MEI_dnam_panel",
            layout_columns(
              col_widths = c(7, 5), 
              fillable = TRUE,
              plotOutput("plots", width = "100%", height = "auto"),
            
            conditionalPanel(
              condition = "input.one_emi_boxplot == true",
              plotOutput("boxplot_ui", width = "500px", height = "690px")
              ))),
  nav_panel("Overall DNAm (all MEIs)", value = "overall_dnam_panel",
            layout_columns(
              col_widths = c(5, 7), 
              fillable = TRUE,
              plotOutput("overall_boxplot_compare1_ui", height = "100%"),
              plotOutput("overall_boxplot_compare2_ui", height = "100%")))
)


## Make The Plot Interactive #####################################################################

server <- function(input, output, session) {
  
  # change page tab and show plot when click the button and 
  # show plot when click 2nd tab
  # A reactive flag to control when to show the all EMI boxplot
  show_overall <- reactiveVal(FALSE)
  
  # Button click toggles tab and flag
  observeEvent(input$all_emi_boxplot, {
    if (!show_overall()) {
      bslib::nav_select("main_navbar", selected = "overall_dnam_panel")
      show_overall(TRUE)
      updateActionButton(session, "all_emi_boxplot", 
                         label = "Return to selected INS view")
      shinyjs::removeClass("all_emi_boxplot", "btn-primary")
      shinyjs::addClass("all_emi_boxplot", "btn-warning")
    } else {
      bslib::nav_select("main_navbar", selected = "selected_MEI_dnam_panel")
      show_overall(FALSE)
      updateActionButton(session, "all_emi_boxplot", 
                         label = "View overall DNAm summary")
      shinyjs::removeClass("all_emi_boxplot", "btn-warning")
      shinyjs::addClass("all_emi_boxplot", "btn-primary")
    }
  })
  
  # Tab click: sync state and button style
  observeEvent(input$main_navbar, {
    if (input$main_navbar == "overall_dnam_panel") {
      show_overall(TRUE)
      updateActionButton(session, "all_emi_boxplot",
                         label = "Return to selected INS view")
      shinyjs::removeClass("all_emi_boxplot", "btn-primary")
      shinyjs::addClass("all_emi_boxplot", "btn-warning")
    } else if (input$main_navbar == "selected_MEI_dnam_panel") {
      show_overall(FALSE)
      updateActionButton(session, "all_emi_boxplot",
                         label = "View overall DNAm summary")
      shinyjs::removeClass("all_emi_boxplot", "btn-warning")
      shinyjs::addClass("all_emi_boxplot", "btn-primary")
    }
  })
  
  
  # 🟡 Plot 1: The main plot, heatmaps for selected INS
  output$plots <- renderPlot({
    ############## get data for one MEI  #######################################
    
    one_mei_flank_eg <- mei_ins_flank_DNAm_df %>%
      dplyr::filter(INS_ID == input$INS_id) %>%
      dplyr::mutate(
        flank_type = factor(flank_type, levels = c("flank_up", "flank_down")),
        cpg_aligned = case_when(
          flank_type == "flank_up" ~ -1 * (cpg_order_within_flank+1),
          flank_type == "flank_down" ~ cpg_order_within_flank+1)
      ) %>% 
      tidyr::complete(sample, flank_type, cpg_aligned,
                      fill = list(DNAm = NA)) %>%
      # -1 used to mark missing values (NA)
      mutate(DNAm_dummy = ifelse(is.na(DNAm), -1, DNAm),
             sample = factor(sample, levels = rev(sample_order))) 
    
    
    ############## Plot p1 -- MEI flank DNAm ###################################
    
    p1 <- ggplot(one_mei_flank_eg, aes(x = cpg_aligned, y = sample)) +
      # Layer 1: Tiles for dummy values (originally NA) → grey tiles
      geom_tile(
        data = dplyr::filter(one_mei_flank_eg, 
                             DNAm_dummy == -1,
                             (flank_type == "flank_up" & cpg_aligned <= 0) |
                               (flank_type == "flank_down" & cpg_aligned >= 0)),
        fill = "grey85", color = "grey52"
      ) +
      # Layer 2: Tiles for non-NA methylation values → blue gradient
      geom_tile(
        data = dplyr::filter(one_mei_flank_eg, DNAm_dummy != -1),
        aes(fill = DNAm_dummy),
        # color = "#F0FFFF"
        ) +
      # Set color gradient (limits exclude dummy values)
      scale_fill_gradient(
        low = "white",
        high = "navy", #"#1A1A2E" "#4169E1"   
        name = "DNAm",
        limits = c(0, 1)  # Only show color gradient for real values between 0 and 1
      ) +
      facet_wrap(~flank_type, scales = "free_x", labeller = as_labeller(c(
        "flank_up" = "Upstream",
        "flank_down" = "Downstream")), drop = FALSE) +
      scale_x_continuous(breaks = function(x) {
        breaks <- pretty(x)
        breaks[breaks == floor(breaks)]  # keep only integer breaks
      }) +
      scale_y_discrete(limits = rev(sample_order))+
      labs(x = "", y = "")+
      theme_pubr(base_size = 20)+
      theme(axis.title.x = element_blank(),
            panel.background = element_rect(fill = "grey95"),
            strip.background = element_rect(fill = "#F0FFFF", color = "black")) 
    # "#FAF0E6" "#FAEBD7" "#ACE1AF" "#98FB98" "#A8E4A0" "#DBECD4"
    
    
    ############## Plot p2 -- ppt image only ###################################
    
    # Read the PDF image using magick
    ppt_image <- image_read_pdf(file.path("test_dat/ins_symbol.pdf"), density = 600)
    ppt_image_grob <- grid::rasterGrob(ppt_image, interpolate = TRUE)
    
    ppt_plot <- ggplot() +
      annotation_custom(ppt_image_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      theme_void() +
      coord_cartesian(expand = FALSE)
    
    
    ############## Plot p3 -- MEI DNAm  ########################################
    
    one_mei_eg <- mei_ins_only_DNAm_df %>%
      dplyr::filter(INS_ID == input$INS_id) %>%
      dplyr::mutate(cpg_aligned = cpg_order_within_ins+1) %>% 
      tidyr::complete(sample, cpg_aligned,
                      fill = list(DNAm = NA)) %>%
      # -1 used to mark missing values (NA)
      mutate(DNAm_dummy = ifelse(is.na(DNAm), -1, DNAm),
             sample = factor(sample, levels = rev(sample_order))) 
     
    # Plot in layers
    p3 <- ggplot() +
      # Layer 1: Tiles for dummy values (originally NA) → grey tiles
      geom_tile(
        data = dplyr::filter(one_mei_eg, DNAm_dummy == -1),
        aes(x = cpg_aligned, y = sample),
        fill = "grey90", color = "grey25"   #"#CCC5B9"#"#D4CFC9" "#B2BEB5"
      ) +
      # Layer 2: Tiles for non-NA methylation values → blue gradient
      geom_tile(
        data = dplyr::filter(one_mei_eg, DNAm_dummy != -1),
        aes(x = cpg_aligned, y = sample, fill = DNAm_dummy),
        # color = "#F0FFFF"
      ) +
      # Set color gradient (limits exclude dummy values)
      scale_fill_gradient(
        low = "white",
        high = "navy", #"#1A1A2E" "#4169E1"
        name = "DNAm",
        limits = c(0, 1)  # Only show color gradient for real values between 0 and 1
      ) +
      scale_x_continuous(breaks = function(x) {
        breaks <- pretty(x)
        breaks[breaks == floor(breaks)]  # keep only integer breaks
      }) +
      scale_y_discrete(limits = rev(sample_order))+
      labs(x = "CpGs ordered by relative position \n (NAs shown in grey)", y = "")+
      theme_pubr(base_size = 20)+
      theme(panel.background = element_rect(fill = "grey95", size = 1))
    
    
    ############## Plot p1+p2+p3 and save ######################################
    # Combine all three plots using patchwork or cowplot
    # Center ppt_plot and p3 using spacers
    ppt_row <- plot_spacer() + ppt_plot + plot_spacer()
    p3_row   <- plot_spacer() + p3 + plot_spacer()
    p3_row  <- p3_row  + plot_layout(widths = c(0.1, 0.8, 0.1))
    
    # get MEI group
    mei_group <- unique(one_mei_eg$mei_info)
    mei_group <- str_remove(mei_group, "rest_")
    mei_group <- str_remove(mei_group, "_MEI..+")
    
    # Combine
    combined_plot <- (p1 / ppt_plot / p3_row) +
      plot_layout(heights = c(1, 0.2, 1), guides = "collect") +
      plot_annotation(
        title = paste0(input$INS_id, ", ", 
                       unique(one_mei_eg$coord), ", ",
                       mei_group)) &
      theme(
        # Plot title
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12.5, colour = "black"),
        strip.text = element_text(colour = "black"),
        
        # Legend title
        legend.title = element_text(size = 20, face = "bold"),
        
        # Legend labels (tick values)
        legend.text = element_text(size = 16, face = "bold", colour = "black"),
        
        # Optional: adjust spacing between items
        legend.key.size = unit(1, "cm"),  # size of boxes
        legend.margin = margin(5, 5, 5, 5),
        
        # Optional: place inside or below plot
        legend.position = "right"
      )
    
    # show the final plot
    print(combined_plot)
  }, res = 100)
  
  
  # 🔵 Plot 2: Boxplot for single selected INS (if checkbox is checked)
  output$boxplot_ui <- renderPlot({
    req(input$one_emi_boxplot, input$INS_id, input$compare_group, input$compare_location)
    
    ## get selected mei df
    one_mei_flank_eg <- mei_ins_flank_DNAm_df %>%
      dplyr::filter(INS_ID == input$INS_id) %>%
      dplyr::mutate(location = flank_type) %>% 
      dplyr::select(sample:DNAm, location, INS_ID, group)
    
    one_mei_eg <- mei_ins_only_DNAm_df %>%
      dplyr::filter(INS_ID == input$INS_id) %>%
      dplyr::mutate(location = coord) %>% 
      dplyr::select(sample:DNAm, location, INS_ID, group)
    
    mei_box_df <- rbind(one_mei_flank_eg, one_mei_eg) %>% 
      dplyr::mutate(
        region = ifelse(!location %in% c("flank_up", "flank_down"), "INS", location),
        region = factor(region, levels = c("flank_up", "INS", "flank_down")),
        location = factor(location, 
                          levels = c("flank_up", 
                                     unique(one_mei_eg$location), 
                                     "flank_down"))) %>% 
      dplyr::filter(group %in% input$compare_group,         # filter by selected groups
                    region %in% input$compare_location  # filter by selected locations
                    )
    
    ## plot
    
    # my_comparisons <- list(c("flank_up", unique(one_mei_eg$location)), 
    #                          c(unique(one_mei_eg$location), "flank_down"),
    #                          c("flank_up", "flank_down"))
    
    # --- Build dynamic comparison list using the *actual* INS coordinate ---
    ins_coord <- unique(one_mei_eg$location)
    my_comparisons <- list(
      c("flank_up", ins_coord),
      c(ins_coord, "flank_down"),
      c("flank_up", "flank_down")
    )
    my_comparisons <- Filter(function(pair) all(pair %in% mei_box_df$location), my_comparisons)
    
    ggplot(mei_box_df,
           aes(x = location, y = DNAm, fill = region)) +
      geom_boxplot() +
      # Add pairwise comparisons p-value
      stat_compare_means(comparisons = my_comparisons, size = 8, 
                         # label.y = c(1, 1.25, 1.5), 
                         label.y = seq(1, 1 + 0.25 * (length(my_comparisons) - 1), by = 0.25),
                         method.args = list(exact = FALSE),)+ 
      stat_compare_means(label.y = 1.9, size = 8, method.args = list(exact = FALSE))+     # Add global p-value
      facet_wrap(~group, nrow = 2)+
      ylim(c(0, 2))+
      theme_pubr(base_size = 20) +
      labs(title = paste0("Compare mean DNAm: ", input$INS_id),
           subtitle = "Flank regions: 1000 bp upstream and\ndownstream of EMI")+
      theme(
        # Plot title
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 18,  hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 16, face = "bold", colour = "black"),
        
        # Legend title
        legend.title = element_text(size = 16, face = "bold"),
        
        # Legend labels (tick values)
        legend.text = element_text(size = 16, face = "bold", colour = "black"),
        
        # Optional: adjust spacing between items
        legend.key.size = unit(0.8, "cm"),  # size of boxes
        
        # Optional: place inside or below plot
        legend.position = "bottom"
      )
  })

  
  # 🔴 Plot 3.1 for tab2: Overall boxplot (all INS, switched via button)
  output$overall_boxplot_compare1_ui <- renderPlot({
    # only generate after user clicks
    # req(input$all_emi_boxplot)  
    
    # Add a progress bar wrapper
    withProgress(message = "Rendering overall boxplot...", value = 0.1, {

    ## plot for all samples and all MEIs in this MEI group
    my_comparisons <- list(c("flank_up", "INS"), 
                           c("INS", "flank_down"),
                           c("flank_up", "flank_down"))
    
    # Optionally simulate some work:
    incProgress(0.2, detail = "Preparing data...")  # progress can be incremented
    
    all_mei_plot <- ggplot(all_mei_box_df,
           aes(x = region, y = DNAm)) +
      # geom_boxplot(fill = "cornsilk3") + # outlier.shape = NA
      geom_violin(fill = "cornsilk3") + # outlier.shape = NA
      stat_summary(fun = median, geom = "point", size = 5, color = "navy")+
      # Add pairwise comparisons p-value
      stat_compare_means(comparisons = my_comparisons, size = 8, 
                         label.y = c(1, 1.25, 1.5), method.args = list(exact = FALSE),)+ 
      # Add global p-value
      # stat_compare_means(label = "p.format", label.y = 1.9, size = 8, 
      #                    method.args = list(exact = FALSE))+     
      facet_wrap(~group, ncol = 2)+
      ylim(c(0, 1.8))+
      theme_pubr(base_size = 20) +
      labs(title = paste0("Compare mean DNAm (all MEIs): ", mei_group_all),
           subtitle = "Flank regions: 1000 bp upstream and downstream of EMIs")+
      theme(
        # Plot title
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, colour = "black"),
        plot.subtitle = element_text(size = 18,  hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 16, face = "bold", colour = "black"),
        
        # Legend title
        legend.title = element_text(size = 16, face = "bold"),
        
        # Legend labels (tick values)
        legend.text = element_text(size = 16, face = "bold", colour = "black"),
        
        # Optional: adjust spacing between items
        legend.key.size = unit(0.8, "cm"),  # size of boxes
        
        # Optional: place inside or below plot
        legend.position = "bottom")
    
    incProgress(0.9, detail = "Plotting...")  # progress can be incremented
    all_mei_plot
    })  # end withProgress
  })
  
    
  # 🔴 Plot 3.2: Overall boxplot (all INS, switched via button)
  output$overall_boxplot_compare2_ui <- renderPlot({
    # only generate after user clicks
    # req(input$all_emi_boxplot) 
    
    # Add a progress bar wrapper
    withProgress(message = "Rendering overall boxplot...", value = 0.1, {
      
    # Optionally simulate some work:
    incProgress(0.2, detail = "Preparing data...")  # progress can be incremented
 
    sample_all_mei_plot <- ggplot(all_mei_box_df,
           aes(x = sample, y = DNAm, fill = region)) +
      geom_boxplot(position = position_dodge(width = 0.75)) +
      stat_summary(fun = median, geom = "point", 
                   size = 3, color = "navy", 
                   position = position_dodge(width = 0.75))+
      # facet_wrap(~region, nrow=3)+
      theme_pubr(base_size = 20) +
      labs(title = paste0("All sample DNAm (all MEIs): ", mei_group_all),
           subtitle = "Flank regions: 1000 bp upstream and downstream of EMIs")+
      theme(
        # Plot title
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, colour = "black"),
        plot.subtitle = element_text(size = 18,  hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", colour = "black",
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 16, face = "bold", colour = "black"),
        
        # Legend title
        legend.title = element_text(size = 16, face = "bold"),
        
        # Legend labels (tick values)
        legend.text = element_text(size = 16, face = "bold", colour = "black"),
        
        # Optional: adjust spacing between items
        legend.key.size = unit(0.8, "cm"),  # size of boxes
        
        # Optional: place inside or below plot
        legend.position = "bottom"
      )
    
    incProgress(0.9, detail = "Plotting...")  # progress can be incremented
    sample_all_mei_plot
    })  # end withProgress
  })
    
}

shinyApp(ui = ui, server = server)








