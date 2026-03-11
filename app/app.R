# ==============================================================================
# Script: app.R (Shiny Clinician Dashboard)
# Purpose: Interactive web interface to predict Immunotherapy Response
# Author: Multi-Omics Pipeline
# ==============================================================================

if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if (!requireNamespace("bslib", quietly = TRUE)) install.packages("bslib")
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidymodels", quietly = TRUE)) install.packages("tidymodels")
if (!requireNamespace("vip", quietly = TRUE)) install.packages("vip")

library(shiny)
library(bslib)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidymodels)
library(vip)

# ------------------------------------------------------------------------------
# 1. Global Setup
# ------------------------------------------------------------------------------
# Load the pre-trained workflow
model_path <- "models/rf_model.rds"
model_loaded <- FALSE
rf_workflow <- NULL

if (file.exists(model_path)) {
  rf_workflow <- readRDS(model_path)
  model_loaded <- TRUE
}

# ------------------------------------------------------------------------------
# 2. UI Definition
# ------------------------------------------------------------------------------
ui <- page_navbar(
  title = "ImmunoPredict: Biomarker Response Dashboard",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2c3e50"),
  
  # Tab 1: Patient Prediction
  nav_panel(title = "Clinical Prediction",
    layout_sidebar(
      sidebar = sidebar(
        width = 400,
        title = "Patient Biomarker Profile",
        p("Enter patient omics values below to predict response to Immune Checkpoint Blockade."),
        hr(),
        textInput("patient_id", "Patient ID", value = "NEW-PT-001"),
        numericInput("tmb_input", "Tumor Mutational Burden (Mutations/Mb)", 
                     value = 10.5, min = 0, max = 150, step = 0.5),
        helpText("Higher TMB (>10) is typically associated with better response."),
        
        sliderInput("ifng_input", "IFN-gamma Signature Score (logCPM)", 
                    min = -5, max = 15, value = 5.0, step = 0.5),
        helpText("Reflects inflamed tumor microenvironment."),
        
        hr(),
        actionButton("predict_btn", "Predict Response", class = "btn-primary w-100", icon = icon("stethoscope"))
      ),
      
      # Main Content Panel
      card(
        card_header("Prediction Results"),
        card_body(
          uiOutput("prediction_ui"),
          br(),
          plotlyOutput("gauge_plot")
        )
      )
    )
  ),
  
  # Tab 2: Model Interpretability
  nav_panel(title = "Model Explainability (SHAP/VIP)",
    card(
      card_header("Global Feature Importance"),
      card_body(
        p("This plot shows which multi-omics features the Random Forest relies on most heavily to make predictions across the training cohort."),
        plotOutput("vip_plot", height = "500px")
      )
    )
  )
)

# ------------------------------------------------------------------------------
# 3. Server Logic
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---------------------------------------------------
  # Prediction Logic
  # ---------------------------------------------------
  observeEvent(input$predict_btn, {
    
    if (!model_loaded) {
      output$prediction_ui <- renderUI({
        div(class = "alert alert-danger", 
            "Error: Model file (models/rf_model.rds) not found. Please run 03_train_models.R first.")
      })
      return()
    }
    
    # Construct a single-row dataframe mimicking the format used in training
    new_patient <- tibble(
      Patient = input$patient_id,
      TMB = input$tmb_input,
      IFNg_Score = input$ifng_input
    )
    
    # Validate Inputs
    if (is.na(new_patient$TMB) || is.na(new_patient$IFNg_Score)) {
      output$prediction_ui <- renderUI({
        div(class = "alert alert-warning", "Please ensure all numeric fields are filled correctly.")
      })
      return()
    }
    
    # Predict Probability
    pred_prob <- predict(rf_workflow, new_data = new_patient, type = "prob")
    prob_responder <- pred_prob$.pred_Responder[1]
    
    # Predict Class
    pred_class <- predict(rf_workflow, new_data = new_patient, type = "class")
    class_result <- as.character(pred_class$.pred_class[1])
    
    # ---------------------------------------------------
    # Render UI Feedback
    # ---------------------------------------------------
    alert_class <- ifelse(class_result == "Responder", "alert-success", "alert-danger")
    
    output$prediction_ui <- renderUI({
      div(class = sprintf("alert %s", alert_class),
          h4(sprintf("Prediction: %s", class_result)),
          p(sprintf("Probability of Response (CR/PR): %.1f%%", prob_responder * 100)),
          hr(),
          p(sprintf("Patient %s has a %s likelihood of benefiting from immunotherapy based on current biomarkers.", 
                    new_patient$Patient, 
                    ifelse(prob_responder > 0.5, "high", "low")))
      )
    })
    
    # ---------------------------------------------------
    # Render Gauge Plot
    # ---------------------------------------------------
    output$gauge_plot <- renderPlotly({
      plot_ly(
        domain = list(x = c(0, 1), y = c(0, 1)),
        value = prob_responder * 100,
        title = list(text = "Response Probability"),
        type = "indicator",
        mode = "gauge+number",
        gauge = list(
          axis = list(range = list(NULL, 100)),
          bar = list(color = ifelse(prob_responder > 0.5, "forestgreen", "darkred")),
          steps = list(
            list(range = c(0, 50), color = "lightgray"),
            list(range = c(50, 100), color = "whitesmoke")
          ),
          threshold = list(
            line = list(color = "black", width = 4),
            thickness = 0.75,
            value = 50
          )
        )
      ) %>% layout(margin = list(t = 50, b = 0, l = 0, r = 0))
    })
  })
  
  # ---------------------------------------------------
  # Model Explainability (VIP Plot)
  # ---------------------------------------------------
  output$vip_plot <- renderPlot({
    if (model_loaded) {
      # Extract underlying parsnip model
      extracted_model <- extract_fit_parsnip(rf_workflow)
      
      # Plot variable importance
      vip(extracted_model, 
          geom = "col", 
          aesthetics = list(fill = "#2c3e50", alpha = 0.8)) +
        theme_minimal(base_size = 15) +
        labs(
          title = "Random Forest Feature Importance",
          subtitle = "Gini Impurity Decrease",
          x = "Biomarker Feature",
          y = "Importance Score"
        )
    } else {
      # Dummy plot if model isn't loaded
      ggplot() + 
        annotate("text", x=0, y=0, label="Model not loaded. Please train the model first.", size=6) +
        theme_void()
    }
  })
}

# ------------------------------------------------------------------------------
# 4. App Launch
# ------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
