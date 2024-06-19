library(shiny)
library(bslib)
library(ReturnCurves)
library(evd)
library(ismev)
library(latex2exp)
library(DT)
library(ggplot2)
library(dplyr)
library(gridExtra)

source("plots_eda.R")
source("plot_rc.R")
source("plot_adf.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Return Curve Estimation"),
  fileInput("data", label = "File input", accept = c(".csv", ".rds", ".txt"),
            buttonLabel = "Browse...", placeholder = "No file selected"),
  uiOutput("select_column_x"),
  uiOutput("select_column_y"),
  hr(),
  sidebarLayout(
    sidebarPanel(
      selectInput("analysis", "Choose Analysis:", 
                  choices = c("Exploratory Data Analysis", "Return Curve Estimation", "Angular Dependence Function")),
      uiOutput("analysisinputs")#,
      # tags$style(HTML('.irs--shiny .irs-line { background: linear-gradient(to right, red, green);}
      #               '))
    ),
    mainPanel(
      uiOutput("analysis")
    )
  )
)



server <- function(input, output, session) {

  output$analysisinputs <- renderUI({
    if (input$analysis == "Return Curve Estimation") {
      fluidRow(
        column(6,
          withMathJax(),
          sliderInput("lengthw", "Number of angles \\(\\omega\\)",
                      min = 101, max = 1001, step = 100, value = 101),
          hr(),
          numericInput("probability", "Curve survival probability \\(p\\)", value = 0.001,
                       min = 0, max = 1),
          hr(),
          radioButtons("method", "Method to estimate \\(\\lambda(\\omega)\\)",
                      choiceValues = list("hill", "cl"), choiceNames = list("Hill", "Composite Likelihood")),
          hr(),
          sliderInput("qmarg1", "Marginal quantile for the Marginal transformation for the first variable",
                      min = 0.01, max = 0.99, step = 0.01, value = 0.95),
          sliderInput("qmarg2", "Marginal quantile for the Marginal transformation for the second variable",
                      min = 0.01, max = 0.99, step = 0.01, value = 0.95),
          radioButtons("constrainedshape", "Constrained the shape parameter of the GPD fit",
                      choices = c(TRUE, FALSE), inline = T)
        ),
        column(6,
               withMathJax(),
               sliderInput("q", "Marginal quantile for the min-projection variable and/or Hill estimator",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               sliderInput("qalphas1", "Marginal quantile used for the conditional extremes model for the first variable",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               sliderInput("qalphas2", "Marginal quantile used for the conditional extremes model for the second variable",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               hr(),
               numericInput("k", "Polynomial degree", value = 7),
               hr(),
               radioButtons("constrained", "Incorporate knowledge of conditional extremes parameters",
                            choices = c(FALSE, TRUE), inline = T)
        ),
        hr(),
        column(12, style="background-color:#ff9b9b",
          withMathJax(),
          numericInput("tol", "Convergence tolerance for the composite maximum likelihood procedure",
                       value = 0.0001),
          hr(),
          numericInput("parinit", "Initial values for the parameters \\(\\beta\\)",
                       value = 0)
        ),
        column(6,
               hr(),
               actionButton("rcgof", "Goodness-of-fit"),
               hr(),
               uiOutput("rcgof_inputs")
        ),
        column(6,
               hr(),
               actionButton("unc", "Uncertainty"),
               hr(),
               uiOutput("uncertainty_inputs")
        )
      )
    }
    else if(input$analysis == "Angular Dependence Function") {
      fluidRow(
        column(6,
               withMathJax(),
               sliderInput("lengthw", "Number of angles \\(\\omega\\)",
                           min = 101, max = 1001, step = 100, value = 101),
               hr(),
               radioButtons("method", "Method to estimate \\(\\lambda(\\omega)\\)",
                            choiceValues = list("hill", "cl"), choiceNames = list("Hill", "Composite Likelihood")),
               hr(),
               sliderInput("qmarg1", "Marginal quantile for the Marginal transformation for the first variable",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               sliderInput("qmarg2", "Marginal quantile for the Marginal transformation for the second variable",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               radioButtons("constrainedshape", "Constrained the shape parameter of the GPD fit",
                            choices = c(TRUE, FALSE), inline = T)
        ),
        column(6,
               withMathJax(),
               sliderInput("q", "Marginal quantile for the min-projection variable and/or Hill estimator",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               sliderInput("qalphas1", "Marginal quantile used for the conditional extremes model for the first variable",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               sliderInput("qalphas2", "Marginal quantile used for the conditional extremes model for the second variable",
                           min = 0.01, max = 0.99, step = 0.01, value = 0.95),
               hr(),
               numericInput("k", "Polynomial degree", value = 7),
               hr(),
               radioButtons("constrained", "Incorporate knowledge of conditional extremes parameters",
                            choices = c(FALSE, TRUE), inline = T)
        ),
        hr(),
        column(12, style="background-color:#ff9b9b",
               withMathJax(),
               numericInput("tol", "Convergence tolerance for the composite maximum likelihood procedure",
                            value = 0.0001),
               hr(),
               numericInput("parinit", "Initial values for the parameters \\(\\beta\\)",
                            value = 0)
        ),
        hr(),
        actionButton("adfgof", "Goodness-of-fit"),
        hr(),
        uiOutput("adfgof_inputs")
      )
    }
  })
  
  rcunctoggleState <- reactiveVal(FALSE)
  rcplotOutput <- reactiveVal(NULL)
  
  observeEvent(input$unc, {
    rcunctoggleState(!rcunctoggleState())
    output$uncertainty_inputs <- renderUI({
      if (rcunctoggleState()) {
        tagList(
          withMathJax(),
          numericInput("blocksize", "Size of blocks for block bootstrap",
                       value = 1),
          hr(),
          numericInput("nboot", "Number of bootstrap samples",
                       value = 50),
          hr(),
          numericInput("nangles", "Number of angles in \\((0, \\pi/2)\\)",
                       value = 150),
          hr(),
          sliderInput("alpha", "Significance level for the \\((1-\\alpha)\\)% CI",
                      min = 0.01, max = 0.99, step = 0.05, value = 0.05)
        )
      }
    })
  })
  
  rcgoftoggleState <- reactiveVal(FALSE)
  
  observeEvent(input$rcgof, {
    rcgoftoggleState(!rcgoftoggleState())
    output$rcgof_inputs <- renderUI({
      if (rcgoftoggleState()) {
        tagList(
          withMathJax(),
          numericInput("blocksize", "Size of blocks for block bootstrap",
                       value = 1),
          hr(),
          numericInput("nboot", "Number of bootstrap samples",
                       value = 250),
          hr(),
          numericInput("nangles", "Number of angles in \\((0, \\pi/2)\\)",
                       value = 150),
          hr(),
          sliderInput("alpha", "Significance level for the \\((1-\\alpha)\\)% CI",
                      min = 0.01, max = 0.99, step = 0.05, value = 0.05)
        )
      }
    })
    output$rcunc <- renderPlot({
      if (rcunctoggleState()) {
        req(rcplotOutput())
        
        uncrcplot(rcplotOutput(), input$blocksize, input$nboot, input$nangles, input$alpha)
      }
    })
    output$rcgof <- renderPlot({
      if (rcgoftoggleState()) {
        req(rcplotOutput())
        
        gofrcplot(rcplotOutput(), input$blocksize, input$nboot, input$nangles, input$alpha)
      }
    })
  })
  
  adfgoftoggleState <- reactiveVal(FALSE)
  adfplotOutput <- reactiveVal(NULL)
  
  observeEvent(input$adfgof, {
    adfgoftoggleState(!adfgoftoggleState())
    output$adfgof_inputs <- renderUI({
      if (adfgoftoggleState()) {
        tagList(
          withMathJax(),
          numericInput("ray", "Ray \\(\\omega\\)",
                       value = 0.5, min = 0, max = 1),
          hr(),
          numericInput("blocksize", "Size of blocks for block bootstrap",
                       value = 1),
          hr(),
          numericInput("nboot", "Number of bootstrap samples",
                       value = 250),
          hr(),
          sliderInput("alpha", "Significance level for the \\((1-\\alpha)\\)% tolerance intervals",
                      min = 0.01, max = 0.99, step = 0.05, value = 0.05),
          hr()
        )
      }
    })
    output$adfgof <- renderPlot({
      if (adfgoftoggleState()) {
        req(adfplotOutput())
        
        gofadfplot(adfplotOutput(), input$ray, input$blocksize, input$nboot, input$alpha)
      }
    })
  })
  
  output$analysis <- renderUI({
    if (input$analysis == "Exploratory Data Analysis") {
      fluidRow(
        column(6,
               card(card_header("Histogram"),
                    plotOutput("hist")
               ),
               card(card_header("Autocorrelation Plots"),
                    plotOutput("acf")
               )

        ),
        column(6,
               card(card_header("Time Series"),
                    plotOutput("timeseries")
               ),
               card(card_header("Joint"),
                    plotOutput("joint")
               )
        )
      )
    }else if (input$analysis == "Return Curve Estimation") {
      column(12,
      card(card_header("Return Curve Estimation"),
           plotOutput("rc")
      ),
      if (rcunctoggleState()) {
        card(card_header("Uncertainty"),
             plotOutput("rcunc")
             )
          },
      if (rcgoftoggleState()) {
        card(card_header("Goodness-of-fit"),
             plotOutput("rcgof")
        )
      }
    )
    }else if (input$analysis == "Angular Dependence Function") {
      column(12, 
        card(card_header("Angular Dependence Function"),
             plotOutput("adf")
        ),
        if (adfgoftoggleState()) {
          card(card_header("Goodness-of-fit"),
               plotOutput("adfgof")
          )
        }
      )
    }
  })
  
  data <- reactive({
    req(input$data)
    ext <- tools::file_ext(input$data$name)
    switch(ext,
           csv = read.csv(input$data$datapath),
           rds = readRDS(input$data$datapath),
           txt = read.table(input$data$datapath, header = TRUE, sep = ";", dec = "."),
           stop("Unsupported file type")
    )
  })
  
  observe({
    req(data())
    numeric_cols <- names(data())[sapply(data(), is.numeric)]
    updateSelectInput(session, "colX", choices = numeric_cols)
    updateSelectInput(session, "colY", choices = numeric_cols)
  })
  
  output$select_column_x <- renderUI({
    req(data())
    selectInput("colX", "Select first variable", choices = NULL)
  })
  
  output$select_column_y <- renderUI({
    req(data())
    selectInput("colY", "Select second variable", choices = NULL)
  })
  
  output$hist <- renderPlot({
    # req(data())
    req(data(), input$colX, input$colY)
    validate(
      need(input$colX %in% names(data()), "Select a valid first variable"),
      need(input$colY %in% names(data()), "Select a valid second variable"),
      need(is.numeric(data()[[input$colX]]), "First variable must be numeric"),
      need(is.numeric(data()[[input$colY]]), "Second variable must be numeric")
    )
    histplot(data(), input$colX, input$colY)
  })
  
  output$timeseries <- renderPlot({
    # req(data())
    req(data(), input$colX, input$colY)
    validate(
      need(input$colX %in% names(data()), "Select a valid first variable"),
      need(input$colY %in% names(data()), "Select a valid second variable"),
      need(is.numeric(data()[[input$colX]]), "First variable must be numeric"),
      need(is.numeric(data()[[input$colY]]), "Second variable must be numeric")
    )
    timeseriesplot(data(), input$colX, input$colY)
  })
  
  output$acf <- renderPlot({
    # req(data())
    req(data(), input$colX, input$colY)
    validate(
      need(input$colX %in% names(data()), "Select a valid first variable"),
      need(input$colY %in% names(data()), "Select a valid second variable"),
      need(is.numeric(data()[[input$colX]]), "First variable must be numeric"),
      need(is.numeric(data()[[input$colY]]), "Second variable must be numeric")
    )
    acfplot(data(), input$colX, input$colY)
  })
  
  output$joint <- renderPlot({
    # req(data())
    req(data(), input$colX, input$colY)
    validate(
      need(input$colX %in% names(data()), "Select a valid first variable"),
      need(input$colY %in% names(data()), "Select a valid second variable"),
      need(is.numeric(data()[[input$colX]]), "First variable must be numeric"),
      need(is.numeric(data()[[input$colY]]), "Second variable must be numeric")
    )
    jointplot(data(), input$colX, input$colY)
  })
  
  output$rc <- renderPlot({
    # req(data())
    req(data(), input$colX, input$colY)
    validate(
      need(input$colX %in% names(data()), "Select a valid first variable"),
      need(input$colY %in% names(data()), "Select a valid second variable"),
      need(is.numeric(data()[[input$colX]]), "First variable must be numeric"),
      need(is.numeric(data()[[input$colY]]), "Second variable must be numeric")
    )
    
    rc_data <- rcplot(data(), input$colX, input$colY, input$qmarg1, input$qmarg2, input$constrainedshape,
           input$lengthw, input$probability, input$method, input$q, input$qalphas1, 
           input$qalphas2, input$k, input$constrained, input$tol, input$parinit)
    rcplotOutput(rc_data)
    plot(rc_data)
  })
  
  output$adf <- renderPlot({
    # req(data())
    req(data(), input$colX, input$colY)
    validate(
      need(input$colX %in% names(data()), "Select a valid first variable"),
      need(input$colY %in% names(data()), "Select a valid second variable"),
      need(is.numeric(data()[[input$colX]]), "First variable must be numeric"),
      need(is.numeric(data()[[input$colY]]), "Second variable must be numeric")
    )
    
    adf_result <- adfplot(data(), input$colX, input$colY, input$qmarg1, input$qmarg2, input$constrainedshape,
           input$lengthw, input$method, input$q, input$qalphas1, 
           input$qalphas2, input$k, input$constrained, input$tol, input$parinit)
    
    adfplotOutput(adf_result)
    plot(adf_result)
  })
  
}

shinyApp(ui = ui, server = server)



