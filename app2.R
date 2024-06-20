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
library(shinydashboard)

source("plots_eda.R")
source("plot_rc.R")
source("plot_adf.R")

ui <- dashboardPage(
  dashboardHeader(title = "Return Curves"),
  dashboardSidebar(
    tags$style(HTML('.irs--shiny .irs-line { background: linear-gradient(to right, #ff0e0e 0%, #37b61e 100%);}
                       .irs--shiny .irs-bar { background: none;}
                    ')),
    sidebarMenu(
      fileInput("data", label = "File input", accept = c(".csv", ".rds", ".txt"), buttonLabel = "Browse...", placeholder = "No file selected"),
      uiOutput("select_column_x"),
      uiOutput("select_column_y"),
      menuItem("Exploratory Data Analysis",
               tabName = "eda", icon = icon("eye")),
      menuItem("Return Curve Estimation",
               tabName = "retcurve", icon = icon("chart-simple")),
      menuItem("Angular Dependence Function",
               tabName = "adf", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "eda",
              fluidRow(box(dataTableOutput("edatable"), width = 12)),
              fluidRow(box(plotOutput("hist")),
                       box(plotOutput("acf"))),
              fluidRow(box(plotOutput("timeseries")),
                       box(plotOutput("joint")))
      ),
      tabItem(tabName = "retcurve",
              withMathJax(),
              fluidRow(
                column(6,
                       box(sliderInput("rcqmarg1", "Marginal quantile for the Marginal transformation for the first variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           sliderInput("rcqmarg2", "Marginal quantile for the Marginal transformation for the second variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           radioButtons("rcconstrainedshape", "Constrained the shape parameter of the GPD fit",
                                        choices = c(TRUE, FALSE), inline = T),
                           numericInput("probability", "Curve survival probability \\(p\\)", value = 0.001,
                                        min = 0, max = 1, step = 0.001)),
                       box(numericInput("rclengthw", "Number of rays \\(\\omega\\)",
                                        min = 51, max = 1001, value = 101, step = 100),
                           radioButtons("rcmethod", "Method to estimate \\(\\lambda(\\omega)\\)",
                                        choiceValues = list("hill", "cl"), choiceNames = list("Hill", "Composite Likelihood")),
                           sliderInput("rcq", "Marginal quantile for the min-projection variable and/or Hill estimator",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           numericInput("rck", "Polynomial degree", value = 7, min = 1))
                ),
                column(6,
                       box(radioButtons("rcconstrained", "Incorporate knowledge of conditional extremes parameters",
                                        choices = c(FALSE, TRUE), inline = T),
                           sliderInput("rcqalphas1", "Marginal quantile used for the conditional extremes model for the first variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           sliderInput("rcqalphas2", "Marginal quantile used for the conditional extremes model for the second variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95)),
                       box(title = "Recommended not to change", status = "warning", solidHeader = T,
                           numericInput("rctol", "Convergence tolerance for the composite maximum likelihood procedure",
                                        value = 0.0001, min = 0),
                           numericInput("rcparinit", "Initial values for the parameters \\(\\beta\\)",value = 0))
                )),
                box(plotOutput("rc"), width = 12),
                fluidRow(
                  column(12,
                         box(actionButton("unc", "Uncertainty"),
                             uiOutput("uncertainty_inputs"),
                             plotOutput("rcunc")),
                         box(actionButton("rcgof", "Goodness-of-fit"),
                             uiOutput("rcgof_inputs"),
                             plotOutput("rcgof")))
                )
      ),
      tabItem(tabName = "adf",
              withMathJax(),
              fluidRow(
                column(6,
                       box(sliderInput("qmarg1", "Marginal quantile for the Marginal transformation for the first variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           sliderInput("qmarg2", "Marginal quantile for the Marginal transformation for the second variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           radioButtons("constrainedshape", "Constrained the shape parameter of the GPD fit",
                                        choices = c(TRUE, FALSE), inline = T)),
                       box(numericInput("lengthw", "Number of rays \\(\\omega\\)",
                                        min = 51, max = 1001, value = 101, step = 100),
                           radioButtons("method", "Method to estimate \\(\\lambda(\\omega)\\)",
                                        choiceValues = list("hill", "cl"), choiceNames = list("Hill", "Composite Likelihood")),
                           sliderInput("q", "Marginal quantile for the min-projection variable and/or Hill estimator",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           numericInput("k", "Polynomial degree", value = 7, min = 1))
                ),
                column(6,
                       box(radioButtons("constrained", "Incorporate knowledge of conditional extremes parameters",
                                        choices = c(FALSE, TRUE), inline = T),
                           sliderInput("qalphas1", "Marginal quantile used for the conditional extremes model for the first variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95),
                           sliderInput("qalphas2", "Marginal quantile used for the conditional extremes model for the second variable",
                                       min = 0.01, max = 0.99, step = 0.01, value = 0.95)),
                       box(title = "Recommended not to change", status = "warning", solidHeader = T,
                           numericInput("tol", "Convergence tolerance for the composite maximum likelihood procedure",
                                        value = 0.0001, min = 0),
                           numericInput("parinit", "Initial values for the parameters \\(\\beta\\)",value = 0))
                )),
              fluidRow(
                box(plotOutput("adfplot")),
                box(actionButton("adfgof", "Goodness-of-fit"),
                    uiOutput("adfgof_inputs"),
                    plotOutput("adfgof"))
              )
      )
    )
  )
)


server <- function(input, output, session) {
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
  
  output$edatable <- renderDataTable({
    req(data())
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
  
  
  rcunctoggleState <- reactiveVal(FALSE)
  rcplotOutput <- reactiveVal(NULL)
  
  observeEvent(input$unc, {
    rcunctoggleState(!rcunctoggleState())
    output$uncertainty_inputs <- renderUI({
      if (rcunctoggleState()) {
        tagList(
          withMathJax(),
          numericInput("rcblocksize", "Size of blocks for block bootstrap",
                       value = 1, min = 1),
          numericInput("rcnboot", "Number of bootstrap samples",
                       value = 50, min = 1),
          numericInput("rcnangles", "Number of rays in \\((0, \\pi/2)\\)",
                       value = 150, min = 1),
          sliderInput("rcalpha", "Significance level for the \\((1-\\alpha)\\)% CI",
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
          numericInput("rcgofblocksize", "Size of blocks for block bootstrap",
                       value = 1, min = 1),
          numericInput("rcgofnboot", "Number of bootstrap samples",
                       value = 250, min = 1),
          numericInput("rcgofnangles", "Number of rays in \\((0, \\pi/2)\\)",
                       value = 150, min = 1),
          sliderInput("rcgofalpha", "Significance level for the \\((1-\\alpha)\\)% CI",
                      min = 0.01, max = 0.99, step = 0.05, value = 0.05)
        )
      }
    })
    output$rcunc <- renderPlot({
      if (rcunctoggleState()) {
        req(rcplotOutput())
        uncrcplot(rcplotOutput(), input$rcblocksize, input$rcnboot, input$rcnangles, input$rcalpha)
      }
    })
    output$rcgof <- renderPlot({
      if (rcgoftoggleState()) {
        req(rcplotOutput())
        gofrcplot(rcplotOutput(), input$rcgofblocksize, input$rcgofnboot, input$rcgofnangles, input$rcgofalpha)
      }
    })
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
    
    rc_data <- rcplot(data(), input$colX, input$colY, input$rcqmarg1, input$rcqmarg2, input$rcconstrainedshape,
                      input$rclengthw, input$rcprobability, input$rcmethod, input$rcq, input$rcqalphas1, 
                      input$rcqalphas2, input$rck, input$rcconstrained, input$rctol, input$rcparinit)
    rcplotOutput(rc_data)
    plot(rc_data)
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
          numericInput("blocksize", "Size of blocks for block bootstrap",
                       value = 1),
          numericInput("nboot", "Number of bootstrap samples",
                       value = 250),
          sliderInput("alpha", "Significance level for the \\((1-\\alpha)\\)% tolerance intervals",
                      min = 0.01, max = 0.99, step = 0.05, value = 0.05),
          )
      }
    })
  })
      
  output$adfplot <- renderPlot({
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

