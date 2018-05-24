library(shiny)
library(glmnet)
source("./helpers.R")
source_data <- readRDS("./shiny_metareg_data.rds")


ui <- fluidPage(

  # Application title
  titlePanel("Interacting with Lasso Regressions"),

  # Sidebar with a slider input for a particular lambda value
  sidebarLayout(
    sidebarPanel(
      sliderInput("lambda",
        "Lambda:",
        min = 0,
        max = 100,
        value = 20,
        animate = TRUE),
      # add selectInput for dataset
      # Put helptext in here when it runs
      selectInput("dataset",
        "Choose a study:",
        choices = c("toy example", "Lightart", "Murphy",
          "Liu", "Ripke", "Ripke - NULL"),
        selected = "Murphy")
    ),

    # Show a plot of the generated distribution
    mainPanel(

      # next row is parallel coordinates plots which can be brushed ----
      fluidRow(
        column(6,
          plotOutput("yvsyhat")),
        column(6,
          plotOutput("pc1vspc2"))
      ),
      fluidRow(

        column(2,
          h2("Lambda"),
          textOutput("actualLam")),

        column(1,
          h2("R^2"),
          textOutput("R2")),
        column(4,
          h2("Selected Variables"),
          textOutput("selVars")
        ),
        column(4,
          h2("Rank of Active Set:"),
          textOutput("actRank")),
        fluidRow(
          column(6,
            plotOutput("lasso_path")),
          column(6,
            plotOutput("pcp_extra"))
        ),
        fluidRow(
          column(10,
            tableOutput("extra_summary")
          )
        )

      )
    )
  )
)


# Define server logic required to list the selected variables
server <- function(input, output) {

  datasetChoice <-reactive({
    switch(input$dataset,
      "toy example" = "toy",
      "Lightart" = "lightart",
      "Murphy" = "murphy",
      "Liu" = "liu",
      "Ripke" = "ripke",
      "Ripke - NULL" = "ripke_null")
    })


  datasetX <- reactive(source_data[[datasetChoice()]][["X"]])

  datasetY <- reactive(source_data[[datasetChoice()]][["Y"]])

  datasetNames <- reactive(source_data[[datasetChoice()]][["varnames"]])

  lasso_Fit <-  reactive(glmnet(x = datasetX(), y = datasetY()))
#

  runlambda <- reactive({

    ifelse((input$lambda > length(lasso_Fit()$lambda)), length(lasso_Fit()$lambda),
      input$lambda)

  })
#
  active_set <- reactive({
    activeset_names(lasso_Fit(), lambda_pos = runlambda(), datasetX())
  })
#
  summary_corr <- reactive({
    glmnet_summextra(lasso_Fit(), runlambda(), datasetX(), datasetY())
  })
#
   output$R2 <- renderText(lasso_Fit()$dev.ratio[[runlambda()]])
   output$actualLam <- renderText(lasso_Fit()$lambda[[runlambda()]])
#
   output$yvsyhat <- renderPlot({
     plot(datasetY(), glmnet_resids(lasso_Fit(), lambda_pos = runlambda(),
       datasetX(), datasetY()))
     abline(0,1)
     abline(0,0)
   })

    output$selVars <- renderText(active_set())
    output$actRank <- renderText(
       Matrix::rankMatrix(datasetX()[, active_set()])[[1]])

     output$pc1vspc2 <- renderPlot({
       active_pcs(datasetX(), active_set = active_set(),
         study_names = datasetNames())
     })
#   # #
#   #
#   #
   output$extra_summary <- renderTable(summary_corr()$extra_summary)
   output$pcp_extra <- renderPlot(summary_corr()$pcp_extra)
   output$lasso_path <- renderPlot({
     plot(lasso_Fit(), label = TRUE, xvar = "lambda")
     abline(v = log(lasso_Fit()$lambda[runlambda()]), lwd = 2, lty = 3,
       col = "orange")
   })
}

# Run the application
shinyApp(ui = ui, server = server)

