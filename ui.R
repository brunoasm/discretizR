
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("discretizR"),

  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      fileInput("rawtable", label = "Upload measurement table"), #upload table with measurements
      fileInput("chartable", label = "Upload character table"), #upload table with measurements
      uiOutput('body_size_options'),
      uiOutput('select_variable'),
      uiOutput('output_for_mrbayes'),
      uiOutput('output_for_tnt')
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel(
          'Regression Plots',
          fluidRow(
            uiOutput('select_graph_color')
          ),
          fluidRow(textOutput('sd_message')),
          fluidRow(plotOutput('regression_size')),
          fluidRow(textOutput('lm_message')),
          fluidRow(plotOutput('residuals_species'))
        ),
        tabPanel(
          'Gaussian mixture disgnostics',
          fluidRow(
            column(12,
                   fluidRow(plotOutput('densiplot'))),
            column(12,
                   'Check the help of package ',tags$a(href='https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html','mclust'),'for explanation on plots below.'),
            column(6,
              fluidRow(
                plotOutput('mclust_density'),
                plotOutput('BICplot')
              )
            ),
            column(6,
              fluidRow(
                plotOutput('qqplot'),
                plotOutput('cdfplot')
              )
            )
          )
        ),
        tabPanel(
          'Character correlation',
          fluidRow(column(6, sliderInput('mincorr',
                               'Choose the minimum (absolute) correlation to cluster characters:',
                               min=0,
                               max=1,
                               value=0.1)),
                   column(6, uiOutput('corr_group_selector'))),
          fluidRow(column(12, uiOutput('corr_char_selector'))),
          fluidRow(column(12, textOutput('zero_message'))),
          fluidRow(column(12, tableOutput('corr_table'))),
          fluidRow(column(12, tableOutput('pvalue_table')))
        )
      )
    )
  )
))
