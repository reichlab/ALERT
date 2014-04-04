shinyUI(pageWithSidebar(
    
    headerPanel("The ALERT Algorithm"),
    
    sidebarPanel(
        
        h4("Upload Data"),
        
        fileInput('file1', 'Choose CSV File (Must consist of only two columns, the first for the date and the second for cases)', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
        checkboxInput('header', 'Header', TRUE),
        radioButtons('sep', 'Separator', c(Comma=',',Semicolon=';', Tab='\t'),'Comma'),
        
        tags$hr(),
        
        h4("Choose Parameters"),
        
        sliderInput("firstMonth", "First Month of Flu Season", min=1, max=12, value=8),
        sliderInput("minWeeks", "Shortest Possible Flu Season (Weeks)", min=1, max=10, value=8),
        sliderInput("k", "+/- k Weeks", min=0, max=5, value=2),
        sliderInput("lag", "Days Between Cases Reports and Policy Action", min=1, max=21, value=7),
        sliderInput("target.dur", "Target Duration (Weeks)", min=10, max=20, value=12),
        sliderInput("target.pct", "Target % of Cases Covered", min=0, max=100, value=85, step=5),
        checkboxInput("allThresholds", "Use More Thresholds?", FALSE),
        
        tags$hr(),
        
        helpText("Created by Stephen A Lauer and Nicholas G Reich"),
        
        helpText(a("Send us your comments or feedback!", href="mailto:slauer@schoolph.umass.edu", target="_blank")),
        
        helpText(a("ALERT on GitHub", href="https://github.com/nickreich/ALERT", target="_blank"))
    ),
    
    mainPanel(
        tabsetPanel(
            tabPanel("Data Summary", dataTableOutput("summary"), 
                     tags$style(type="text/css", '#summary tfoot {display:none;}'), 
                     plotOutput("dataplot")),
            tabPanel("Performance Graphs", plotOutput("durplot"),
                     plotOutput("pctplot"))
        )
    )
))