require(ggplot2)
require(ALERT)
data(maskData)
#source("ALERT.R")
#load("maskData.RData")

shinyServer(function(input, output) {
    data <- reactive({
        inFile <- input$file1
        if(is.null(inFile)){
            ddd <- maskData
        }

        else
            ddd <- read.csv(inFile$datapath, header=input$header, sep=input$sep)
        
        colnames(ddd) <- c("Date", "Cases")
        ddd$Date <- as.Date(ddd$Date, "%m/%d/%y")
        as.data.frame(ddd)
    })
    
    output$dataplot <- renderPlot({
        data.plot <- ggplot(data=data()) + 
            geom_bar(aes(x=Date, y=Cases), stat="identity", fill="#0072B2", 
                     color="#0072B2") + 
            theme_bw() + 
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_text(face="bold", size=18),
                  axis.text.x  = element_text(size=16),
                  axis.text.y  = element_text(size=16)) +
            ylab("Influenza A case counts")
        print(data.plot)
    })
    
    threshdat <- reactive({
        tmp <- createALERT(data=data(), 
                           firstMonth=input$firstMonth, 
                           minWeeks=input$minWeeks, 
                           k=input$k, 
                           lag=input$lag, 
                           target.pct=input$target.pct/100, 
                           allThresholds=input$allThresholds
        )$out
        colnames(tmp) <- c("threshold", "average duration", "mean % cases captured", 
                           "min % cases captured", "max % cases captured", 
                           "sd % cases captured", "% peaks captured", 
                           "% peaks +/- k weeks captured", "mean # of weeks below threshold", 
                           "mean # of weeks longer than optimal")
        as.data.frame(round(tmp,1))
    })
    
    output$summary = renderDataTable({
        threshdat()
    }, options=list(bFilter=0, bSortClasses=TRUE, bProcessing=0, bPaginate=0, bInfo=0))        
    
    #output$threshplot <- renderPlot({
    #    p <- qplot(data=threshdat(), x=as.factor(threshold), 
    #               y=threshdat()[,"average duration"], geom="bar", stat="identity",
    #               xlab="Threshold", ylab="average duration", color="#999999") +
    #        theme_bw()
    #    y.bar <- ifelse(input$variable=="average duration", input$target.dur, 
    #                    input$target.pct)
    #    p <- p + geom_hline(aes(yintercept=y.bar), color="990000", linetype="dashed")
    #    print(p)
    #})
    
    output$threshplot <- renderPlot({
        threshdat <- threshdat()
        p <- qplot(data=threshdat, x=as.factor(threshold), 
                   y=threshdat[,input$variable], geom="bar", stat="identity",
                   xlab="Threshold", ylab=input$variable)
        y.bar <- ifelse(input$variable=="average duration", input$target.dur, input$target.pct)
        p <- p + geom_hline(aes(yintercept=y.bar), color="990000", linetype="dashed")
        print(p)
    })
})