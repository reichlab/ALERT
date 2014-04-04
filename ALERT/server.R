require(ggplot2)
require(ALERT)
require(reshape2)
data(maskData)
#source("ALERT.R")
#load("maskData.RData")
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
            geom_bar(aes(x=Date, y=Cases), stat="identity", fill="#0072B2", color="#0072B2") + 
            scale_y_continuous(name="Case Counts") +
            theme_bw() + 
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_text(face="bold", size=18),
                  axis.text.x  = element_text(size=16),
                  axis.text.y  = element_text(size=16))
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
        as.data.frame(round(tmp,1))
    })
    
    output$summary = renderDataTable({
        threshdat <- threshdat()
        colnames(threshdat) <- c("Threshold", "Average Duration", "Mean % Cases Captured", 
                           "Min % Cases Captured", "Max % Cases Captured", 
                           "SD % Cases Captured", "% Peaks Captured", 
                           "% Peaks +/- k Weeks Captured", "Mean # of Weeks Below Threshold", 
                           "Mean # of Weeks Longer Than Optimal")
        threshdat
    }, options=list(bFilter=0, bSortClasses=TRUE, bProcessing=0, bPaginate=0, bInfo=0))        
    
    output$durplot <- renderPlot({
        threshdat <- threshdat()
        
        p <- qplot(data=threshdat, x=threshold, y=mean.dur) + 
            geom_point(color=cbPalette[6], size=4) +
            geom_line(color=cbPalette[6], size=2) +
            scale_x_continuous(name = "", breaks=threshdat$threshold) + 
            scale_y_continuous(name = "Average Duration (Weeks)", breaks=seq(0,20,4), limits=c(0,20)) +
            geom_hline(aes(yintercept=input$target.dur), color=cbPalette[7], linetype="dashed", size=1) +
            #ggtitle("Average Duration and Mean Percentage of Cases Captured by Threshold") +
            theme_bw() + theme(axis.text.x = element_text(size=16),
                               axis.text.y = element_text(size=16),
                               axis.title.y = element_text(size=18, face="bold"),
                               plot.title = element_text(size=20, face="bold"))
        
        print(p)   
    })
    
    output$pctplot <- renderPlot({
        threshdat <- threshdat()
        
        p <- qplot(data=threshdat, x=threshold, y=mean.pct.cases.captured) +
            geom_point(color=cbPalette[4], size=4) +
            geom_line(color=cbPalette[4], size=2) +
            scale_x_continuous(name = "Threshold", breaks=threshdat$threshold) + 
            scale_y_continuous(name = "Mean % Cases Captured", breaks=seq(0,100,20), limits=c(0,100)) +
            geom_hline(aes(yintercept=input$target.pct), color=cbPalette[7], linetype="dashed", size=1) +
            theme_bw() + theme(axis.title.x = element_text(size=18, face="bold"),
                               axis.text.x = element_text(size=16),
                               axis.text.y = element_text(size=16),
                               axis.title.y = element_text(size=18, face="bold"))
        
        print(p)   
    })
})