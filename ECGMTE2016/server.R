
# This is the server logic for a Shiny web application.

library(shiny)
library(grid)
library(ggplot2)


shinyServer(function(input, output) {

        # Default Parameters
        # ==========
        n <- 100                        ## resolution of ECG
        x <- seq(0, 6, by=0.01)         ## ECG strip size (standard)
        x2 <- seq(0, 24, by=0.01)       ## ECG strip size (long lead)
        leadnames <- c("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")
        
        small_h <- 0.2                 ## small square height
        small_w <- 0.2                 ## small square width
        big_h <- 5
        big_w <- 5
        
        strip_x <- 6                    ## individual strip length
        strip_y <- 4                    ## individual strip height

        ## set/ change default numbers here (here is taken to be lead II)
        ecgp <- c(0.25,0.09,0.16)
        ecgq <- c(0.025,0.066,0.166)
        ecgqrs <- c(1.6,0.11,0)
        ecgs <- c(0.025,0.66,0.09)
        ecgt <- c(0.35,0.144,0.2)
        ecgu <- c(0.015,0.0476,0.433)

        ## correct the time (t)
        ## ==========
        ##
        ## leadII is the variable that contains the default values of the ECG that all
        ## the other ECG leads are based upon with the time t corrected
        ## it has a total of 18 variables: 6 waves (p, q, qrs, s, t, u) 
        ## their 3 individual wave characteristics (a, d, t)
        ##
        leadII <- c(ecgp, ecgq, ecgqrs, ecgs, ecgt, ecgu)
        leadII[c(12,15,18)] <- -1*(leadII[c(12,15,18)])
        leadII[15] <- leadII[15] - 0.045

        ## corrected_axis is the default for a,d,t of all the waves of all 12 leads
        default_axis <- cbind(leadII, leadII, leadII, leadII, leadII, leadII, leadII, leadII, leadII, leadII, leadII, leadII)
        colnames(default_axis) <- leadnames
        rownames(default_axis) <- c("a_p","d_p","t_p","a_q","d_q","t_q","a_qrs","d_qrs","t_qrs",
                                      "a_s","d_s","t_s","a_t","d_t","t_t","a_u","d_u","t_u")


        ## Default ECG Themes
        ecgtheme <- theme(axis.text.y=element_blank(),
                          axis.text.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.title.x=element_blank(),
                          plot.margin = unit(c(0,-0.6,-0.6,-0.6),"cm"),
                          panel.background=element_rect(fill="white"), 
                          panel.grid.major=element_line(colour="red"),
                          panel.grid.minor=element_line(colour="pink"))      
        
        
        # The Functions used to Draw the ECG
        # =========
        
        
        # FUNCTION: theaxisfactor(front_axis, coronal_axis)
        # =======
        # function returns the axis factors for 12 leads given the,
        # front_axis (qrs_axis) AND the coronal_axis (uni_axis)
        #
        # I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6
        # 0, 60, 120, -150, -30, 90, 
        
        theaxisfactor <- function(front_axis, coronal_axis) {
                return( c(cos(pi*(front_axis-0)/180),cos(pi*(front_axis-60)/180),cos(pi*(front_axis-120)/180),
                          cos(pi*(front_axis+150)/180),cos(pi*(front_axis+30)/180),cos(pi*(front_axis-90)/180),
                          sin(pi*(coronal_axis-90)/180),sin(pi*(coronal_axis-80)/180),sin(pi*(coronal_axis-70)/180),
                          sin(pi*(coronal_axis-60)/180),sin(pi*(coronal_axis-30)/180),sin(pi*(coronal_axis-0)/180)))                
        }

        
        ## Function: Generic WAVE1 (p, t, u) function
        ## -----------
        ## draws the generic wave1 in the ECG
        ## 
        type1wav <- function(x,a,d,t,l) {
                x=x+t
                b=(2*l)/d
                w1=1/l
                w2=0
                
                for(i in 1:n) {
                        harm1=(((sin((pi/(2*b))*(b-(2*i))))/(b-(2*i))+(sin((pi/(2*b))*(b+(2*i))))/(b+(2*i)))*(2/pi))*cos((i*pi*x)/l)             
                        w2=w2+harm1
                }
                
                wav1 <- w1+w2
                return(a*wav1)
        }
        
        ## Function: Generic WAVE2 (q, qrs, s) function
        ## -----------
        ## draws the generic wave2 in the ECG
        ##
        type2wav <- function(x,a,d,t,l) {
                b=(2*l)/d
                w1=(a/(2*b))*(2-b)
                w2=0
                for(i in 1:n) {
                        harm2=(((2*b*a)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l)
                        w2=w2+harm2
                }
                return(w1+w2)
        }

        
        ##
        ## ecgpoints function
        ## ==========
        ecgpoints <- function(x, li,
                              a_p, d_p, t_p,
                              a_q, d_q, t_q,
                              a_qrs, d_qrs, t_qrs,
                              a_s, d_s, t_s,
                              a_t, d_t,t_t,
                              a_u, d_u, t_u) {
                
                if (input$show_p) {pwav=type1wav(x,a_p,d_p,t_p,li)}
                else {pwav=type1wav(x,0,d_p,t_p,li)}
                
                if (input$show_q) {qwav=type2wav(x,a_q,d_q,t_q,li)
                qrswav=type2wav(x,a_qrs,d_qrs,t_q, li)
                swav=type2wav(x,a_s,d_s,t_s,li)}
                else {qwav=type2wav(x,0,d_q,t_q,li)
                qrswav=type2wav(x,0,d_qrs,t_qrs, li)
                swav=type2wav(x,0,d_s,t_s,li)}
                
                if (input$show_t) {twav=type1wav(x,a_t,d_t,t_t,li)}
                else {twav=type1wav(x,0,d_t,t_t,li)}
                
                if (input$show_u) {uwav=type1wav(x,a_u,d_u,t_u,li)}
                else {uwav=type1wav(x,0,d_u,t_u,li)}
                
                ecg=pwav+qrswav+twav+swav+qwav+uwav
                return(ecg)
                
        }
        
        # multiplot function (used to combine (twelve) plots to create
        # a 12-lead ECG)
        # ==========
        
        multiplot <- function(..., plotlist = NULL, file, cols = 4, layout = NULL) {
                require(grid)
                plots <- c(list(...), plotlist)
                numPlots = length(plots)
                if (is.null(layout)) {
                        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                         ncol = cols, nrow = ceiling(numPlots/cols))
                }
                if (numPlots == 1) {
                        print(plots[[1]])
                } else {
                        grid.newpage()
                        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                        for (i in 1:numPlots) {
                                matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                                print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                                layout.pos.col = matchidx$col))
                }}}

        
        
        
        ## generate the data when inputs change
        ## (reactive process)
        
        ## =======================
        ## REACTIVE FUNCTION: data_shortlead
        ## =======================
        
        axis <- reactive({

                # calculate the ECG plots with initial values
                #
                # this for loop corrects the waves based on the axis values,
                # and feeds it back to the corrected_axis
                #
                # FUTURE USE: change the axis inputs for theaxisfactor
                temp_axis <- default_axis
                temp_axis["a_qrs",] <- temp_axis["a_qrs",]*theaxisfactor(input$qaxis, input$qaxis)
                temp_axis["a_p",] <- temp_axis["a_p",]*theaxisfactor(input$paxis, input$paxis)
                temp_axis["a_t",] <- temp_axis["a_t",]*theaxisfactor(input$taxis, input$taxis)
                temp_axis
        })
        
        temp_min_y <- function(lead) {
                floor(mean(ecg_data()[,lead])) - strip_y/2
        }
        
        temp_max_y <- function(lead) {
                temp_min_y(lead) + strip_y
        }
        
        zl1 <- function() {
                ggplot(ecg_data(), aes(x=x,y=I)) + geom_line() + ecgtheme +
                        annotate("text", label="I", x=0.5, y=temp_max_y(2), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(2)-0.5,temp_max_y(2)+2.5,small_h),
                                           breaks=seq(temp_min_y(2)-0.5,temp_max_y(2)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(2),temp_max_y(2)))
        }
        
        zl2 <- function() {
                ggplot(ecg_data(), aes(x=x,y=II)) + geom_line() + ecgtheme +
                        annotate("text", label="II", x=0.5, y=temp_max_y(3), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(3)-0.5,temp_max_y(3)+2.5,small_h),
                                           breaks=seq(temp_min_y(3)-0.5,temp_max_y(3)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(3),temp_max_y(3)))
        }
        
        zll2 <- function() {
                ggplot(ecg_data_long(), aes(x=x2,y=II)) + geom_line() + ecgtheme +
                        annotate("text", label="II", x=0.5, y=temp_max_y(3), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x*4,small_w),
                                           breaks=seq(0,strip_x*4, small_w*big_w),
                                           limits=c(0,strip_x*4)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(3)-0.5,temp_max_y(3)+2.5,small_h),
                                           breaks=seq(temp_min_y(3)-0.5,temp_max_y(3)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(3),temp_max_y(3)))
        }
        
        
        zl3 <- function() {
                ggplot(ecg_data(), aes(x=x,y=III)) + geom_line() + ecgtheme +
                        annotate("text", label="III", x=0.5, y=temp_max_y(4), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(4)-0.5,temp_max_y(4)+2.5,small_h),
                                           breaks=seq(temp_min_y(4)-0.5,temp_max_y(4)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(4),temp_max_y(4)))
        }
        
        zvr <- function() {
                ggplot(ecg_data(), aes(x=x,y=aVR)) + geom_line() + ecgtheme +
                        annotate("text", label="aVR", x=0.5, y=temp_max_y(5), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(5)-0.5,temp_max_y(5)+2.5,small_h),
                                           breaks=seq(temp_min_y(5)-0.5,temp_max_y(5)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(5),temp_max_y(5)))
        }
        
        zvl <- function() {
                ggplot(ecg_data(), aes(x=x,y=aVL)) + geom_line() + ecgtheme +
                        annotate("text", label="aVL", x=0.5, y=temp_max_y(6), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(6)-0.5,temp_max_y(6)+2.5,small_h),
                                           breaks=seq(temp_min_y(6)-0.5,temp_max_y(6)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(6),temp_max_y(6)))
        }
        
        zvf <- function() {
                ggplot(ecg_data(), aes(x=x,y=aVF)) + geom_line() + ecgtheme +
                        annotate("text", label="aVF", x=0.5, y=temp_max_y(7), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(7)-0.5,temp_max_y(7)+2.5,small_h),
                                           breaks=seq(temp_min_y(7)-0.5,temp_max_y(7)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(7),temp_max_y(7)))
        }

        
        zv1 <- function() {
                ggplot(ecg_data(), aes(x=x,y=V1)) + geom_line() + ecgtheme +
                        annotate("text", label="V1", x=0.5, y=temp_max_y(8), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(8)-0.5,temp_max_y(8)+2.5,small_h),
                                           breaks=seq(temp_min_y(8)-0.5,temp_max_y(8)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(8),temp_max_y(8)))
        } 
        
        zv2 <- function() {
                ggplot(ecg_data(), aes(x=x,y=V2)) + geom_line() + ecgtheme +
                        annotate("text", label="V2", x=0.5, y=temp_max_y(9), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(9)-0.5,temp_max_y(9)+2.5,small_h),
                                           breaks=seq(temp_min_y(9)-0.5,temp_max_y(9)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(9),temp_max_y(9)))
        }
        
        zv3 <- function() {
                ggplot(ecg_data(), aes(x=x,y=V3)) + geom_line() + ecgtheme +
                        annotate("text", label="V3", x=0.5, y=temp_max_y(10), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(10)-0.5,temp_max_y(10)+2.5,small_h),
                                           breaks=seq(temp_min_y(10)-0.5,temp_max_y(10)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(10),temp_max_y(10)))
        }
        
        zv4 <- function() {
                ggplot(ecg_data(), aes(x=x,y=V4)) + geom_line() + ecgtheme +
                        annotate("text", label="V4", x=0.5, y=temp_max_y(11), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(11)-0.5,temp_max_y(11)+2.5,small_h),
                                           breaks=seq(temp_min_y(11)-0.5,temp_max_y(11)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(11),temp_max_y(11)))
        }
        
        zv5 <- function() {
                ggplot(ecg_data(), aes(x=x,y=V5)) + geom_line() + ecgtheme +
                        annotate("text", label="V5", x=0.5, y=temp_max_y(12), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(12)-0.5,temp_max_y(12)+2.5,small_h),
                                           breaks=seq(temp_min_y(12)-0.5,temp_max_y(12)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(12),temp_max_y(12)))
        }
        
        
        zv6 <- function() {
                ggplot(ecg_data(), aes(x=x,y=V6)) + geom_line() + ecgtheme +
                        annotate("text", label="V6", x=0.5, y=temp_max_y(13), size = 8, colour="black") +
                        scale_x_continuous(minor_breaks=seq(0,strip_x,small_w),
                                           breaks=seq(0,strip_x, small_w*big_w),
                                           limits=c(0,strip_x)) + 
                        scale_y_continuous(minor_breaks=seq(temp_min_y(13)-0.5,temp_max_y(13)+2.5,small_h),
                                           breaks=seq(temp_min_y(13)-0.5,temp_max_y(13)+2.5,small_h*big_h),
                                           limits=c(temp_min_y(13),temp_max_y(13)))
        }
        
        ecg_data_long <- reactive({
                li <- 30/input$rate
                temp_ecg_data <- data.frame(x2)
                temp_ecg_data <- cbind(temp_ecg_data, ecgpoints(x2, li, 
                          axis()["a_p",2],axis()["d_p",2],axis()["t_p",2],
                          axis()["a_q",2],axis()["d_q",2],axis()["t_q",2],
                          axis()["a_qrs",2],axis()["d_qrs",2],axis()["t_qrs",2],
                          axis()["a_s",2],axis()["d_s",2],axis()["t_s",2],
                          axis()["a_t",2],axis()["d_t",2],axis()["t_t",2],
                          axis()["a_u",2],axis()["d_u",2],axis()["t_u",2]
                ))
                names(temp_ecg_data) <- c("x2","II")
                temp_ecg_data
        }) 
                
        ecg_data <- reactive({
                li <- 30/input$rate             ## parameters correct for HR
                temp_ecg_data <- data.frame(x)
                for (elead in 1:12) {
                        temp_ecg_data <- cbind(temp_ecg_data,
                        ecgpoints(x, li, 
                                axis()["a_p",elead],axis()["d_p",elead],axis()["t_p",elead],
                                axis()["a_q",elead],axis()["d_q",elead],axis()["t_q",elead],
                                axis()["a_qrs",elead],axis()["d_qrs",elead],axis()["t_qrs",elead],
                                axis()["a_s",elead],axis()["d_s",elead],axis()["t_s",elead],
                                axis()["a_t",elead],axis()["d_t",elead],axis()["t_t",elead],
                                axis()["a_u",elead],axis()["d_u",elead],axis()["t_u",elead]
                        ))}
                names(temp_ecg_data) <- c("x", "I","II","III","aVR","aVL","aVF",
                                          "V1","V2","V3","V4","V5","V6")
                temp_ecg_data
        })
        
        
        output$display12leads <- renderPlot({
             
        # draw the 12-lead ECG with the specified variables
        multiplot(zl1(), zvr(), zv1(), zv4(), zl2(), zvl(), zv2(), zv5(), zl3(), zvf(), zv3(), zv6(), zll2(), 
                  layout=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,13),
                                ncol=4,byrow=TRUE))
        })
        
        output$display_lead_I <- renderPlot({ zl1() })
        output$display_lead_II <- renderPlot({ zl2() })
        output$display_lead_III <- renderPlot({ zl3() })
        output$display_lead_aVR <- renderPlot({ zvr() })
        output$display_lead_aVL <- renderPlot({ zvl() })
        output$display_lead_aVF <- renderPlot({ zvf() })
        output$display_lead_V1 <- renderPlot({ zv1() })
        output$display_lead_V2 <- renderPlot({ zv2() })
        output$display_lead_V3 <- renderPlot({ zv3() })
        output$display_lead_V4 <- renderPlot({ zv4() })
        output$display_lead_V5 <- renderPlot({ zv5() })
        output$display_lead_V6 <- renderPlot({ zv6() })
##        output$display_debug <- renderText({ paste() })

})
