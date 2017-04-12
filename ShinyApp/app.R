#########################################################################
# Copyright 2016 The Pennsylvania State University
#
# This app was developed by Kelsey Ruckert to showcase the results in 
# Ruckert K.L., Oddo P.C., and Keller K. Accounting for sea-level rise uncertainty increases flood risk 
# area: An example from San Francisco Bay. (In prep.).  This work was supported by the National
# Science Foundation through the Network for Sustainable Climate Risk
# Management (SCRiM) under NSF cooperative agreement GEO-1240507.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#########################################################################
library("shiny")
library("shinyRGL")
library(plotrix)
library(zoo)
library(ggvis)
library(RColorBrewer)

mycolors = brewer.pal(6,"RdYlBu")

# Load functions
#source("Data/return_level_plot.R")
source("Data/plot_sf.r")

# Load data
source("Data/Input_preanalysed_data.R")
colors_grays = colorRampPalette(c("gray20", "gray100"))

colors_blues = colorRampPalette(c("steelblue", "aliceblue"))

ui <- fluidPage(
  tags$style(type = "text/css", "label { font-size: 16px; align:center}"),
  br(),
  fluidRow(
    column(12,
           h1("Future Flood Height Probability Tutorial"),
           h3('Introduction', style="font-style: italic"),
           p("Predicting future flood height (storm surge including sea-level rise) for a particular area is a difficult task for decision-makers for several reasons. 
             First, not all regions have extensive tide data. Second, future changes in sea level are deeply uncertain, especially 
             at the local level, so the choice of how to represent future sea-level rise is a difficult, yet important issue. For 
             example, overestimating future sea-level rise could lead to unnecessary high costs for protection measures, whereas 
             underestimating future sea-level rise could put infrastructure and lives at risk."),
           h3('Objective', style="font-style: italic"),
           p("This interactive tutorial illustrates how changes in sea level affect flooding events and their 
              associated probability of occurence using San Francisco, California 
              as an example. More importantly, this tutorial illustrates how accounting for uncertainty in sea-level 
              projections increases the return level (the value associated with a probability of occurrence) 
              compared to accounting for the mean sea-level projection."),
           h3('Directions', style="font-style: italic"),
           p("In this tutorial, you will determine (Step #1) what year to re-evaluate the flood probability in San Francisco by clicking on one of the button choices.  
             Next, in Step #2 you can choose how much sea-level rise to account for in that year by moving the slider. In Step #3, you can choose to 
             compare the future flood probability calculated from your choices to (i) nothing, (ii) the future flood probability accounting for the mean sea-level 
             anomaly, and to (iii) the future flood probability accounting for uncertainty in future sea-level rise. When you're satisfied with your choices, 
             click the 'Calculate flood height!' button. This action will cause the application to calculate the future flood probability using your choices in 
             Step #1 and #2, generate a flood probability plot, and generate a table for comparison to the choice made in Step #3. The 'Calculate flood height!' button, 
             only needs to be clicked on once. Afterwards, the flood probability plot and table will be modified automatically using the options in the three steps."),
           wellPanel(style = "background-color: #FFFFFF;",
           helpText("Please consider sending feedback to: ", tags$b("klr324@psu.edu."), " This helps us tremendously  
             with improving the tutorial and our future tutorials.",
             style="font-style: italic; text-align: center"),
           
           helpText('The source code is avaliable on ', a(href="https://github.com/scrim-network/Ruckertetal_SFB_SLR_2016", "Github", 
                                                                  target="_blank"), '  along with the code for the associated paper, ', 
                    a(href="http://dx.doi.org/10.1371%2Fjournal.pone.0174666", "Ruckert et al. (2017)", target="_blank"), '.',
           style="font-style: italic; text-align: center")
           ),
           hr()
           ),
           
           conditionalPanel(condition = "input.checkbox != 1",
                            column(8, offset = 4, 
                                   div(checkboxInput("checkbox", label = "Agree to the terms and conditions below", value = 0), style="text-align: center")),
                            
                            column(12,
                                   wellPanel(
             h3('GNU General Public License', style="font-style: italic; text-align: center"),
             p('Copyright 2016 Kelsey Ruckert. Penn State University: klr324@psu.edu.', style="font-style: italic; text-align: center"),
             p('This app was developed by Kelsey Ruckert and this work was supported 
               by the National Science Foundation through the Network for Sustainable Climate Risk 
               Management (SCRiM) under NSF cooperative agreement GEO-1240507.', style="font-style: italic; text-align: center"),
             p('Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
               and associated documentation files (the "Software"), to deal in the Software without restriction, 
               including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
               and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
               subject to the following conditions:', style="font-style: italic; text-align: center"),
             p('The above copyright notice and this permission notice shall be included in all copies or substantial 
               portions of the Software.', style="font-style: italic; text-align: center"),
             p('THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANITY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING 
               BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTIBILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
               NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
               DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
               OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.', style="font-style: italic; text-align: center")
                                   ),
             hr())
           )
    ),
  conditionalPanel(condition = "input.checkbox == 1",
  fluidRow(
    column(4,
           h3("Step #1: Choose a re-evalutation year"),
           helpText("In what year do you plan to re-evaluate the flood probability in San Francisco, California?"),
           radioButtons("year", "Projection year:",
                        c("2020" = "Year_2020",
                          "2050" = "Year_2050",
                          "2080" = "Year_2080",
                          "2100" = "Year_2100"),
                        selected = "Year_2100")
    ),
    column(4,
           h3("Step #2: Account for sea-level rise"),
           helpText("How much sea level do you want to take into consideration when re-evaluating flood probability?"),
           sliderInput(inputId = "num",
                       label = "Sea-level anomaly (m):",
                       value = 0, min = 0, max = 2.5, step = 0.01)
           ),
    column(4,
           h3("Step #3: Compare future flood height probability"),
           helpText("Do you want to compare your estimated future flood height (baseline storm surge + the accounted future sea-level rise) 
                    to the flood height accounting for the mean sea-level estimate or the flood height accounting 
                    for uncterain sea-level rise?"),
           radioButtons("add_line", "Compare:",
                        c("None" = "no_lines",
                          "Mean sea level" = "my_mean",
                          "Account for uncertain sea level" = "account_unc"),
                        selected = "no_lines"),
           actionButton("submitstep3", "Calculate flood height!")
           )

  ),
  fluidRow(
    column(12,
           hr()
    )
  ),

  fluidRow(
    column(6,
           plotOutput("dens", height=500),
           br(),
           p("Above is the probability density function of the projected sea-level rise in the year you specified. 
             The yellow verticle line displays the amount of sea-level rise you have chosen to account 
             for when re-evaluating flood risks and where this estimate is located along the probability density function. If you choose to compare 
             with accounting for the mean, then a blue verticle line will appear indicating the location of the mean. If you choose to compare with 
             accounting for uncertain sea-level rise, then the probability density function will be colored with a blue gradient indicating that each 
             sea-level anomaly is accounted for when considering future flood probability. It's represented as a gradient because each anomaly has a different 
             probability of occuring.")
    ),
    conditionalPanel(condition = "input.submitstep3 == 1",
    column(6,
           plotOutput("survf", height=500),
           br(),
           p("Above is the survival function displaying the flood height probability for San Francisco. The baseline flood probability is 
             shown as the light blue line, the historical observations are black dots, and the gray envelope is the 99% credible interval 
             of the baseline storm surge plus the future sea-level rise. The 99% credible interval is shown as a gradient from more likely 
             in gray to less likely in white. The yellow line is the user's predicted future flood height probability. This is estimated 
             by adding the sea-level anomaly specified in Step #2 to the baseline flood probability. The 
             blue line and the dark blue line is the flood height accounting for the mean sea-level and the flood height accounting 
             sea-level rise uncertainty, respectively.")
  ),
  fluidRow(
    column(12,
           hr()
    )
  ),
  fluidRow(
    column(12,
           p("Step #3 (b): The table compares the estimates calculated from the user's specifications and the results from accounting for the mean 
             sea-level or the results from accounting for sea-level rise uncertainty."),
           tableOutput("comparison.table"), align = 'center'
           )
  ))),
  fluidRow(
    column(12,
           hr()
    )
  )
  ),
  fluidRow(
    column(12,
           p('[Figures and text modified from ', a(href="http://dx.doi.org/10.1371%2Fjournal.pone.0174666", "Ruckert et al. (2017)", target="_blank"), '.]', 
             style="font-style: italic; text-align: center"),
           h6('Shiny app developed by Kelsey Ruckert. Penn State University: klr324@psu.edu.', style="font-style: italic; text-align: center"),
           helpText('Note: this illustrative tutorial is not to be used to assess actual coastal hazards.', style="font-style: italic; text-align: center"),
           hr()
    )
  )
)



server <- function(input, output, session) {

# # Sea-level plot
# #=================================
data <- reactive({
  year <- switch(input$year,
                 Year_2020 = sea_data$X2020,
                 Year_2050 = sea_data$X2050,
                 Year_2080 = sea_data$X2080,
                 Year_2100 = sea_data$X2100)

})

output$dens <- renderPlot({
  den_x = density(data())

  if (input$add_line == "my_mean")  {
    par(mgp=c(1.5,.5,0),mar=c(3.5, 3, 3.5, 2))
    plot(den_x, col = mycolors[6], main="",
         lwd=3, xlab="",
         sub="With respect to the year 2000",
         ylab="", yaxt="n")
    
    title(xlab="Sea-level anomaly (m)", ylab="Probability density", main="Step #1 & Step #2:", cex.lab=1.5)
    
    abline(v = mean(data()), lwd = 4, col = mycolors[5])
    
  } else if (input$add_line == "account_unc") {
    par(mgp=c(1.5,.5,0),mar=c(3.5, 3, 3.5, 2))
    plot(den_x, col = mycolors[6], main="",
         lwd=3, xlab="",
         sub="With respect to the year 2000",
         ylab="", yaxt="n")
    
    poly_x_axis=c(den_x$x, seq(max(den_x$x), min(den_x$x), length.out=length(den_x$x))); poly_y_axis=c(den_x$y, rep(100, length(den_x$y)))
    gradient.rect(min(den_x$x+0.1), 0, max(den_x$x), max(den_x$y),col=colors_blues(100), gradient="x",border=NA)
    polygon(poly_x_axis, poly_y_axis, col="white",border=NA)
    lines(den_x, col= mycolors[6], lwd=2)
    box()
    
    title(xlab="Sea-level anomaly (m)", ylab="Probability density", main="Step #1 & Step #2:", cex.lab=1.5)
    
  } else { 
    par(mgp=c(1.5,.5,0),mar=c(3.5, 3, 3.5, 2))
    plot(den_x, col = mycolors[6], main="",
         lwd=3, xlab="",
         sub="With respect to the year 2000",
         ylab="", yaxt="n")
    
    title(xlab="Sea-level anomaly (m)", ylab="Probability density", main="Step #1 & Step #2:", cex.lab=1.5)
    }
        abline(v = input$num, lwd = 4, col = mycolors[3])
  })

# Future Storm Surge plot
#=================================
output$survf <- renderPlot({

  if (input$submitstep3 == 0)
    return()

  quant_point5 = quantile(data(),0.005) #99% credible interval
  quant_995 = quantile(data(),0.995)

  MinMaxsf[2:length(flood_meters),1] <- flood_meters[2:length(flood_meters)] + quant_point5
  MinMaxsf[2:length(flood_meters),2] <- flood_meters[2:length(flood_meters)] + quant_995
  
  freq_x_axis=c(MinMaxsf[,1], rep(0, length(MinMaxsf[,1]))); freq_y_axis=c(frequency, rev(frequency))
  freq_x_axis_2=c(MinMaxsf[,2], rep(5, length(MinMaxsf[,1]))); freq_y_axis=c(frequency, rev(frequency))

  slr.plotted[2:length(flood_meters)] <- flood_meters[2:length(flood_meters)] + input$num

  par(mgp=c(1.5,.5,0),mar=c(3.5, 3, 3.5, 3))
  plot.sf(coredata(year.res.max)/100, pch = 20, bg = "white",
          ylab = "",
          xlab = "",
          yaxt = "n", yaxs = 'i',
          ylim = c(10^-2.5, 10^0+0.25),
          xlim = c(flood_meters[2], MinMaxsf[10000,2]+0.25))
  title(xlab="Return level (m) in San Francisco Bay", ylab="Survival function [1 - cumulative frequency]", main="Step #3 (a):", cex.lab=1.5)
  axis(4, at=c(10^0, 10^-1, 10^-2), label=c("1", "10", "100"), cex=1.5)
  mtext(4, text="Return period (years)", line=2, cex=1.5)
  
  if (input$year == "Year_2020")  {
    gradient.rect(1,10e-5,2,1,col=colors_grays(100),gradient="x",border=NA)
    
  } else if (input$year == "Year_2050") {
    gradient.rect(1,10e-5,2.2,1,col=colors_grays(100),gradient="x",border=NA)
    
  } else if (input$year == "Year_2080") {
    gradient.rect(1,10e-5,2.75,1,col=colors_grays(100),gradient="x",border=NA)
    
  } else {
    gradient.rect(1,10e-5,3.25,1,col=colors_grays(100),gradient="x",border=NA)
  }

  polygon(freq_x_axis, freq_y_axis, col="white",border=NA)
  polygon(freq_x_axis_2, freq_y_axis, col="white",border=NA)
  points(obs.x[order.x], sf, pch=21, bg="black")
  box()

  lines(flood_meters, frequency, type="l",lwd=4, col = mycolors[4])
  lines(slr.plotted, frequency, type="l",col = mycolors[3],lwd=4)
  points(slr.plotted[year100prob], frequency[year100prob], pch=21, bg=mycolors[3])
  text(slr.plotted[year100prob]-0.05, 0.007, paste(round(slr.plotted[year100prob],2),' ','m', sep=''), cex=1.1)

  if (input$add_line == "my_mean")  {
    new_line[2:length(flood_meters)] = flood_meters[2:length(flood_meters)] + mean(data())
    lines(new_line, frequency, type="l",lwd=4, col = mycolors[5])
    points(new_line[year100prob], frequency[year100prob], pch=21, bg=mycolors[5])
    text(new_line[year100prob]-0.05, 0.007, paste(round(new_line[year100prob],2),' ','m', sep=''), cex=1.1)
  } else {
  }

  #Add in 100-yr value lines
  abline(h=0.01, lty=2) # add in the 1 in 100
  axis(2, at=10^(-4:-2), label=parse(text=paste("10^", -4:-2, sep="")))
  text(1.25, 0.013, "100-yr
return period", cex=1.1)

  # Plot points at 1:100 value
  points(flood_meters[year100prob], frequency[year100prob], pch=21, bg=mycolors[4])
  points(slr.plotted[year100prob], frequency[year100prob], pch=21, bg=mycolors[3])

  if (input$add_line == "account_unc" && input$year == "Year_2020")  {
    lines(returnlevels2020, returnperiods2020, type="l", col=mycolors[6], lwd=3.5)
    points(returnlevels2020[return100_yr2020], returnperiods2020[return100_yr2020], pch=21, bg=mycolors[6])
    text(1.75, 0.013, paste(round(returnlevels2020[return100_yr2020],2),' ','m', sep=''), cex=1.1)

  } else if (input$add_line == "account_unc" && input$year == "Year_2050") {
    lines(returnlevels2050, returnperiods2050, type="l", col=mycolors[6], lwd=3.5)
    points(returnlevels2050[return100_yr2050], returnperiods2050[return100_yr2050], pch=21, bg=mycolors[6])
    text(2, 0.013, paste(round(returnlevels2050[return100_yr2050],2),' ','m', sep=''), cex=1.1)

  } else if (input$add_line == "account_unc" && input$year == "Year_2080") {
    lines(returnlevels2080, returnperiods2080, type="l", col=mycolors[6], lwd=3.5)
    points(returnlevels2080[return100_yr2080], returnperiods2080[return100_yr2080], pch=21, bg=mycolors[6])
    text(2.5, 0.013, paste(round(returnlevels2080[return100_yr2080],2),' ','m', sep=''), cex=1.1)

  } else if (input$add_line == "account_unc" && input$year == "Year_2100") {
    lines(returnlevels2100, returnperiods2100, type="l", col=mycolors[6], lwd=3.5)
    points(returnlevels2100[return100_yr2100], returnperiods2100[return100_yr2100], pch=21, bg=mycolors[6])
    text(3, 0.013, paste(round(returnlevels2100[return100_yr2100],2),' ','m', sep=''), cex=1.1)
  } else {
  }
#   Table
# #================
output$comparison.table <- renderTable(rownames = TRUE, {

  user = c(round(flood_meters[year100prob],2), input$num, round(slr.plotted[year100prob],2), "-", "-")

  if (input$add_line == "my_mean")  {
    increase = new_line[year100prob] - slr.plotted[year100prob]
    percent_increase = round(increase/slr.plotted[year100prob],2)*100
    compare = c(round(flood_meters[year100prob],2), round(mean(data()),2), round(new_line[year100prob],2),
                round(increase,2), percent_increase)

    m <- matrix(c(user[1], compare[1], user[2], compare[2], user[3], compare[3], user[4], compare[4], user[5], compare[5]), nrow=2, ncol=5)
    rownames(m) <- c("User", "Mean")
    colnames(m) <- c("Baseline 100-yr storm surge (m)", "Sea-level anomaly (m)", "Estimated 100-yr return level (m)",
                     "Increase over user estimated return level (m)", "Increase over user estimated return level (%)")
    m

  } else if (input$add_line == "account_unc" && input$year == "Year_2020")  {
    increase = returnlevels2020[return100_yr2020] - slr.plotted[year100prob]
    percent_increase = round(increase/slr.plotted[year100prob],2)*100
    compare = c(round(flood_meters[year100prob],2), "-", round(returnlevels2020[return100_yr2020],2),
                round(increase,2), percent_increase)

    m <- matrix(c(user[1], compare[1], user[2], compare[2], user[3], compare[3], user[4], compare[4], user[5], compare[5]), nrow=2, ncol=5)
    rownames(m) <- c("User", "Uncertainty")
    colnames(m) <- c("Baseline 100-yr storm surge (m)", "Sea-level anomaly (m)", "Estimated 100-yr return level (m)",
                     "Increase over user estimated return level (m)", "Increase over user estimated return level (%)")
    m

  } else if (input$add_line == "account_unc" && input$year == "Year_2050")  {
    increase = returnlevels2050[return100_yr2050] - slr.plotted[year100prob]
    percent_increase = round(increase/slr.plotted[year100prob],2)*100
    compare = c(round(flood_meters[year100prob],2), "-", round(returnlevels2050[return100_yr2050],2),
                round(increase,2), percent_increase)

    m <- matrix(c(user[1], compare[1], user[2], compare[2], user[3], compare[3], user[4], compare[4], user[5], compare[5]), nrow=2, ncol=5)
    rownames(m) <- c("User", "Uncertainty")
    colnames(m) <- c("Baseline 100-yr storm surge (m)", "Sea-level anomaly (m)", "Estimated 100-yr return level (m)",
                     "Increase over user estimated return level (m)", "Increase over user estimated return level (%)")
    m

  } else if (input$add_line == "account_unc" && input$year == "Year_2080")  {
    increase = returnlevels2080[return100_yr2080] - slr.plotted[year100prob]
    percent_increase = round(increase/slr.plotted[year100prob],2)*100
    compare = c(round(flood_meters[year100prob],2), "-", round(returnlevels2080[return100_yr2080],2),
                round(increase,2), percent_increase)

    m <- matrix(c(user[1], compare[1], user[2], compare[2], user[3], compare[3], user[4], compare[4], user[5], compare[5]), nrow=2, ncol=5)
    rownames(m) <- c("User", "Uncertainty")
    colnames(m) <- c("Baseline 100-yr storm surge (m)", "Sea-level anomaly (m)", "Estimated 100-yr return level (m)",
                     "Increase over user estimated return level (m)", "Increase over user estimated return level (%)")
    m

  } else if (input$add_line == "account_unc" && input$year == "Year_2100")  {
    increase = returnlevels2100[return100_yr2100] - slr.plotted[year100prob]
    percent_increase = round(increase/slr.plotted[year100prob],2)*100
    compare = c(round(flood_meters[year100prob],2), "-", round(returnlevels2100[return100_yr2100],2),
                round(increase,2), percent_increase)

    m <- matrix(c(user[1], compare[1], user[2], compare[2], user[3], compare[3], user[4], compare[4], user[5], compare[5]), nrow=2, ncol=5)
    rownames(m) <- c("User", "Uncertainty")
    colnames(m) <- c("Baseline 100-yr storm surge (m)", "Sea-level anomaly (m)", "Estimated 100-yr return level (m)",
                     "Increase over user estimated return level (m)", "Increase over user estimated return level (%)")
    m
  } else {
    m <- matrix(c(round(flood_meters[year100prob],2), input$num, round(slr.plotted[year100prob],2)), nrow=1, ncol=3)
    rownames(m) <- c("User")
    colnames(m) <- c("Baseline 100-yr storm surge (m)", "Sea-level anomaly (m)", "Estimated 100-yr return level (m)")
    m
  }
  })
  })
##########################################
# Must accept terms and conditions before accessing app #
##########################################

# # we need to have a quasi-variable flag to indicate when license and conditions are accepted.
output$input.checkbox <- reactive({FALSE})
outputOptions(output, 'input.checkbox', suspendWhenHidden = FALSE)

output$input.checkbox <- reactive({ TRUE })

}

shinyApp(ui = ui, server = server)
