library(shiny)
library(zoo)

load("Data/year_res_max.RData")
obs.x = coredata(year.res.max)/100
num.x <- length(obs.x)
sf <- seq(1,1/num.x,by=-1/num.x)
order.x <- order(obs.x)

source("Data/return_level_plot.R")
source("Data/plot_sf.R")

sea_data = read.csv("Data/Sea_level_proj_anomalies.csv")
#sea_data = sea_data$x

storm_data = read.csv("Data/Storm_surge_freq.csv")
frequency = storm_data$V1
flood_meters = storm_data$V2
remove(storm_data)

# Find which q is the 100-yr value
my_num = which(frequency >= 0.01)
my_num.max = which.max(my_num)
year100prob <- my_num.max +1
# check to make sure the year100prob value is the 100-yr value (10^-2)
#print(frequency[year100prob])

stormFreqAcc2100 = read.csv("Data/StormFreqAccountUNC2100.csv")
returnperiods2100 <- stormFreqAcc2100[,2]
returnlevels2100 <- stormFreqAcc2100[,3]
return100_yr2100 <- stormFreqAcc2100[1,4]

stormFreqAcc2080 = read.csv("Data/StormFreqAccountUNC2080.csv")
returnperiods2080 <- stormFreqAcc2080[,2]
returnlevels2080 <- stormFreqAcc2080[,3]
return100_yr2080 <- stormFreqAcc2080[1,4]

stormFreqAcc2050 = read.csv("Data/StormFreqAccountUNC2050.csv")
returnperiods2050 <- stormFreqAcc2050[,2]
returnlevels2050 <- stormFreqAcc2050[,3]
return100_yr2050 <- stormFreqAcc2050[1,4]

stormFreqAcc2020 = read.csv("Data/StormFreqAccountUNC2020.csv")
returnperiods2020 <- stormFreqAcc2020[,2]
returnlevels2020 <- stormFreqAcc2020[,3]
return100_yr2020 <- stormFreqAcc2020[1,4]

med_x = mean(sea_data$X2050)
#min_x = min(sea_data$X2020)
max_x = max(sea_data$X2100)

MinMaxsf <- mat.or.vec(length(flood_meters), 2)
MinMaxsf[1, 1:2] = flood_meters[2]

new_line <- rep(NA, length(flood_meters))
new_line[1] = flood_meters[2]

slr.plotted <- rep(NA, length(flood_meters))
slr.plotted[1] = flood_meters[2]

ui <- fluidPage(
  # Application title
      titlePanel(
        title=h1('Future Storm Frequency Tool', style="position: center")),
  fixedRow(
    # MAIN space
    column(12,
  
  fluidRow(
    column(12,
           tabsetPanel(
             tabPanel(p("Welcome", style = "color: #000000"), value=1), 
             tabPanel(p("Sea-level rise", style = "color: #000000"), value=2),
             tabPanel(p("Sea level cumulative density", style = "color: #000000"), value=3),
             tabPanel(p("Storm surge", style = "color: #000000"), value=4),
             id = "conditionedPanels"
             ))),
    fluidRow(
      column(12, 
             conditionalPanel(condition="input.conditionedPanels==1",
                              h3('Welcome to the online Future Storm Frequency Tool!', style = "padding: 0px 100px 10px"),
                              br(),
                              #h4('Note'),
                              p(em("The information presented in this app. is a representation of the results in Ruckert et al. (in prep). This information is 
privilaged and hence not intended to be public. Once the manuscript is published, 
                                 it will of course be open to the public."), style = "padding: 0px 100px 10px"),
                              br(),
                              h4('Basic information', style = "padding: 0px 100px 10px"),
                              p('"The warming climate is causing sea-levels to rise around the globe. As sea-levels rise, settlements and 
                                ecosystems in low-lying coastal area become more vulnerable to flooding" (Ruckert et al. (in prep.)). In 
                                the United States, strategies and infrastructures that manage flood risk are often designed for the annual 
                                flooding probability of 1 in 100 (also know as the 100-yr event) (Ruckert et al. (in prep.)).', style = "padding: 0px 100px 10px"),
                              p('The FSF tool is a tool to predict future 100-yr storm surge based on sea-level rise. The online tool 
                                corresponds to the San Francisco Bay area and is a pilot version of the results in Ruckert et al. (in prep.). 
                                This tool illustrates a simple point and should', strong(' NOT'),  ' be used to make decisions.', style = "padding: 0px 100px 10px"),
                              p('This tool illustrates that increasing sea-level increases the return level and reduces the flood frequency. 
                                Additionally, accounting for sea-level rise uncertainty increases the return level and reduces 
                                the flood frequency in comparison to using the mean sea-level projection. 
                                More details can be found in Ruckert et al. (in prep.).', style = "padding: 0px 100px 10px"),
                              br(),
                              h4('Key features', style = "padding: 0px 100px 10px"),
                              p('In this online version, we give the user the ability to manipulate sea-level to determine future storm surge 
                                frequency. Our online version includes the ability to:', style = "padding: 0px 100px 10px"),
                              tags$ol(
                                tags$li("Select a projection year (i.e., 2020, 2050, 2080, and 2100)"), 
                                tags$li("Increase sea-level by 0 to 3.0 m"), 
                                tags$li("Quantify the 100-yr storm surge based on the chosen sea-level anomaly"),
                                tags$li("Add lines for comparision such as the mean and mean ± 1 standard deviation"), 
                                tags$li("Compare the generated storm surge frequency to the frequency generated by accounting for probabilistic
                                        sea-level rise."), style = "padding: 0px 100px 10px"
                                ),
                              br(),
                              h4('Other details', style = "padding: 0px 100px 10px"),
                              p('Based on our analysis, the current the 100-yr storm surge for San Francisco Bay is 1.59 m. In the survival function plot, the current 
                                storm surge is represented as the light blue line and the black dots are the annual block maxima observations 
                                (the largest storm surge recorded in a particular year). Also the gray envelope represents the 99% credible interval of 
                                storm surge/sea-level rise for the selected year.', style = "padding: 0px 100px 10px"),
                              p('The sea-level projections are with respect to the year 2000 and are global estimates. Additionally, we increase sea-level rise by 4, 
                                21, 40, and 55 cm to account for changes in dams and reserviors in the selected years (Heberger et al. 2009).', style = "padding: 0px 100px 10px"),
                              br(),
                              h4('References', style = "padding: 0px 100px 10px"),
                              tags$ul(
                                tags$li('Heberger M., Cooley H., Herrera P., Gleick P., Moore E. 2009. The impacts of sea-level rise on the California Coast.'), 
                                tags$li('Ruckert K.L., Oddo P.C., and Keller K. Accounting for sea-level rise uncertainty increases flood risk area: An example from San Francisco Bay. (In prep.).'), 
                                style = "padding: 0px 100px 10px"),
                              # p(, style = "padding: 0px 100px 10px"),
                              # p(
                              #   ', style = "padding: 0px 100px 10px"),
                              br(),
                              h5('Pilot developed by Kelsey Ruckert. 
                                 Penn State University: klr324@psu.edu.')
                              )
             )
      ),
    fluidRow(
      column(8, 
             conditionalPanel(condition="input.conditionedPanels==2",
                              plotOutput("dens", height=500)),
             conditionalPanel(condition="input.conditionedPanels==3",
                              plotOutput("cdf", height=500)),
             conditionalPanel(condition="input.conditionedPanels==4",
                              plotOutput("survf", height=500))
      ),
      column(4, wellPanel(
        conditionalPanel(condition="input.conditionedPanels==1",
                         h4("Find out more", style="position: center"), 
                         a(href="http://www.scrimhub.org", "Click here for more information on research conducted at SCRiM (the Network for Sustainable Climate Risk 
                           Management).")
        ),
        conditionalPanel(condition="input.conditionedPanels!=1",
                         h4('Sea-level Rise in San Francisco Bay'),
                         p('Current 100-yr storm surge is 1.59 m.'), 
                         p('Gray envelope represents the 99% credible interval.'),
                         p('The sea-level anomaly and the associated storm surge are printed in the plot title of the Storm surge tab'),
                         
                         br(),
                         
                         radioButtons("year", "Projection year:",
                                      c("2020" = "Year_2020",
                                        "2050" = "Year_2050",
                                        "2080" = "Year_2080",
                                        "2100" = "Year_2100")),
                         br(),
                         
                         sliderInput(inputId = "num", 
                                     label = "Sea-level anomaly (m):", 
                                     value = round(med_x,2), min = 0, max = round(max_x,2), step = 0.05),
                         
                         br(),
                         p(strong('Add informative lines:')),
                         checkboxInput(inputId = "my_mean",
                                       label = strong("Mean sea level + storm surge"),
                                       value = FALSE),
                         
                         checkboxInput(inputId = "p_std",
                                       label = "Mean +1 standard deviation",
                                       value = FALSE),
                         
                         checkboxInput(inputId = "m_std",
                                       label = "Mean -1 standard deviation",
                                       value = FALSE),
                         
                         checkboxInput(inputId = "account_unc",
                                       label = strong("Account for probabilistic sea-level"),
                                       value = FALSE),
                         
                         h6('Pilot developed in Spring 2016 by Kelsey Ruckert. 
                            Penn State University: klr324@psu.edu.')
                         )
        )
      )
    )
  )
  )
)

server <- function(input, output) {
  
  data <- reactive({
    year <- switch(input$year,
                   Year_2020 = sea_data$X2020,
                   Year_2050 = sea_data$X2050,
                   Year_2080 = sea_data$X2080,
                   Year_2100 = sea_data$X2100,
                   sea_data$X2100)

  })

  output$dens <- renderPlot({
    den_x = density(data())

    par(mgp=c(1.5,.5,0),mar=c(3.5, 3, 3.5, 2))
    plot(den_x, col="gray", main="",
         lwd=3, xlab="",
    sub="With respect to the year 2000",
    ylab="", yaxt="n")
    title(xlab="Sea-level anomaly (m)", ylab="Probability density", cex.lab=1.5)
    abline(v = input$num, lwd = 4, col="royalblue")
    
    if (input$my_mean)  {
      abline(v = mean(data()), lwd = 4, col="mediumseagreen")
    } else {
    }
    
    if (input$p_std)  {
      abline(v = mean(data()) + sd(data()), lwd = 4, col="darkgreen")
    } else {
    }
    
    if (input$m_std)  {
      abline(v = mean(data()) - sd(data()), lwd = 4, col="darkseagreen1")
    } else {
    }

  })

  output$cdf <- renderPlot({
    cdf_x = ecdf(data())
    par(mgp=c(1.5,.5,0),mar=c(3.5, 3, 3.5, 2))
    plot(cdf_x, lwd=3, ylab="", xlab="", col="gray",
         sub="With respect to the year 2000", main="")
    title(xlab="Sea-level anomaly (m)", ylab="Cumulative density", cex.lab=1.5)

    abline(h=0.5, lty=2) # add in the 1 in 100

    abline(v = input$num, lwd = 4, col="royalblue")
    
    if (input$my_mean)  {
      abline(v = mean(data()), lwd = 4, col="mediumseagreen")
    } else {
    }
    
    if (input$p_std)  {
      abline(v = mean(data()) + sd(data()), lwd = 4, col="darkgreen")
    } else {
    }
    
    if (input$m_std)  {
      abline(v = mean(data()) - sd(data()), lwd = 4, col="darkseagreen1")
    } else {
    }
    
  })

  output$survf <- renderPlot({

    quant_point5 = quantile(data(),0.005) #99% credible interval
    quant_995 = quantile(data(),0.995)

    MinMaxsf[2:length(flood_meters),1] <- flood_meters[2:length(flood_meters)] + quant_point5
    MinMaxsf[2:length(flood_meters),2] <- flood_meters[2:length(flood_meters)] + quant_995
    freq_x_axis=c(MinMaxsf[,1], rev(MinMaxsf[,2])); freq_y_axis=c(frequency, rev(frequency))

    slr.plotted[2:length(flood_meters)] <- flood_meters[2:length(flood_meters)] + input$num

    par(mgp=c(1.5,.5,0),mar=c(3.5, 3, 3.5, 2))
    plot.sf(coredata(year.res.max)/100, pch = 20, bg = "white",
            ylab = "",
            xlab = "",
            yaxt = "n", yaxs = 'i',
            ylim = c(10^-2.5, 10^0+0.25),
            xlim = c(flood_meters[2], MinMaxsf[10000,2]+0.25),
            main = paste('Sea-level rise =',' ', input$num,' m ', '& 1:100 design level =', ' ', round(slr.plotted[year100prob],2),' ','m', sep=''))
    title(xlab="Return Level (m) in San Francisco Bay", ylab="Survival function [1 - cdf]", cex.lab=1.5)

    polygon(freq_x_axis, freq_y_axis, col="gray")
    points(obs.x[order.x], sf, pch=21, bg="black")

    lines(flood_meters, frequency, type="l",lwd=4, col = "skyblue")
    lines(slr.plotted, frequency, type="l",col = "royalblue",lwd=4)
    
    if (input$my_mean)  {
      new_line[2:length(flood_meters)] = flood_meters[2:length(flood_meters)] + mean(data())
      lines(new_line, frequency, type="l",lwd=4, col = "mediumseagreen")
      points(new_line[year100prob], frequency[year100prob], pch=21, bg="mediumseagreen")
      text(new_line[year100prob]-0.05, 0.007, paste(round(new_line[year100prob],2),' ','m', sep=''), cex=1.1)
    } else {
    }
    
    if (input$p_std)  {
      new_line[2:length(flood_meters)] = flood_meters[2:length(flood_meters)] + (mean(data()) + sd(data()))
      lines(new_line, frequency, type="l",lwd=4, col = "darkgreen")
    } else {
    }
    
    if (input$m_std)  {
      new_line[2:length(flood_meters)] = flood_meters[2:length(flood_meters)] + (mean(data()) - sd(data()))
      lines(new_line, frequency, type="l",lwd=4, col = "darkseagreen1")
    } else {
    }

    #Add in 100-yr value lines
    abline(h=0.01, lty=2) # add in the 1 in 100
    axis(2, at=10^(-4:-2), label=parse(text=paste("10^", -4:-2, sep="")))
    text(1.25, 0.013, "1:100 level", cex=1.1)

    # Plot points at 1:100 value
    points(flood_meters[year100prob], frequency[year100prob], pch=21, bg="skyblue")
    points(slr.plotted[year100prob], frequency[year100prob], pch=21, bg="royalblue")

    if (input$account_unc && input$year == "Year_2020")  {
      lines(returnlevels2020, returnperiods2020, type="l", col="slateblue", lwd=3.5)
      points(returnlevels2020[return100_yr2020], returnperiods2020[return100_yr2020], pch=21, bg="slateblue")
      text(1.75, 0.013, paste(round(returnlevels2020[return100_yr2020],2),' ','m', sep=''), cex=1.1)

    } else if (input$account_unc && input$year == "Year_2050") {
      lines(returnlevels2050, returnperiods2050, type="l", col="slateblue", lwd=3.5)
      points(returnlevels2050[return100_yr2050], returnperiods2050[return100_yr2050], pch=21, bg="slateblue")
      text(2.15, 0.013, paste(round(returnlevels2050[return100_yr2050],2),' ','m', sep=''), cex=1.1)

    } else if (input$account_unc && input$year == "Year_2080") {
      lines(returnlevels2080, returnperiods2080, type="l", col="slateblue", lwd=3.5)
      points(returnlevels2080[return100_yr2080], returnperiods2080[return100_yr2080], pch=21, bg="slateblue")
      text(2.85, 0.013, paste(round(returnlevels2080[return100_yr2080],2),' ','m', sep=''), cex=1.1)

    } else if (input$account_unc && input$year == "Year_2100") {
   lines(returnlevels2100, returnperiods2100, type="l", col="slateblue", lwd=3.5)
   points(returnlevels2100[return100_yr2100], returnperiods2100[return100_yr2100], pch=21, bg="slateblue")
   text(3.5, 0.013, paste(round(returnlevels2100[return100_yr2100],2),' ','m', sep=''), cex=1.1)
    } else {
}

  })
  
}

shinyApp(ui = ui, server = server)