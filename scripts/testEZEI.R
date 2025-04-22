data("pt01EcoG")
timeWindow <- c(-10, 20)
epoch <- Epoch(pt01EcoG)
fs=1000
sozIndex <- attr(pt01EcoG, "sozIndex")
timeNum <- ncol(epoch)
windowParams = c(1, 0.2) 

epoch <- Epoch(pt01EcoG)
visuIEEGData(epoch)

EIPt01<-computeEpileptogenicIndex(epoch, windowParams)

ermaster<-EIPt01$energyRatio
  

plotER<-plotERHeatmap(ermaster,sozIndex=sozIndex)
plotER

#################

data("pt01EcoG")

## sozIndex is the index of the electrodes we assume are in the SOZ
sozIndex <- attr(pt01EcoG, "sozIndex")
## precomputed Epileptogenic Index object
data("pt01EI")

## plot the mean power heatmap
plotER<-plotERHeatmap(ei=pt01EI,sozIndex=sozIndex)
plotER<-plotER+ggplot2::ggtitle("Energy ratio heatmap for patient pt01")

plotER






