# test EZEI  with Epoch package
# data download one patient one seizure

library(EZEI)
library(Epoch)
library(ggplot2)
library(ggtext)
library(gsignal)

visualSOZ <- function(epoch, sozNames) {
  p <- plot(epoch)

  elecColor <- rep("black", nrow(epoch))
  elecColor[rownames(epoch) %in% sozNames] <- 'red'
  elecColor <- rev(elecColor) # match the electrode order in the plot

  p +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1)+
    theme(axis.text.y = element_markdown(colour = elecColor))
}

dl <- EpochDownloader()
pt01sz1<-dl$Retrostudy_subpt01_1

pt01sz1Clipped<-crop(pt01sz1, start=-10, end=20)

pt01sozName <- rownames(pt01sz1Clipped)[rowData(pt01sz1Clipped)$soz]
pt01Display <- c(pt01sozName, "MLT1", "MLT2", "MLT3", "MLT4")
pt01sz1Reordered <- pt01sz1Clipped[pt01Display, ]
visualSOZ(pt01sz1Reordered, pt01sozName)


windowParams = c(1, 0.2)

#fs=1000
#epoch=pt01sz1m10p20s14e

EIPt01<-computeEpileptogenicIndex(epoch=pt01sz1Reordered, windowParams)


ermaster<-EIPt01$energyRatio


plotER<-plotERHeatmap(ei=EIPt01,sozIndex=sozIndex)
plotER
