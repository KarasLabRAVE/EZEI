## Patient K006 with TBI
library(gsignal)
library(R.matlab)
library(readxl)
data <- readMat('/Users/aclesage/Documents/RAVEProjects/TBIDatam30sp30s/K006sz1_m30sp30s.mat')
#K006Epoch<- data$a
dumm<-t(data$a)  
fs<-512
dumm2<-gsignal::resample(dumm,1000,fs)
elec_resamp<-t(dumm2)
K006Epoch<-elec_resamp[,20001:50001]

fs<-1000

## add channel names to the rows
goodChannels <- c(1:26,29:59,61:177,179:242)
sozChannels<-c(14,23:26,37:38,47:49)
channelNames <- read_excel('/Users/aclesage/Documents/RAVEProjects/TBIDatam30sp30s/K006IctalSEEGChannels.xls')
rownames(K006Epoch) <- channelNames$name[goodChannels]


sozIndex<-which(goodChannels%in%sozChannels==TRUE)
sozNames<-channelNames$name[sozChannels]

attr(K006Epoch, "sozIndex") <- sozIndex
attr(K006Epoch, "sozNames") <- sozNames


## Add time stamps to the columns
times <- seq(-10, 20, length.out=ncol(K006Epoch))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
colnames(K006Epoch)<-times_with_sign

epoch <- Epoch(K006Epoch)
windowParams = c(1, 0.2) 
display<-c(sozIndex)
visuIEEGData(epoch[display,])

K006EI<-computeEpileptogenicIndex(epoch[1:60,], windowParams)

plotER<-plotERHeatmap(ei=K006EI,sozIndex=sozIndex)
plotER
