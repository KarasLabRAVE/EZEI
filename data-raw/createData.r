
## load pt01epochdata.mat
## Patient PT01 from the Fragility data set

library(R.matlab)
library(readxl)
data <- readMat('data-raw/subpt01_seizure1m10sp20s.mat')
pt01EpochRaw <- data$data

## add channel names to the rows
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
sozChannels<-c(33:34,62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01EpochRaw) <- channelNames$name[goodChannels]

sozIndex<-which(goodChannels%in%sozChannels==TRUE)
sozNames<-channelNames$name[sozChannels]

## Add time stamps to the columns
times <- seq(-10, 20, length.out=ncol(pt01EpochRaw))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
colnames(pt01EpochRaw)<-times_with_sign


#pt01EcoG<-pt01EpochRaw[,8001:15000]
#pt01EcoG<-pt01EpochRaw
display <- c(sozIndex, 77:80)
pt01EcoG<-pt01EpochRaw[display,]
sozIndex<-c(1:10)
attr(pt01EcoG, "sozIndex") <- sozIndex
attr(pt01EcoG, "sozNames") <- sozNames
usethis::use_data(pt01EcoG, overwrite = TRUE)




