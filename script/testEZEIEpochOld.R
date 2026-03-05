# test EZEI  with Epoch package
# data download one patient one seizure




dl <- EpochDownloader()
pt01sz1<-dl$Retrostudy_subpt01_1

#pt01sz1ieeg<-tblData(pt01sz1)

# crop
pt01sz1m10p20s<-crop(pt01sz1, start=-10, end=20)
sozIndex<-which(pt01sz1m10p20s@rowData$soz==TRUE)

display <- c(sozIndex, 77:80)
# subset on 14 electrodes
pt01sz1m10p20s14e<-pt01sz1m10p20s[display,]
sozIndex<-which(pt01sz1m10p20s14e@rowData$soz==TRUE)

windowParams = c(1, 0.2)

#fs=1000
#epoch=pt01sz1m10p20s14e

EIPt01<-computeEpileptogenicIndex(epoch=pt01sz1m10p20s14e, windowParams)


ermaster<-EIPt01$energyRatio


plotER<-plotERHeatmap(ei=EIPt01,sozIndex=sozIndex)
plotER
