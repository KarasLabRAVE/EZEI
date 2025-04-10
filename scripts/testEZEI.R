data("pt01EcoG")

## Visualize a subject of electrodes
sozIndex <- attr(pt01EcoG, "sozIndex")
display <- c(sozIndex, 77:80)

epoch <- Epoch(pt01EcoG)
visuIEEGData(epoch = epoch[display, ])

timeWindow <- c(-2, 5)

timeBandwidth = 3  # Set time-half bandwidth
numTapers = 5  # Set number of tapers (optimal is time_bandwidth*2 - 1)
windowParams = c(1, 0.2)  # Window size is 1s with step size of 0.2s
minNfft = 0  # No minimum nfft
weighting = 'unity'  # weight each taper at 1
detrendOpt = 'off'  # detrend each window by subtracting the average
parallel = FALSE  # use multiprocessing
numWorkers = 3  # use 3 cores in multiprocessing
plotOn = FALSE  # plot spectrogram
verbose = FALSE  # print extra info
xyflip = FALSE  # do not transpose spect output matrix


# Set spectrogram params
frequencyRange = c(0.5, 250)  # Frequency Band
thetaBand<-c(3.5,7.4)
alphaBand<-c(7.4,12.4)
betaBand<-c(12.4,24) #[13-30]Hz
gammaBand<-c(24,140) #[30-90]Hz


# This number sort of serves as a built-in "time cost". High-frequency
# activity needs to be robust enough to surpass this v number and still
# overcome the threshold.
# It's arbitrary but I think 0.5 is what Bartolomei used.
v=0.5

# Lambda also serves as a somewhat arbitrary threshold value.
lambda=15

elecNum <- nrow(epoch)
timeNum <- ncol(epoch)
elecNames <- epoch$electrodes
timePoints <- epoch$times
dataMat <- epoch$data

fs=1000

nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1

#nwt=295
data   <- vector(mode="numeric", length=timeNum)
data[1:timeNum]<-dataMat[1,1:timeNum]
# Compute the multitaper spectrogram
results = multitaperSpectrogramR(data, fs, frequencyRange, timeBandwidth, numTapers, windowParams, minNfft, weighting, detrendOpt, parallel, numWorkers,
                                   plotOn, verbose, xyflip)

stimes = results[[2]]
nwt=length(stimes)
ermaster=matrix(0,elecNum,nwt)
unmaster=matrix(0,elecNum,nwt)


#ComputeERMaster

maxspec=0

stimesgp=stimes+timeWindow[1]

for(ie in 1:elecNum){
  
  data[1:timeNum]<-dataMat[ie,1:timeNum]
  # Compute the multitaper spectrogram
  results = multitaperSpectrogramR(data, fs, frequencyRange, timeBandwidth, numTapers, windowParams, minNfft, weighting, detrendOpt, parallel, numWorkers,
                                   plotOn, verbose, xyflip)
  spectrogram=nanpow2db(results[[1]])
  if(max(spectrogram)>maxspec){
    maxspec=max(spectrogram)
  }
}

for(ie in 1:elecNum){
  #ie=40
  #print(ie)
  data[1:timeNum]<-dataMat[ie,1:timeNum]
  # Compute the multitaper spectrogram
  results = multitaperSpectrogramR(data, fs, frequencyRange, timeBandwidth, numTapers, windowParams, minNfft, weighting, detrendOpt, parallel, numWorkers,
                                   plotOn, verbose, xyflip)
  spec=results[[1]]
  sfreq=results[[3]]
  # spec=spec[-c(57,67),]
  # sfreq=sfreq[-c(57:67)]
  
  spectrogram=nanpow2db(results[[1]])
  
  tfdf<-data.frame(spectrogram)
  colnames(tfdf)<-stimesgp
  rownames(tfdf)<-sfreq
  tfmap_data <- expand.grid(Time = stimesgp, Frequency=sfreq)
  tfmap_data$Value <- c(t(spectrogram))
  # 
  # 
  # titlepng=paste(subject_code,'Seizure',as.character(j),channelDefs$name[ie],sep=" ")
  # 
  # ggplot(tfmap_data, aes(x = Time, y = Frequency, fill = Value)) +
  #   geom_tile() +
  #   ggtitle(titlepng)+
  #   labs(x = "Time (s)", y = "Frequency") +
  #   scale_fill_viridis(option = "turbo",limits=c(20,maxspec)) +  #
  #   theme_minimal()
  # 
  # resfile=paste(pathres,'TF_',subject_code,'_Seizure',as.character(j),'_e',as.character(electrodes[ie]),'_',channelDefs$name[ie],'.png',sep="")
  # ggsave(resfile)
  
  
  freqStart=thetaBand[1]
  freqEnd=thetaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  specttheta = results[[1]][freqStartIndex:freqEndIndex,]
  ertheta=colSums(specttheta)
  
  freqStart=alphaBand[1]
  freqEnd=alphaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectalpha = results[[1]][freqStartIndex:freqEndIndex,]
  eralpha=colSums(spectalpha)
  
  freqStart=betaBand[1]
  freqEnd=betaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectbeta = results[[1]][freqStartIndex:freqEndIndex,]
  erbeta=colSums(spectbeta)
  
  freqStart=gammaBand[1]
  freqEnd=gammaBand[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectgamma = results[[1]][freqStartIndex:freqEndIndex,]
  ergamma=colSums(spectgamma)
  
  er=(erbeta+ergamma)/(eralpha+ertheta)
  ern=cumsum(er)/c(1:nwt)
  un=er-ern-v
  un=cumsum(un)
  
  ermaster[ie,]=er
  unmaster[ie,]=un
  
}


stimes = results[[2]]
sfreq=results[[3]]

Nd   <- vector(mode="numeric", length=elecNum)
Na   <- vector(mode="numeric", length=elecNum)

Nd[1:elecNum]=10*nwt
Na[1:elecNum]=nwt

# plot(stimes,ermaster[31,],type='lines')
# plot(stimes,unmaster[31,],type='lines')

for(it in 2:nwt){
  
  unt=unmaster[,1:it]
  un=apply(unt,1,FUN=min)
  ind=apply(unt,1,which.min)
  
  undiff=unmaster[,it]-un
  pastThreshold=undiff>lambda
  pastThreshold[Nd!=10*nwt]=FALSE
  
  ie=which(pastThreshold==TRUE)
  Nd[ie]=ind[ie]
  Na[ie]=it
  
}


tau=1
H=5

hspan=which.min(abs(stimes-H))

Nd[is.na(Nd)]=nwt-hspan
Nd[Nd>nwt-hspan]=nwt-hspan
#
EI<- vector(mode="numeric", length=elecNum)
N0=min(Nd)
t0=stimes[N0]

for(ie in 1:elecNum){
  EI[ie]=mean(ermaster[ie,Nd[ie]:Nd[ie]+hspan])/(stimes[Nd[ie]]-t0+tau)
}

t0=t0+timeWindow[1]

maxei=max(EI)
EI=EI/maxei


