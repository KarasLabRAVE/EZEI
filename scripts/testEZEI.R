data("pt01EcoG")
timeWindow <- c(-10, 20)
epoch <- Epoch(pt01EcoG)
fs=1000
sozIndex <- attr(pt01EcoG, "sozIndex")
timeNum <- ncol(epoch)
windowParams = c(1, 0.2) 

nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1
data   <- vector(mode="numeric", length=timeNum)
data[1:timeNum]<-epoch$data[sozIndex[1],1:timeNum]
# Compute the multitaper spectrogram
results = multitaperSpectrogramR(data, fs)
## Visualize a subject of electrodes
#display <- c(sozIndex, 77:80)
stimes = results[[2]]
epoch <- Epoch(pt01EcoG)
visuIEEGData(epoch)


# Set spectrogram params

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
elecNames <- epoch$electrodes
timePoints <- epoch$times

nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1



#ComputeERMaster

stimesgp=stimes+timeWindow[1]
ermaster=matrix(0,elecNum,nwt)
unmaster=matrix(0,elecNum,nwt)



for(ie in 1:elecNum){
  data[1:timeNum]<-epoch$data[ie,1:timeNum]
  # Compute the multitaper spectrogram
  results = multitaperSpectrogramR(data, fs)
  spec=results[[1]]
  sfreq=results[[3]]

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

rownames(ermaster)<-rownames(pt01EcoG)
colnames(ermaster)<-stimesgp

plotER<-plotERHeatmap(ermaster,sozIndex=sozIndex)
plotER


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

EI<- vector(mode="numeric", length=elecNum)
N0=min(Nd)
t0=stimes[N0]+timeWindow[1]

for(ie in 1:elecNum){
  EI[ie]=mean(ermaster[ie,Nd[ie]:Nd[ie]+hspan])/(stimes[Nd[ie]-t0+tau])
}

maxei=max(EI)
EI=EI/maxei

#channelname<-channelDefs$name
votethres   <- vector(mode="numeric", length=elecNum)

for(ie in 1:elecNum){
  #print(EI[ie])
  if(EI[ie]>0.3){
    votethres[ie]=1
  }
}



