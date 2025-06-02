standardizeIEEG <- function(data) {
  scaling <- 10^floor(log10(max(data)))
  plotData <- data / scaling
}

#' compute the Epileptogenic Index Object
#'
#' @param epoch Matrix or Epoch object. iEEG data matrix or Epoch object. If matrix, the row names are the electrode names and the column names are the time points
#' @param fs Numeric. frequency of signal iEEG acquisition
#' @param windowParams  Numeric. Window parameters of the multitaper method
#'
#' @return A power analysis matrix
#' @export
#'
#' @examples
#' data("pt01EcoG")
#' windowParams<-c(0.25,0.1)
#' epoch <- Epoch(pt01EcoG)
#' EIPt01<-computeEpileptogenicIndex <- function(epoch, windowParams)
computeEpileptogenicIndex <- function(epoch,  windowParams, fs=1000){

  # Set spectrogram parameters
  # define frequency bands
  thetaBand<-c(3.5,7.4)
  alphaBand<-c(7.4,12.4)
  betaBand<-c(12.4,24)
  gammaBand<-c(24,140)

  timeSeries<-tblData(epoch)

  # This number sort of serves as a built-in "time cost". High-frequency
  # activity needs to be robust enough to surpass this v number and still
  # overcome the threshold.
  # It's arbitrary but I think 0.5 is what Bartolomei used.
  v=0.5

  # Lambda also serves as a somewhat arbitrary threshold value.
  lambda=15

  timeNum <- ncol(timeSeries)
  elecNum <- nrow(timeSeries)

  nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1

  rangeBand<-c(3.5,140)
  data   <- vector(mode="numeric", length=timeNum)
  data[1:timeNum]<-timeSeries[1,1:timeNum]
  # Compute the multitaper spectrogram
  results = multitaperSpectrogramR(data=data, fs=fs, windowParams = windowParams, frequencyRange=rangeBand)
  stimes = results[[2]]

  times<-as.numeric(colnames(timeSeries))
  timesOnset<-results[[2]]+times[1]

  ermaster=matrix(0,elecNum,nwt)
  unmaster=matrix(0,elecNum,nwt)

  for(ie in 1:elecNum){

    data[1:timeNum]<-timeSeries[ie,1:timeNum]
    # Compute the multitaper spectrogram
    results = multitaperSpectrogramR(data, fs, windowParams = windowParams, frequencyRange=rangeBand)
    spec=results[[1]]
    sfreq=results[[3]]

    eralpha=computeMeanPowBand(results[[1]],alphaBand,sfreq)
    ertheta=computeMeanPowBand(results[[1]],thetaBand,sfreq)
    erbeta=computeMeanPowBand(results[[1]],betaBand,sfreq)
    ergamma=computeMeanPowBand(results[[1]],gammaBand,sfreq)

    er=(erbeta+ergamma)/(eralpha+ertheta)
    ern=cumsum(er)/c(1:nwt)
    un=er-ern-v
    un=cumsum(un)

    ermaster[ie,]=er
    unmaster[ie,]=un

  }

  stimes = results[[2]]

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

  maxNd<-nwt-hspan
  Nd[is.na(Nd)]=maxNd
  Nd[Nd>nwt-hspan]=maxNd

  EI<- vector(mode="numeric", length=elecNum)
  N0=min(Nd)
  t0=stimes[N0]+times[1]

  for(ie in 1:elecNum){
    EI[ie]=mean(ermaster[ie,Nd[ie]:Nd[ie]+hspan])/(stimes[Nd[ie]-t0+tau])
  }

  maxei=max(EI)
  EI=EI/maxei

  idNoDetect<-which(Nd==maxNd)
  td<-stimes[Nd]+times[1]
  td[idNoDetect]<-times[timeNum]


  EpileptogenicIndex(
    energyRatio= ermaster,
    epileptogenicIndex=EI,
    timeDetect<-td,
    startTimes = timesOnset,
    electrodes = rownames(timeSeries)
  )

}

computeMeanPowBand <- function(spect, band, sfreq){

  freqStart=band[1]
  freqEnd=band[2]
  freqStartIndex <- which.min(abs(sfreq - freqStart))
  freqEndIndex  <- which.min(abs(sfreq - freqEnd))
  spectBand = spect[freqStartIndex:freqEndIndex,]
  erBand=colSums(spectBand)

  erBand

}
