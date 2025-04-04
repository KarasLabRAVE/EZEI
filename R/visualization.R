
#' Visualization of ictal iEEG
#'
#' @param ieegts Matrix or Fragility object. Either a matrix of iEEG time
#' series x(t), with time points as rows and electrodes names as columns,
#' or a Fragility object from \code{calcAdjFrag}
#' @param title String. Figure title
#' @param display Integer or string. Vector electrodes to display
#' @return plot raw signal
#'
#' @examples
#'data("pt01Epoch")
#'sozIndex <- attr(pt01Epoch,"sozIndex")
#'display <- c(sozIndex,77:80)
#'timeRange <- c(-10,10)
#'iEEGplot<-visuIEEGData(ieegts=pt01Epoch,timeRange=timeRange,display=display)
#'iEEGplot
#' @export
visuIEEGData<-function(ieegts, timeRange=NULL, title = "Patient name seizure number", display=NULL){

  titlepng<- title
  if(is.null(display)){
    display<-1:ncol(ieegts)
  }


  elecName<-colnames(ieegts)
  if(typeof(display)=="integer"){

    displayTot<-1:ncol(ieegts)
    diffDisplayTot<-setdiff(display,displayTot)

    if(length(diffDisplayTot)!=0){
      listDisplayMissing<-paste(as.character(diffDisplayTot),collapse=" ")
      message<-paste("Display electrodes indices. Numbers ",listDisplayMissing,"are out of electrode number limit")
      warning(message)
      display<-display[!display%in%diffDisplayTot]
      displayCor<-paste(as.character(display),collapse=" ")
      message<-paste("Keeping indices.",displayCor)
      warning(message)

    }

    displayid<-display

  }else{

    diffDisplayTot<-setdiff(display,elecName)

    if(length(diffDisplayTot)!=0){
      listDisplayMissing<-paste(diffDisplayTot,collapse=" ")
      message<-paste("Display electrodes names. Names ",listDisplayMissing,"are out of name list")
      warning(message)
      display<-display[!display%in%diffDisplayTot]
      displayCor<-paste(display,collapse=" ")
      message<-paste("Keeping names.",displayCor)
      warning(message)

    }

    displayid<-which(elecName%in%display)
  }

  scaling <- 10^floor(log10(max(ieegts)))
  plotData<-ieegts[,displayid]/scaling
  gaps<-2
  displayNames<-colnames(ieegts)[displayid]
  nElec<-length(displayid)
  nt<-nrow(plotData)
  if(is.null(timeRange)){
    xlabel<-"Time Index"
    stimes<-seq_len(nt)
  }else{
    xlabel<-"Time (s)"
    stimes<-seq(timeRange[1],timeRange[2],length.out=nt)
  }


  for(i in 1:ncol(plotData)){
    plotData[, i] <- (plotData[, i]- mean(plotData[, i]))+
      (ncol(plotData)-i)*gaps
  }
  plotData<-data.frame(plotData)
  breakplot<-(c(1:nElec)-1)*gaps

  p<-ggplot2::ggplot(data=plotData,ggplot2::aes(x=stimes,y=plotData))+
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = xlabel, y = "Electrode",size=2)+
    ggplot2::geom_vline(xintercept =0,
                        color = "black", linetype = "dashed", size = 1)

  for(i in 1:nElec){
    p<-p+ggplot2::geom_line(ggplot2::aes_string(y=names(plotData)[i]))
  }
  displayNames<-rev(displayNames)
  p<-p+ggplot2::scale_y_continuous(labels=displayNames,breaks=breakplot)

  return(p)

}
