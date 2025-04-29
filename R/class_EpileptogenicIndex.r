.EpileptogenicIndex<- setClass(
    "EpileptogenicIndex",
    slots = list(
        energyRatio= "matrix",
        epileptogenicIndex="numeric",
        timeDetect="numeric",
        startTimes = "numeric",
        electrodes = "character"
    )
)


EpileptogenicIndex<- function(energyRatio, epileptogenicIndex, timeDetect, startTimes, electrodes) {

  .EpileptogenicIndex(
    energyRatio = energyRatio,
    epileptogenicIndex=epileptogenicIndex,
    timeDetect=timeDetect,
    startTimes = startTimes,
    electrodes = electrodes
  )
}



#' Print the EpileptogenicIndex object
#' @param object A EpileptogenicIndex object
#' @rdname show-EpileptogenicIndex-method
#' @export
setMethod("show", "EpileptogenicIndex", function(object) {
  cat("\nEpileptogenicIndex object\n")
    slots <- c("energyRatio","epileptogenicIndex","timeDetect","startTimes","electrodes")
  printSlots(object, slots = slots)
  cat("Use '$attr' to access the data\n")
  invisible(object)
})


#' Get the number of rows or columns of a EpileptogenicIndex object
#'
#' @param x A EpileptogenicIndex object
#'
#' @rdname dim-EpileptogenicIndex-method
setMethod("nrow", "EpileptogenicIndex", function(x) {
  nrow(x@energyRatio)
})

#' @rdname dim-EpileptogenicIndex-method
setMethod("ncol", "EpileptogenicIndex", function(x) {
  ncol(x@energyRatio)
})


#' Subset a EpileptogenicIndex object
#'
#' @param x A EpileptogenicIndex object
#' @param i A logical vector or a numeric vector of indices to subset the electrodes
#' @param j A logical vector or a numeric vector of indices to subset the time windows
#' @param ... Additional arguments (not used)
#' @param drop Additional arguments (not used)
#' 
#' @rdname subset-EpileptogenicIndex-method
setMethod("[", "EpileptogenicIndex", function(x, i, j, ..., drop = FALSE) {
  
  if (!missing(i)){
    i <- checkIndex(i, x$electrodes)
  }else{
    i <- TRUE
  }
  if(missing(j)){
    j <- TRUE
  }
  
  energyRatio_subset <- x@energyRatio[i, j, drop = FALSE]
  epileptogenicIndex_subset<-x@epileptogenicIndex[i]
  timeDetect_subset<-x@timeDetect[i]
  startTimes_subset <- x@startTimes[j]
  electrodes_subset <- x@electrodes[i]
  .EpileptogenicIndex(
    energyRatio = energyRatio_subset,
    epileptogenicIndex=epileptogenicIndex_subset,
    timeDetect = timeDetect_subset,
    startTimes = startTimes_subset,
    electrodes = electrodes_subset,
  )
})

setMethod("$", "EpileptogenicIndex", function(x, name) {
  slot(x, name)
})

setMethod("$<-", "EpileptogenicIndex", function(x, name, value) {
  slot(x, name) <- value
  invisible(x)
})