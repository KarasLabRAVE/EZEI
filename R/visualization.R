
#' Visualization of ictal iEEG
#'
#' @return A ggplot object
#'
#' @examples
#' data("pt01EcoG")
#' 
#' ## Visualize a subject of electrodes
#' sozIndex <- attr(pt01EcoG, "sozIndex")
#' display <- c(sozIndex, 77:80)
#' 
#' epoch <- Epoch(pt01EcoG)
#' visuIEEGData(epoch = epoch[display, ])
#' @export
visuIEEGData <- function(epoch) {
  if (is(epoch, "matrix")){
    epoch <- Epoch(epoch)
  }
  
  gaps <- 2
  
  elecNames <- epoch$electrodes
  data <- epoch$data
  elecNum <- nrow(data)
  timesNum <- ncol(data)
  
  plotData <- standardizeIEEG(data)
  
  times <- epoch$times
  if (is.null(times)) {
    xlabel <- "Time Index"
    timeTicks <- seq_len(timesNum)
  } else {
    xlabel <- "Time (s)"
    timeTicks <- times
  }
  
  plotData <- apply(plotData, 1, function(x) x - mean(x))
  plotData <- as.data.frame(plotData)
  plotData$timeTicks <- timeTicks
  breakplot <- (seq_len(elecNum) - 1) * gaps
  
  elecNamesReversed <- rev(elecNames)
  
  ## add gaps between electrodes
  for (i in seq_along(elecNamesReversed)) {
    elec <- elecNamesReversed[i]
    plotData[[elec]] <- plotData[[elec]] + (i-1) * gaps
  }
  
  
  p <- ggplot2::ggplot(data = plotData)
  for (i in seq_along(elecNamesReversed)) {
    elec <- elecNamesReversed[i]
    p <- p + ggplot2::geom_line(ggplot2::aes(x = .data$timeTicks, y = .data[[elec]]))
  }
  
  p +
    ggplot2::labs(x = xlabel, y = "Electrode", size = 2) +
    ggplot2::scale_y_continuous(labels = elecNamesReversed, breaks = breakplot)
}
