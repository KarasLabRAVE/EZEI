# A plot function that takes a data frame and returns a heatmap plot
makeHeatMap <- function(df, xTicksNum = 10){
  xLabels <- colnames(df)
  yLabels <- rownames(df)
  
  if(is.null(xLabels)){
    xLabels <- seq_len(ncol(df))
  }
  if(is.null(yLabels)){
    yLabels <- seq_len(nrow(df))
  }
  
  df$y <- yLabels
  df_long <- reshape2::melt(df, id.vars = "y", variable.name = "x", value.name = "value")
  colnames(df_long) <- c("y", "x", "value")
  
  ## sort df_long by rownames(ERMatReorderd)
  df_long$x <- factor(df_long$x, levels = xLabels)
  df_long$y <- factor(df_long$y, levels = rev(yLabels))
  ## show 10 time points on x-axis at most
  if (length(xLabels) > xTicksNum){
    step <- ceiling(length(xLabels) / xTicksNum)
    breaksIdx <- seq(1, length(xLabels), by = step)
    breaks <- xLabels[breaksIdx]
  } else {
    breaks <- xLabels
  }
  
  ggplot2::ggplot(df_long) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)) +
    ggplot2::scale_x_discrete(labels = breaks, breaks = breaks) +
    ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5)) +
    viridis::scale_fill_viridis(option = "turbo") +
    ggplot2::theme_minimal()
}

#' Visualization functions ( ER matrix)
#'
#' @description `plotERHeatmap`: plot epileptogenic ratio heatmaps with electrodes marked as soz colored
#'
#' @param ermaster ermaster object from 
#' @param sozIndex Integer or string. A group of electrodes to mark as in the Seizure Onset Zone (SOZ)
#' 
#' @return A ggplot object
#'
#' @examples
#' 
#' data("pt01EcoG")
#' 
#' ## sozIndex is the index of the electrodes we assume are in the SOZ
#' sozIndex <- attr(pt01EcoG, "sozIndex")
#' 
#' ## precomputed fragility object
#' data("pt01Frag")
#' 
#' ## plot the fragility heatmap
#' plotFragHeatmap(frag = pt01Frag, sozIndex = sozIndex)
#' 
#' @rdname plotERHeatmap
#' @export
plotERHeatmap <- function(
    ER,
    sozIndex = NULL) {
  ## TODO: make sozID an optional
  ## TODO: add plot support to ER
  ERMat <- ER
  elecNum <- nrow(ERMat)
  windowNum <- ncol(ERMat)
  
  elecNames <- rownames(ERMat)
  sozIndex <- checkIndex(sozIndex, elecNames)
  
  group1 <- sozIndex
  group2 <- setdiff(seq_len(elecNum), sozIndex)
  
  elecColor <- rep("blue", elecNum)
  elecColor[seq_along(group2)] <- "black"
  
  startTime <- colnames(ERMat)
  if (is.null(startTime)) {
    xlabel <- "Time Index"
    stimes <- seq_len(windowNum)
  } else {
    xlabel <- "Time (s)"
    stimes <- startTime
  }
  
  colnames(ERMat) <- stimes
  
  ## prepare the data.frame for visualization
  allIndex <- c(group1, group2)
  df <- as.data.frame(ERMat[allIndex, ])
  
  
  
  makeHeatMap(df) +
    ggplot2::labs(x = xlabel, y = "Electrode", size = 2) +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(size = 6, colour = elecColor), # Adjust depending on electrodes
    )
}



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
