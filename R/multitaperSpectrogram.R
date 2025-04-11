#' Multitaper Spectrogram functions 
#'
#' Compute multitaper spectrogram of timeseries data
#' 
#' This R script contains the R implementations of the multitaper spectrogram analysis described in the paper "Sleep Neurophysiological Dynamics Through the 
#' Lens of Multitaper Spectral Analysis". 
#' Multitaper spectral estimation, which was developed in the early 1980
#' s by David Thomson2 and has been shown to have superior statistical properties compared with single-taper spectral estimates
#' 
#' https://github.com/preraulab/multitaper_toolbox
#'
#' Results tend to agree with Prerau Lab python implementation of multitaper spectrogram with precision on the order of at most
#' 10^-7 with SD of at most 10^-5
#'
#' @param data (numeric vector): time series data -- required
#' @param fs (numeric): sampling frequency in Hz  -- required
#' @param frequencyRange (numeric vector): c(<min frequency>, <max frequency>) (default: NULL, adjusted to
#'                                           c(0, nyquist) later)
#' @param timeBandwidth (numeric): time-half bandwidth product (window duration*half bandwidth of main lobe)
#'                                   (default: 5 Hz*s)
#' @param numTapers (numeric): number of DPSS tapers to use (default: NULL [will be computed
#'                                                               as floor(2*time_bandwidth - 1)])
#' @param windowParams  (numeric vector): c(window size (seconds), step size (seconds)) (default: [5 1])
#' @param minNfft (numeric) minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
#' @param weighting (char) weighting of tapers ('unity' (default), 'eigen', 'adapt')
#' @param detrendOpt (char): detrend data window ('linear' (default), 'constant', 'off')
#' @param parallel (logical): use parallel processing to speed up calculation (default: FALSE). Note: speedup is faster on
#                             unix-like machines (Mac, Linux) because they allow fork processes while Windows does not.
#         num_workers (numeric): number of cpus/workers to dedicate to parallel processing (default: FALSE). Note: Will
#                                be ignored if parallel is FALSE. If parallel is TRUE and num_workers is false (or if num_workers
#                                exceeds available workers), will default to max number of workers available minus 1.
#' @param numWorkers (logical): plot results (default: TRUE)
#' @param plotOn (logical): plot results (default: TRUE)
#' @param verbose (logical): display spectrogram properties (default: TRUE)
#' @param xyflip ((logical): return the transpose of mt_spectrogram)
#'
#' @return
#'        mt_spectrogram (matrix): spectral power matrix
#'         stimes (numeric vector): timepoints (s) in mt_spectrogram
#'         sfreqs (numeric vector): frequency values (Hz) in mt_spectrogram
#' @export
#'
#' @examples
#'data("pt01EcoG")
#'
#'timeWindow <- c(-10, 20)
#'epoch <- Epoch(pt01EcoG)
#' fs=1000
#' timeNum <- ncol(epoch)
#'windowParams = c(1, 0.2) 
#' nwt=floor((timeNum/fs-windowParams[1])/windowParams[2])+1
#' data   <- vector(mode="numeric", length=timeNum)
#' data[1:timeNum]<-dataMat[sozIndex[1],1:timeNum]
#' # Compute the multitaper spectrogram
#' results = multitaperSpectrogramR(data, fs, frequencyRange, timeBandwidth, numTapers, windowParams, minNfft, weighting, detrendOpt, parallel, numWorkers,
#'                                  plotOn, verbose, xyflip)
multitaperSpectrogramR <- function(data, fs, frequencyRange=c(0.5,250), timeBandwidth=3, numTapers=5, windowParams=c(1,0.2),
                                     minNfft=0, weighting='unity', detrendOpt='off', parallel=FALSE, numWorkers=3,
                                     plotOn=FALSE, verbose=FALSE, xyflip=FALSE){
 # returns:
  
  # Process user input
  res <- process_input(data, fs, frequencyRange, timeBandwidth, numTapers, windowParams, minNfft, weighting, detrendOpt,
                       plotOn, verbose)

  data <- res[[1]]
  fs <- res[[2]]
  frequencyRange <- res[[3]]
  timeBandwidth <- res[[4]]
  numTapers <- res[[5]]
  winsize_samples <- res[[6]]
  winstep_samples = res[[7]]
  window_start = res[[8]]
  num_windows <- res[[9]]
  nfft <- res[[10]]
  weighting = res[[11]]
  detrendOpt <- res[[12]]
  plotOn <- res[[13]]
  verbose <- res[[14]]

  # Set up spectrogram parameters
  res <- process_spectrogram_params(fs, nfft, frequencyRange, window_start, winsize_samples)
  window_idxs <- res[[1]]
  stimes <- res[[2]]
  sfreqs <- res[[3]]
  freq_inds <- res[[4]]

  # Display spectrogram parameters if desired
  if(verbose){
    display_spectrogram_properties(fs, timeBandwidth, numTapers, c(winsize_samples, winstep_samples), frequencyRange,
                                   detrendOpt)
  }

  # Split data into window segments
  data_segments <- t(sapply(window_idxs, split_data_helper, data=data))

  # COMPUTE THE MULTITAPER SPECTROGRAM
  #     STEP 1: Compute DPSS tapers based on desired spectral properties
  #     STEP 2: Multiply the data segment by the DPSS Tapers
  #     STEP 3: Compute the spectrum for each tapered segment
  #     STEP 4: Take the mean of the tapered spectra

  tic <- proc.time() # start timer for multitaper

  # Compute DPSS tapers (STEP 1)
  dpss_tapers <- multitaper::dpss(winsize_samples, numTapers, timeBandwidth, returnEigenvalues=TRUE)
  dpss_eigen = dpss_tapers$eigen
  dpss_tapers = dpss_tapers$v

  # pre-compute weights
  if(weighting == 'eigen'){
    wt = dpss_eigen / numTapers;
  }else if(weighting == 'unity'){
    wt = calibrate::ones(numTapers,1) / numTapers;
  }else{
    wt = 0;
  }

  # Compute multitaper #
  if(parallel){ # Check for parallelization
    workers_avail <- parallel::detectCores() - 1  # detect cores available and leave 1 for user
    if(numWorkers==FALSE | numWorkers > workers_avail){
      numWorkers = workers_avail
      warning(paste("Number of workers for parallelization either not specified or greater than workers available. Setting number
              of workers to number available minus 1 (", toString(numWorkers), ")"))
    }
    doParallel::registerDoParallel(cores=numWorkers) # register workers with doParallel

    # Create cluster of workers differently depending on OS
    if(.Platform$OS.type == "windows"){ # windows cannot use FORK argument
      cluster <- parallel::makeCluster(numWorkers) # create cluster of workers without forking
    }else{
      cluster <- parallel::makeCluster(numWorkers, type="FORK") # if not windows, use FORK because it's faster
    }

    mt_spectrogram <- parallel::parApply(cluster, data_segments, 1, calc_mts_segment, dpss_tapers=dpss_tapers, nfft=nfft, freq_inds=freq_inds,
                               weighting=weighting, wt=wt, dpss_eigen=dpss_eigen, numTapers=numTapers, detrendOpt=detrendOpt)

    ParallelLogger::stopCluster(cluster) # stop cluster to give back resources
    registerDoSEQ() # switch back to serial processing
  }else{ # if no parallelization, use normal apply
    mt_spectrogram = apply(data_segments, 1, calc_mts_segment, dpss_tapers=dpss_tapers, nfft=nfft, freq_inds=freq_inds,
                           weighting=weighting, wt=wt, dpss_eigen=dpss_eigen, numTapers=numTapers, detrendOpt=detrendOpt)
  }


  # Compute one-sided PSD spectrum
  DC_select = which(sfreqs==0)
  Nyquist_select = which(sfreqs==fs/2)
  select = setdiff(1:(length(sfreqs)), c(DC_select, Nyquist_select))
  mt_spectrogram = rbind(mt_spectrogram[DC_select,], 2*mt_spectrogram[select,], mt_spectrogram[Nyquist_select,]) / fs

  # End timer and get elapsed time
  toc = proc.time()
  elapsed = toc-tic
  if(verbose){
    print(paste("Multitaper compute time: ", toString(round(elapsed[[3]], digits=5)), " seconds", sep=""))
  }


  if(all(as.vector(mt_spectrogram) == 0)){
    print("Spectrogram calculated as all zeros, no plot shown")
  }else if(plotOn){
    print("plotting...")

    # Saving to PNG, loading back in, and plotting as raster is actually faster than just plotting using image.plot
    png(filename=paste("spectrogram.png")) # save as png
    fields::image.plot(x=stimes, y=sfreqs, nanpow2db(t(mt_spectrogram)), xlab="Time (s)",
               ylab='Frequency (Hz)')
    #image.plot(x=stimes, y=sfreqs, t(mt_spectrogram), xlab="Time (s)",
    #           ylab='Frequency (Hz)')
    dev.off()

    im <- readPNG("spectrogram.png") # load png
    file.remove("spectrogram.png") # remove png file
    plot.new()
    rasterImage(im,0,0,1,1, interpolate=FALSE) # plot as raster image
    print("done plotting")
  }




  if(xyflip){
    mt_spectrogram = t(mt_spectrogram)
  }

  return(list(mt_spectrogram, stimes, sfreqs))
}

split_data_helper <- function(indices, data){ # for sapply when splitting data into windows
  data_seg = data[indices]
  return(data_seg)
}



### Helper Functions ###

# Process user input #
process_input <- function(data, fs, frequencyRange=NULL, timeBandwidth=5, numTapers=NULL,
                          windowParams=c(5,1), minNfft=0, weighting='unity', detrendOpt='linear', plotOn=TRUE,
                          verbose=TRUE){

  # Helper function to process multitaper_spectrogram arguments, mainly checking for validity
  #
  # Params:
  #        data (numeric vector): time series data -- required
  #        fs (numeric): sampling frequency in Hz -- required
  #        frequency range (numeric vector): c(<min frequency>, <max frequency>) (default: c(0 nyquist))
  #        timeBandwidth (numeric): time-half bandwidth product (window duration*half bandwidth of main lobe) (default: 5 Hz*s)
  #        numTapers (numeric): number of DPSS tapers to use (default None [will be computed as floor(2*timeBandwidth - 1)])
  #        windowParams (numeric vector): c(window size (seconds), step size (seconds)) default: c(5,1)
  #        detrendOpt (char): detrend data window ('linear' (default), 'constant', 'off')
  #        minNfft: (numeric): minimum allowable NFFT size, adds zero padding for interpolation (default: 0)
  #        plotOn: (logical): plot results (default: TRUE)
  #        verbose (logical)L display spectrogram properties (default; TRUE)
  #
  # Returns:
  #         data (numeric vector) same as input
  #         fs (numeric): same as input
  #         frequencyRange (numeric vector): same as input or calculated from fs if not given
  #         timeBandwidth (numeric): same as input or default if not given
  #         numTapers (numeric): same as input or calculated from time bandwidth if not given
  #         winsize_samples (numeric): number of samples in a single time window
  #         winstep_samples (numeric): number of samples in a single window step
  #         window_start (numeric vector): matrix of timestamps representing the beginning time for each window
  #         num_windows (numeric): number of total windows
  #         nfft (numeric): length of signal to calculate fft on
  #         detrendOpt (char): same as input or default if not given
  #         plotOn (logical): same as input or default if not given
  #         verbose (logical): same as input or default if not given


  # Make sure data is 1D atomic vector
  if((is.atomic(data) == FALSE) | is.list(data)){
    stop("data must be a 1D atomic vector")
  }

  # Set frequency range if not provided
  if(is.null(frequencyRange)){
    frequencyRange <- c(0, fs/2)
  }

  # Set detrend method
  detrendOpt = tolower(detrendOpt)
  if(detrendOpt != 'linear'){
    if(detrendOpt == 'const' || detrendOpt == 'constant'){
      detrendOpt <- 'constant'
    } else if(detrendOpt == 'none' || detrendOpt == 'false' || detrendOpt == 'off'){
      detrendOpt <- 'off'
    }else{
      stop(paste("'", toString(detrendOpt), "' is not a valid detrendOpt argument. The",
                 " choices are: 'constant', 'linear', or 'off'.", sep=""))
    }
  }

  # Set taper weighting options
  weighting = tolower(weighting)
  if(weighting == 'adaptive' || weighting == 'adapt'){
    weighting = 'adapt'
  } else if(weighting == 'eig' || weighting == 'eigen'){
    weighting = 'eigen'
  } else if(weighting != 'unity'){
    stop(paste("'", toString(weighting), "' is not a valid weighing argument. Choices are: 'unity', 'eigen' or 'adapt'"))
  }


  # Check if frequency range is valid
  if(frequencyRange[2] > fs/2){
    frequencyRange[2] <- fs/2
    warning(paste("Upper frequency range greater than Nyquist, setting range to [",
                  toString(frequencyRange[1]), ",", toString(frequencyRange[2]), "].",
                  sep=""))
  }

  # Set number of tapers if none provided
  optimal_numTapers = floor(2*timeBandwidth) - 1
  if(is.null(numTapers)){
    numTapers <- optimal_numTapers
  }

  # Warn if number of tapers is suboptimal
  if(numTapers != optimal_numTapers){
    warning(paste("Suboptimal number of tapers being used. Number of tapers is optimal at floor(2*TW) - 1 which is ",
                  toString(optimal_numTapers), " in this case.", sep=""))
  }


  # Check if window size is valid, fix if not
  if((windowParams[1]*fs) %% 1 != 0){
    winsize_samples <- round(windowParams[1]*fs)
    warning(paste("Window size is not divisible by sampling frequency. Adjusting window",
                  " size to ", toString(winsize_samples/fs), " seconds.", sep=""))
  } else{
    winsize_samples <- windowParams[1]*fs
  }

  # Check if window step size is valid, fix if not
  if((windowParams[2]*fs) %% 1 != 0){
    winstep_samples <- round(windowParams[2]*fs)
    warning(paste("Window step size is not divisible by sampling frequency. Adjusting window",
                  " step size to ", toString(winstep_samples/fs), " seconds.", sep=""))
  } else{
    winstep_samples <- windowParams[2]*fs
  }

  # Get total data length
  len_data = length(data)

  # Check if length of data is smaller than window (bad)
  if(len_data < winsize_samples){
    stop(paste("Data length (", toString(len_data), ") is shorter than the window size (",
               toString(winsize_samples), "). Either increase data length or decrease",
               " window size.", sep=""))
  }

  # Find window start indices and num of windows
  window_start = seq(1, len_data-winsize_samples+1, by=winstep_samples)
  num_windows = length(window_start)

  # Get num points in FFT
  nfft = max(max(2^ceiling(log2(abs(winsize_samples))), winsize_samples), 2^ceiling(log2(abs(minNfft))))

  return(list(data, fs, frequencyRange, timeBandwidth, numTapers, winsize_samples, winstep_samples,
              window_start, num_windows, nfft, weighting, detrendOpt, plotOn, verbose))
}



#' Process spectrogram inputs.
#' Helper function to create frequency vector and window indices
#'
#' @param fs (numeric): sampling frequency in Hz  -- required
#' @param nfft (numeric): length of signal to calculate fft on -- required
#' @param frequencyRange (numeric vector): c(<min frequency>, <max frequency>) (default: NULL, adjusted to
#'                                           c(0, nyquist) later)
#' @param window_start (numeric vector): timestamps representing the beginning time for each window -- required
#' @param datawin_size datawin_size (numeric): seconds in one window -- required
#'
#' @return
#'        window_idxs (matrix): indices of timestamps for each window (nxm where n=number of windows and m=datawin_size)
#'        stimes (numeric vector): times for the centers of the spectral bins (1xt)
#'        sfreqs (numeric vector): frequency bins for spectrogram (1xf)
#'        freq_inds (logical vector): indicates which frequencies are being analyzed in an array of frequencies from 0 to fs
#'                    with steps of fs/nfft
#' @export
#'
#' @examples
#'# Set up spectrogram parameters
#'res <- process_spectrogram_params(fs, nfft, frequencyRange, window_start, winsize_samples)
#'window_idxs <- res[[1]]
#'stimes <- res[[2]]
#'sfreqs <- res[[3]]
#'freq_inds <- res[[4]]
process_spectrogram_params <- function(fs, nfft, frequencyRange, window_start, datawin_size){
  #
  #
  # Params:
  #         window_start (numeric vector): timestamps representing the beginning time for each window -- required
  #
  # Returns:
  #         window_idxs (matrix): indices of timestamps for each window (nxm where n=number of windows and m=datawin_size)
  #         stimes (numeric vector): times for the centers of the spectral bins (1xt)
  #         sfreqs (numeric vector): frequency bins for spectrogram (1xf)
  #         freq_inds (logical vector): indicates which frequencies are being analyzed in an array of frequencies from 0 to fs
  #                    with steps of fs/nfft


  # Create frequency vector
  df <- fs/nfft
  sfreqs <- seq(0, fs, by=df)

  # Get frequencies for given frequency range
  freq_inds <- (sfreqs >= frequencyRange[1]) & (sfreqs <= frequencyRange[2])
  sfreqs <- sfreqs[freq_inds]

  # Compute times in middle of each spectrum
  window_middle_samples <- window_start + round(datawin_size/2)
  stimes <- (window_middle_samples-1) / fs  # stimes starts from 0

  # Get indices for each window
  window_idxs <- lapply(window_start, window_index_helper, datawin_size=datawin_size) # list of indices for n windows


  return(list(window_idxs, stimes, sfreqs, freq_inds))

}

window_index_helper <- function(start, datawin_size){
  res = seq(start, start+datawin_size-1, by=1)
  return(res)
}



# Display Spectrogram Properties #
display_spectrogram_properties <- function(fs, timeBandwidth, numTapers, data_windowParams, frequencyRange, detrendOpt){
  # Prints spectrogram properties
  #
  # Params:
  #         fs (numeric): sampling frequency in Hz  -- required
  #         timeBandwidth (numeric): time-half bandwidth product (window duration*1/2*frequency_resolution) -- required
  #         numTapers (numeric): number of DPSS tapers to use -- required
  #         data_windowParams (numeric vector): c(window length(s), window step size(s) -- required
  #         frequencyRange (numeric vector): c(<min frequency>, <max frequency>) -- required
  #         detrendOpt (char): detrend data window ('linear' (default), 'constant', 'off')
  #
  # Returns:
  #         This function does not return anythin

  data_windowParams = data_windowParams / fs

  # Print spectrogram properties
  print("Multitaper Spectrogram Properties: ")
  print(paste('     Spectral Resolution: ', toString(2 * timeBandwidth / data_windowParams[1]), 'Hz', sep=""))
  print(paste('     Window Length: ', toString(data_windowParams[1]), 's', sep=""))
  print(paste('     Window Step: ', toString(data_windowParams[2]), 's', sep=""))
  print(paste('     Time Half-Bandwidth Product: ', toString(timeBandwidth), sep=""))
  print(paste('     Number of Tapers: ', toString(numTapers), sep=""))
  print(paste('     Frequency Range: ', toString(frequencyRange[1]), "-", toString(frequencyRange[2]), 'Hz', sep=""))
  print(paste('     Detrend: ', detrendOpt, sep=""))

}


# Convert power to dB #
nanpow2db <- function(y){
  # Power to dB conversion, setting negatives and zeros to NaN
  #
  # params:
  #         y: power --required
  #
  # returns:
  #         ydB: dB (with 0s and negativs set to NaN)

  if(length(y)==1){
    if(y==0){
      return(NaN)
    } else(ydB <- 10*log10(y))
  }else{
    y[y==0] <- NaN
    ydB <- 10*log10(y)
  }
  return(ydB)
}



# Calculate multitpaer spectrum of single segment #
calc_mts_segment <- function(data_segment, dpss_tapers, nfft, freq_inds, weighting, wt, dpss_eigen, numTapers, detrendOpt){
  # Calculate multitaper spectrum for a single segment of data
  #
  # params:
  #         data_segment (numeric vector): segment of the EEG data of length window size (s) * fs -- required
  #         dpss_tapers (numeric matrix): DPSS taper params to multiply signal by. Dims are (numTapers, winsize_samples)
  #                                        -- required
  #         nfft (numeric): length of signal to calculate fft on -- required
  #         freq_inds (logical vector): boolean array indicating frequencies to use in an array of frequenices
  #                                    from 0 to fs with steps of fs/nfft --required
  #         weighting (char): weighting of tapers ('unity' (default), 'eigen', 'adapt') --required
  #         wt (numeric vector or numeric): precomputed taper weights --required
  #         dpss_eigen (numeric vector): --required
  #         numTapers (numeric): number of dpss tapers being used --required
  #         detrendOpt (char): detrend data window ('linear' (default), 'constant', 'off') --required
  #
  # returns:
  #         mt_spectrum (numeric matrix): spectral power for single window

  library(pracma)

  # If segment has all zeros, return vector of zeros
  if(all(data_segment==0)){
    ret <- rep(0, sum(freq_inds))
    return(ret)
  }

  # Optionally detrend data to remove low freq DC component
  if(detrendOpt != 'off'){
    data_segment <- detrend(data_segment, tt=detrendOpt)
  }

  # Multiply data by dpss tapers (STEP 2)
  tapered_data <- sweep(dpss_tapers, 1, data_segment, '*')


  # Manually add nfft zero-padding (R's fft function does not support)
  tapered_padded_data <- rbind(tapered_data, matrix(0, nrow=nfft-nrow(tapered_data), ncol=ncol(tapered_data)))


  # Compute the FFT (STEP 3)
  fft_data <- apply(tapered_padded_data, 2, fft)
  # Compute the weighted mean spectral power across tapers (STEP 4)
  Spower = Im(fft_data)^2 + Re(fft_data)^2;
  if(weighting == 'adapt'){
    # daptive weights - for colored noise spectrum (Percival & Walden p368-p370)
    x = matrix(data_segment, nrow=1)
    Tpower <- x %*% (t(x)/length(x))
    Spower_iter <- rowMeans(Spower[,1:2])
    Spower_iter <- matrix(Spower_iter, ncol=1) # (nfft,1)
    a <- (1 - dpss_eigen) * as.vector(Tpower)
    a <- matrix(a, ncol=1) # (numTapers, 1)
    dpss_eigen = matrix(dpss_eigen, nrow=1)
    for(i in 1:3){ # run 3 iterations
      # calculate the MSE weights
      b = (Spower_iter %*% calibrate::ones(1,numTapers)) / ( (Spower_iter %*% dpss_eigen) + repmat(t(a),nfft,1) )
      # calculate new spectral estimate
      wk = b^2 * (calibrate::ones(nfft,1) %*% dpss_eigen)
      Spower_iter = matrix(colSums(t(wk) * t(Spower)), nrow=1) / rowSums(wk)
      Spower_iter = matrix(Spower_iter, ncol=1)
    }
    mt_spectrum = as.vector(Spower_iter)
  } else{
    # eigenvalue or uniform weights
    mt_spectrum = Spower %*% wt
    mt_spectrum = as.vector(mt_spectrum)
  }

  return(mt_spectrum[freq_inds])
}
