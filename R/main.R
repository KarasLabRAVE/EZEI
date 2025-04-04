calcEI<- function(ieegts){

  # Parameters


  time_bandwidth = 3  # Set time-half bandwidth
  num_tapers = 5  # Set number of tapers (optimal is time_bandwidth*2 - 1)
  window_params = c(1, 0.2)  # Window size is 1s with step size of 0.2s
  min_nfft = 0  # No minimum nfft
  weighting = 'unity'  # weight each taper at 1
  detrend_opt = 'off'  # detrend each window by subtracting the average
  parallel = TRUE  # use multiprocessing
  num_workers = 3  # use 3 cores in multiprocessing
  plot_on = FALSE  # plot spectrogram
  verbose = FALSE  # print extra info
  xyflip = FALSE  # do not transpose spect output matrix


  # Set spectrogram params
  frequency_range = c(0.5, 250)  # Frequency Band
  thetaband<-c(3.5,7.4)
  alphaband<-c(7.4,12.4)
  betaband<-c(12.4,24) #[13-30]Hz
  gammaband<-c(24,140) #[30-90]Hz


  # This number sort of serves as a built-in "time cost". High-frequency
  # activity needs to be robust enough to surpass this v number and still
  # overcome the threshold.
  # It's arbitrary but I think 0.5 is what Bartolomei used.
  v=0.5

  # Lambda also serves as a somewhat arbitrary threshold value.
  lambda=15


  dimieeg=dim(ieegts)
  nt=dimieeg[1]
  nElec=dimieeg[2]

  fs=1000

  nwt=floor((nt/fs-window_params[1])/window_params[2])+1

  data   <- vector(mode="numeric", length=nt)
  data[1:nt]<-ieegts[1:nt,1]
  # Compute the multitaper spectrogram
  results = multitaper_spectrogram_R(data, fs, frequency_range, time_bandwidth, num_tapers, window_params, min_nfft, weighting, detrend_opt, parallel, num_workers,
                                     plot_on, verbose, xyflip)



}
