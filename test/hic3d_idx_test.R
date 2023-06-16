pos_idx <- function(x, y, z){
    idx <- x*y*z - z*x*(x+1)/2 - x*y*(y+1)/2 + x*(x+1)*(x+2)/6
    return(idx)
}

hic3d_brute_force <- function(x, y, z, N){ # N : max length
  it = 0
  for(i in 1:x){
    for(j in (i+1):(N-1)){
      for(k in (j+1):N){
        it = it+1  # index counter
        if (i==x && j==y && k==z){
          return(it)
        }
      }
    }
  }
}

hic3d_idx <- function(x, y, z, N){
  sum_xNN <- (x-1)*(N-1)*N - N*(x-1)*x/2 - (x-1)*(N-1)*N/2 + x*(x-1)*(x+1)/6 # = pos_idx(x-1,N-1,N)
  #sum_yz <- z*(y-x) - y*(y+1)/2 + x*(x+1)/2
  sum_yN <- (y-x-1)*N -((y-1)*y - (x+1)*x)/2
  sum_z <-  (z-y)
  #print(paste(sum_xNN,",",sum_yN,",",sum_z))
  #print(paste('x-1',hic3d_idx(x-1,N,N,N), sum_yz))
  # pos_idx((x-1),N-1, N)
  idx <- sum_xNN + sum_yN + sum_z
  return(idx)
}

test_hic3d_idx <- function(N){
  print(paste("running test_hic3d_idx(N=",N,")..."))
  idx <- 1
  start_time <- Sys.time()
  for(i in 1:(N-2)){
    for(j in (i+1):(N-1)){
      for(k in (j+1):N) {
        # print(paste("idx: ", idx))
        hi <- hic3d_idx(i,j,k,N)
        if (hi != idx){
          print(paste("error hi(",i,",",j,",",k,")=", hi," != ", idx))
          traceback()
          # exit(0)
        } else {
        #  print(paste("OK hi(",i,",",j,",",k,")=", hi," == ", idx))
        }
        idx <- idx+1
      }
    }
  }
  print("Test OK")
}

benchmark_hic3d_idx <- function(N){
  print(paste("running benchmark_hic3d_idx(N=",N,")..."))
  start_time <- Sys.time()
  for(i in 1:(N-2)){
    for(j in (i+1):(N-1)){
      for(k in (j+1):N) {
        hic3d_idx(i,j,k,N)
      }
    }
  }
  end_time <- Sys.time()
  time_lapse = difftime(end_time, start_time) #, units="secs")
  return(time_lapse)
}

benchmark_hic3d_brute_force <- function(N){
  print(paste("running benchmark_hic3d_brute_force(N=",N,")..."))
  start_time <- Sys.time()
  for(i in 1:(N-2)){
    for(j in (i+1):(N-1)){
      for(k in (j+1):N) {
        hic3d_brute_force(i,j,k,N)
      }
    }
  }
  end_time <- Sys.time()
  time_lapse = difftime(end_time, start_time) #, units="secs")
  return(time_lapse)
}

test_hic3d_idx(2024)

hic3d_idx_time = benchmark_hic3d_idx(100)
hic3d_idx_time

hic3d_brute_force_time = benchmark_hic3d_brute_force(100)
hic3d_brute_force_time

performance_gain = as.numeric(hic3d_brute_force_time, units="secs") / as.numeric(hic3d_idx_time, units="secs")
print(paste("Performance gain: ", performance_gain, "times"))
