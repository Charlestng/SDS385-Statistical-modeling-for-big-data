# Function Kalman filtering for one-dimension model
# input:
#    data: a data frame containing the time serie to be denoised
#    pos: the position variable
#    speed: the speed variable -- not used in estimation
#    time_lapse: the variable indicating the time lapse between two observations
#    time_stamp: the variable indicating the date of the observation
# output:
#    x_hat: a data frame containing the smoothed position, smoothed speed, and time_stamp
Kalman_filter_1d <- function(data, pos, speed, time_lapse, time_stamp){
  # convert to numeric format
  for (i in 1:ncol(data)){
    data[,i] <- as.numeric(data[,i])
  }
  
  # Initialize starting point of the car
  x_0 <- data[1,c(pos,speed,time_stamp)]
  
  # Number of observations
  n_obs <- nrow(data)
  
  # bookeeping variable for the estimate x_hat
  x_hat <- t(matrix(x_0))
  
  # Create sigma_a
  sigma_p <- 5*10E-10
  
  # Create sigma_v
  sigma_v <- 10E-11
  
  # Initialise the covariance matrix P
  P_k <- matrix(data = c(sigma_a,0,0,sigma_v), nrow = 2, ncol = 2)
  
  # Initialise H_k
  H_k <- matrix(c(1,0), ncol = 2)
  
  # Identity matrix
  Id <- diag(x = 1, nrow = 2, ncol = 2)
  
  # Create R_k
  R_k <-  5*10E-7
  
  # Create sigma_a
  sigma_a <- 5*10E-10/100
  
  
  # Loop over the observations to estimate the filtered time serie
  for (k in 2:n_obs){
    print(k)
    # Predict
    # define F_k
    F_k <- t(matrix(data = c(1, data[k,time_lapse],0,1),
                    nrow = 2,
                    ncol = 2))
    
    
    # Predicted (a priori) state estimate 
    x_hat_k <- F_k %*% as.numeric(x_hat[k-1,1:2])
    
    # create Q_k based on time lapse
    delta_t <- data[k,time_lapse]
    G_k <- matrix(c(delta_t^2/2,
                    delta_t),
                  nrow = 2)
    
    Q_k <- G_k %*% t(G_k) * sigma_a
    # Predicted (a priori) estimate covariance
    P_k <- F_k %*% P_k %*% t(F_k) + Q_k * delta_t
    
    # Update
    # Innovation or measurement residual
    ytilde_k = data[k, pos] - H_k %*% x_hat_k
    
    # Innovation (or residual) covariance
    S_k = H_k %*% P_k %*% t(H_k) + R_k
    
    # Optimal Kalman gain	
    K_k = P_k %*% t(H_k) %*% solve(S_k)
    
    # Updated (a posteriori) state estimate	
    x_hat_k = x_hat_k + K_k %*% as.numeric(ytilde_k)
    
    # Updated (a posteriori) estimate covariance	
    P_k = (Id - K_k %*% H_k) %*% P_k
    
    # add the new x_hat_k to the filtered time serie
    x_hat <- rbind(x_hat, c(t(x_hat_k),data[k,time_stamp]))
  }
  return(data.frame(x_hat))
}

test_lon <- Kalman_filter_1d(data = GPS_data_campus[GPS_data_campus$portions==389,],
                             pos = "lon", speed = "V_lon",
                             time_lapse = "time_lapse", time_stamp = "timestamp")
test_lat <- Kalman_filter_1d(data = GPS_data_campus[GPS_data_campus$portions==389,],
                             pos = "lat", speed = "V_lat",
                             time_lapse = "time_lapse", time_stamp = "timestamp")



plot(test_lon$X3,test_lon$X1, type = "l", col ="red")
lines(GPS_data_campus[GPS_data_campus$portions==389,]$timestamp,
      GPS_data_campus[GPS_data_campus$portions==389,]$lon)
plot(test_lat$X3,test_lat$X1, type = "l", col ="red")
lines(GPS_data_campus[GPS_data_campus$portions==389,]$timestamp,
      GPS_data_campus[GPS_data_campus$portions==389,]$lat)

plot(GPS_data_campus[GPS_data_campus$portions==389,]$lon,
     GPS_data_campus[GPS_data_campus$portions==389,]$lat, type="l")
lines(test_lon$X1,test_lat$X1, col = "red")
