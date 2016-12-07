# Function Kalman filtering for one-dimension model
# input:
#    data: a data frame containing the time serie to be denoised
#    x_pos: the position variable along the X axis (longitude)
#    y_pos: the position variable along the Y axis (latitude)
#    x_speed: the speed variable along the X axis (longitude) -- not used for estimation
#    y_speed: the speed variable along the Y axis (latitude) -- not used for estimation
#    time_lapse: the variable indicating the time lapse between two observations
#    time_stamp: the variable indicating the date of the observation
# output:
#    x_hat: a data frame containing the smoothed position, smoothed speed, and time_stamp

Kalman_filter_cartesian <- function(data, x_pos, y_pos, x_speed, y_speed, time_lapse, time_stamp){
  # convert to numeric format
  for (i in 1:ncol(data)){
    data[,i] <- as.numeric(data[,i])
  }
  
  # Initialize starting point of the car
  x_0 <- data[1,c(x_pos,y_pos,x_speed,y_speed,time_stamp)]
  
  # Number of observations
  n_obs <- nrow(data)
  
  # bookeeping variable for the estimate x_hat
  x_hat <- t(matrix(x_0))
  
  # Initialise the covariance matrix P
  P_k <- diag( c(3.7, 3.7, 3.7, 3.7), nrow = 4, ncol = 4)
  
  # Initialise H_k
  H_k <- matrix(data = c(1,0,
                         0,1,
                         0,0,
                         0,0), nrow = 2, ncol = 4)
  
  # Identity matrix
  Id <- diag(x = 1, nrow = 4, ncol = 4)
  
  # Loop over the observations to estimate the filtered time serie
  for (k in 2:n_obs){
    print(k)
    
    delta_t <- data[k,time_lapse]
    
    # transform F_k according to the filtered x_hat_k
    F_k <- matrix(c(1,0,0,0,
                    0,1,0,0,
                    delta_t,0,1,0,
                    0,delta_t, 0, 1),
                  nrow = 4,
                  ncol = 4)
    
    # adapt Q_k and R_k with delta t
    R_k <- diag(c(0.23,0.26), nrow = 2, ncol = 2) * 10E-4 * delta_t
    Q_k <- diag(c(1.51,1.95,1.51,1.95), nrow = 4, ncol = 4) * 10E-8 * delta_t
    
    
    # Predicted (a priori) state estimate 
    x_hat_k <- F_k %*% as.numeric(x_hat[k-1,1:4])
    
    # create Q_k based on time lapse
    
    # Predicted (a priori) estimate covariance
    P_k <- F_k %*% P_k %*% t(F_k) + Q_k
    
    # Update
    # Innovation or measurement residual
    ytilde_k = data[k, c(x_pos,y_pos)] - H_k %*% x_hat_k
    
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
  x_hat <- data.frame(x_hat)
  colnames(x_hat) <- c("lon", "lat", "V_lon", "V_lat", "timestamp")
  return(x_hat)
}

test_cart <- Kalman_filter_cartesian(data = GPS_data_campus[GPS_data_campus$portions==389,],
                                     x_pos = "lon", y_pos = "lat", x_speed = "V_lon", y_speed = "V_lat",
                                     time_lapse = "time_lapse", time_stamp = "timestamp")

plot(GPS_data_campus[GPS_data_campus$portions==389,]$lon,
     GPS_data_campus[GPS_data_campus$portions==389,]$lat, type="p")
lines(test_cart$lon,test_cart$lat, col = "red")
