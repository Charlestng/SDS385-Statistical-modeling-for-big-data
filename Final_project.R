#### SDS 385 final project

library(FKF)
library(RJSONIO)
library(sp)
library(ggplot2)
library(ggmap)

# Import the GPS tracking data set
GPS_data2 <- read.csv("C:/Users/Charles/Documents/McCombs/Statistical modeling for big data/final project/campus_gammaradiation.csv")

# sort by timestamp
GPS_data <- GPS_data[order(GPS_data$timestamp),]

# create a variable to have the date in numeric format
GPS_data$time_POSIXct <- as.POSIXct(GPS_data[,"timestamp"], format="%Y-%m-%d  %H:%M:%S")
GPS_data$time_num <- as.numeric(GPS_data$time_POSIXct)
n_obs <- nrow(GPS_data)

# create a time_lapse variable
GPS_data$time_lapse <- 0
GPS_data[2:n_obs,"time_lapse"] <- GPS_data[2:n_obs,"time_num"]-GPS_data[1:(n_obs-1),"time_num"]

# Calculate state vector
# x = longitude
# y = latitude
# x' = (x_k+1 - x_k) / (timestamp_k+1 - timestamp_k)
GPS_data$delta_lon <- 0
GPS_data[2:n_obs,"delta_lon"] <- GPS_data[2:n_obs,"lon"]-GPS_data[1:(n_obs-1),"lon"]

GPS_data$V_lon <- 0
GPS_data[2:n_obs,"V_lon"] <- GPS_data[2:n_obs,"delta_lon"] / GPS_data[2:n_obs,"time_lapse"]

# y' = (y_k+1 - y_k) / (timestamp_k+1 - timestamp_k)
GPS_data$delta_lat <- 0
GPS_data[2:n_obs,"delta_lat"] <- GPS_data[2:n_obs,"lat"]-GPS_data[1:(n_obs-1),"lat"]

GPS_data$V_lat <- 0
GPS_data[2:n_obs,"V_lat"] <- GPS_data[2:n_obs,"delta_lat"] / GPS_data[2:n_obs,"time_lapse"]

# Calculate speed
GPS_data$speed <- 0
GPS_data[2:n_obs,"speed"] <- sqrt(GPS_data[2:n_obs, "V_lon"]^2
                                  + GPS_data[2:n_obs, "V_lat"]^2)

# Calculate theta, the angle between the orientation of the car and the X axis
GPS_data$theta <- 0
GPS_data[2:n_obs, "theta"] <- sign(GPS_data[2:n_obs, "delta_lat"]) *
                              GPS_data[2:n_obs, "delta_lon"] /
                              sqrt(GPS_data[2:n_obs, "delta_lon"]^2
                                    + GPS_data[2:n_obs, "delta_lat"]^2)

# Correct where theta is not defined because no movement between two observations
index_theta <- which(GPS_data$delta_lon==0 & GPS_data$delta_lat==0)
index_theta <- index_theta[2:length(index_theta)]
theta <- GPS_data[,"theta"]
for (i in index_theta){
  print(i)
  theta[i] <- theta[i-1]
}
GPS_data[,"theta"] <- theta

# Calculate cos of theta, the angle between the orientation of the car and the X axis
GPS_data$cos_theta <- 0
GPS_data[2:n_obs, "cos_theta"] <- GPS_data[2:n_obs, "delta_lon"] /
                              sqrt(GPS_data[2:n_obs, "delta_lon"]^2
                              + GPS_data[2:n_obs, "delta_lat"]^2)

# Calculate sin of theta, the angle between the orientation of the car and the X axis
GPS_data$sin_theta <- 0
GPS_data[2:n_obs, "sin_theta"] <- GPS_data[2:n_obs, "delta_lat"] /
                                sqrt(GPS_data[2:n_obs, "delta_lon"]^2
                                + GPS_data[2:n_obs, "delta_lat"]^2)


# The police car drives in straight lines by block
summary(GPS_data$lat)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#29.96   30.28   30.28   30.29   30.28   30.69 
summary(GPS_data$lon)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-99.28  -97.73  -97.73  -97.73  -97.73  -97.54 
GPS_data[1,"timestamp"]
# 2015-06-11 16:34:21
GPS_data[n_obs,"timestamp"]
# 2015-11-18 22:33:32

# The coordinates of campus limits are
# lat: 30.27, 30.30, lon : -97.72, -97.76
GPS_data_campus <- GPS_data[GPS_data$lat < 30.30 &
                              GPS_data$lat > 30.27 &
                              GPS_data$lon < -97.72 &
                              GPS_data$lon > -97.76,]


# Split the data into portions of measurement, 99.7% of the timelapses are under 60s
# I will therefore consider that data points separated by more than 60s belong to different
# samples

# Indicate where two measurements are separated by more than 60 seconds
GPS_data_campus$morethan60 <- GPS_data_campus$time_lapse > 60

# Define function to create a vector that indicates for each observation the group number
# it belongs to
# Fuction : split_timelapse
# input:
#    data: a data frame containing a time serie to be split in different groups
#    time_lapse: the time lapse variable in data
#    truth: the indicator of a data separating two groups
# output:
#    portion: a vector of length = number of observations in data containing the group index
split_timelapse <- function(data, time_lapse, truth){
  data$portion_index <- 0
  portion <- c()
  begin <- 1
  num <- 1
  change_portion <- which(data[,truth])
  for (i in c(change_portion,nrow(data))){
    print(i)
    mat_portion <- matrix(num, nrow = i - begin)
    portion <- c(portion, mat_portion)
    begin <- i
    num <- num +1
  }
  portion <- c(portion, matrix(num, nrow = 1))
  portion <- rbind(portion)
  return(portion)
}

# create a vector to split data groups separated by more than 60 seconds
portions_vector <- split_timelapse(GPS_data_campus, "time_lapse", "morethan60")

# add the grouping variable to the data
GPS_data_campus$portions <- t(portions_vector)

# The kalman filtering technique to estimate a model of the form
# x_k = F_k-1 x_k-1 + B_k u_k + w_k
# x_k is of the form t(x, x', y, y')
# x_k is the unobserved true state of the system
# F_k state transition model
# We assume that the movement of the car is going to be rectilinear by block, because most roads
# in cities are straight, each F is going to be a matrix of the form
# 1    time_lapse    0    0
# 0    1             0    0
# 0    0             1    time_lapse
# 0    0             0    1
# We can pre-cache these matrices for efficiency gains in the code too

# B_k control input model
# Here we don't control any input to the system so B_k = 0

# w_k process noise drawn from a multivariate normal distribution with zero mean and 
#     covariance Q_k
# z_k = H_k x_k + v_k
# Since we only truly measure the lon and lat for each k and not the speed
# z_k is the observed state of the system
# H_k is the observation model 
# H_k will be of the form
# 1    0    0    0
# 0    0    1    0
# v_k zero mean gaussian white noise with covariance R_k

# The Kalman filtering model is estimated in two steps
# Predict
#
# Predicted (a priori) state estimate 
#    xhat_k|k-1 = F_k xhat_k-1|k-1 + B_k-1 u_k-1
# Predicted (a priori) estimate covariance
#    P_k|k-1 = F_k P_k-1|k-1 t(F_k) + Q_k
#
# Update
#
# Innovation or measurement residual
#    ytilde_k = z_k - H_k xhat_k|k-1
# Innovation (or residual) covariance
#    S_k = H_k P_k|k-1 t(H_k) + R_k
# Optimal Kalman gain	
#    K_k = P_k|k-1 t(H_k) inv(S_k)
# Updated (a posteriori) state estimate	
#    xhat_k|k = xhat_k|k-1 + K_k ytilde_k
# Updated (a posteriori) estimate covariance	
#    P_k|k = (I - K_k H_k) P_k|k-1

# the measurement noise variance is
GPS_for_noise <- GPS_data[100:350,c("lon","V_lon","lat","V_lat","timestamp")]
GPS_for_noise$lon <- GPS_for_noise$lon - mean(GPS_for_noise$lon)
GPS_for_noise$V_lon <- GPS_for_noise$V_lon - mean(GPS_for_noise$V_lon)
GPS_for_noise$lat <- GPS_for_noise$lat - mean(GPS_for_noise$lat)
GPS_for_noise$V_lat <- GPS_for_noise$V_lat - mean(GPS_for_noise$V_lat)

R_k <-  var(GPS_for_noise[,c("lon","lat")])
sigma_a <- var(GPS_for_noise$lat)

# the kalman filter to treat everything at the same time, but the model here is too simplistic
# to model the movement of the car correctly
Kalman_filter <- function(data, x_pos, y_pos, x_speed, y_speed, time_lapse){
  # convert to numeric format
  for (i in 1:ncol(data)){
    data[,i] <- as.numeric(data[,i])
  }
  
  # Initialize starting point of the car
  x_0 <- data[1,c(x_pos,x_speed,y_pos,y_speed)]
  
  # Number of observations
  n_obs <- nrow(data)
  
  # bookeeping variable for the estimate x_hat
  x_hat <- t(matrix(x_0))
  
  # Initialise the covariance matrix P
  P_k <- matrix(data = 0, nrow = 4, ncol = 4)
  
  # Initialise H_k
  H_k <- matrix(c(1,0,0,0,
                  0,0,1,0), nrow = 2)
  
  # Identity matrix
  Id <- diag(x = 1, nrow = 4, ncol = 4)
  
  # define F_k
  F_k <- t(matrix(data = c(1, 3,0,0,
                           0,1,0,0,
                           0,0,1,3,
                           0,0,0,1),
                  nrow = 4,
                  ncol = 4))
  
  # create Q_k based on time lapse
  G_k <- matrix(c(3^2/2,
                  3,
                  3^2/2,
                  3),
                nrow = 4)
  
  # Loop over the observations to estimate the filtered time serie
  for (k in 2:n_obs){
    print(k)
    # Predict
    # create F_k
    #F_k <- t(matrix(data = c(1, data[k,time_lapse],0,0,
    #                       0,1,0,0,
    #                       0,0,1,data[k,time_lapse],
    #                       0,0,0,1),
    #              nrow = 4,
    #              ncol = 4))
    
    
    # Predicted (a priori) state estimate 
    x_hat_k <- F_k %*% as.numeric(x_hat[k-1,])
    
    # create Q_k based on time lapse
    #G_k <- matrix(c(data[k,time_lapse]^2/2,
    #                data[k,time_lapse],
    #                data[k,time_lapse]^2/2,
    #                data[k,time_lapse]),
    #              nrow = 4)
    Q_k <- G_k %*% t(G_k) * sigma_a
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
    x_hat <- rbind(x_hat, t(x_hat_k))
  }
  return(data.frame(x_hat))
}


# Function Kalman filtering for one-dimension model
# input:
#    data: a data frame containing the time serie to be denoised
#    pos: the position variable
#    speed: the speed variable
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
  sigma_p <- var(GPS_for_noise[,pos])
  
  # Create sigma_v
  sigma_v <- var(GPS_for_noise[,speed])
  
  # Initialise the covariance matrix P
  P_k <- matrix(data = c(sigma_a,0,0,sigma_v), nrow = 2, ncol = 2)
  
  # Initialise H_k
  H_k <- matrix(c(1,0), ncol = 2)
  
  # Identity matrix
  Id <- diag(x = 1, nrow = 2, ncol = 2)
  
  # Create R_k
  R_k <-  var(GPS_for_noise[,pos])
  
  # Create sigma_a
  sigma_a <- var(GPS_for_noise[,pos])/100
  
  
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
    # create Q_k based on time lapse
    G_k <- matrix(c(data[k,time_lapse]^2/2,
                    data[k,time_lapse]),
                  nrow = 2)
    
    Q_k <- G_k %*% t(G_k) * sigma_a
    # Predicted (a priori) estimate covariance
    P_k <- F_k %*% P_k %*% t(F_k) + Q_k
    
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

colnames(test_lat) <- c("lat", "speed_lat", "time_num")
colnames(test_lon) <- c("lon", "speed_lon", "time_num")
test <- cbind(test_lat,test_lon)

lionmap <- get_map(location = c(-97.733, 30.288), zoom = 15, maptype = "hybrid")

ggmap(lionmap) + geom_path(data = GPS_data_campus[GPS_data_campus$portion == 389,
                                                  c("lat","lon")]
                           , aes(x = lon, y = lat), 
                           size = 1, color = "red", type="l")
ggmap(lionmap) + geom_point(data = test[,c("lat","lon")]
           , aes(x = lon, y = lat), 
           size = 1, color = "blue", type="l")

# Function Kalman filtering for one-dimension model
# input:
#    data: a data frame containing the time serie to be denoised
#    pos: the position variable
#    speed: the speed variable
#    time_lapse: the variable indicating the time lapse between two observations
#    time_stamp: the variable indicating the date of the observation
# output:
#    x_hat: a data frame containing the smoothed position, smoothed speed, and time_stamp
Kalman_filter <- function(data, x_pos, y_pos, theta, speed, time_lapse, time_stamp){
  # convert to numeric format
  for (i in 1:ncol(data)){
    data[,i] <- as.numeric(data[,i])
  }
  
  # Initialize starting point of the car
  x_0 <- data[1,c(x_pos,y_pos,theta,speed,time_stamp)]
  
  # Number of observations
  n_obs <- nrow(data)
  
  # bookeeping variable for the estimate x_hat
  x_hat <- t(matrix(x_0))
  
  # Initialise the covariance matrix P
  P_k <- diag( c(3.7, 6.4, 3.7, 6.7), nrow = 4, ncol = 4)
  
  # Initialise H_k
  H_k <- diag(x = 1, nrow = 4, ncol = 4)
  
  # Identity matrix
  Id <- diag(x = 1, nrow = 4, ncol = 4)
  
  # R_k
  R_k <- diag(c(0.23,0.26,0.01,1.05), nrow = 4, ncol = 4)/100000
  
  # Q_k
  Q_k <- diag(c(1.51,5.58,1.95,1.68), nrow = 4, ncol = 4)/100000
  
  # error bookeping
  error <- c(0)
  
  # Loop over the observations to estimate the filtered time serie
  for (k in 2:n_obs){
    print(k)
    # Predict
    # create F_k
    #F_k <- t(matrix(data = c(1, data[k,time_lapse],0,0,
    #                       0,1,0,0,
    #                       0,0,1,data[k,time_lapse],
    #                       0,0,0,1),
    #              nrow = 4,
    #              ncol = 4))
    
    # precache the time_lapse
    delta_t <- data[k+1,time_lapse]
    
    # transform F_k according to the filtered x_hat_k
    F_k <- matrix(c(1,0,0,0,
                    0,1,0,0,
                    0,0,1,0,
                    delta_t * cos(as.numeric(x_hat[k-1,3])),
                    delta_t * sin(as.numeric(x_hat[k-1,3])), 0, 1),
                  nrow = 4,
                  ncol = 4)
    
    # adapt Q_k and R_k with delta t
    #R_k <- diag(c(0.23,0.26,0.01,1.05), nrow = 4, ncol = 4) * 100 * delta_t
    #Q_k <- diag(c(1.51,5.58,1.95,1.68), nrow = 4, ncol = 4) * 100 * delta_t
    
    
    # Predicted (a priori) state estimate 
    x_hat_k <- F_k %*% as.numeric(x_hat[k-1,1:4])
    
    # create Q_k based on time lapse
    
    # Predicted (a priori) estimate covariance
    P_k <- F_k %*% P_k %*% t(F_k) + Q_k
    
    # Update
    # Innovation or measurement residual
    ytilde_k = data[k, c(x_pos,y_pos,theta,speed)] - H_k %*% x_hat_k
    error <- c(error, ytilde_k)
    
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
  x_hat <- cbind(x_hat,error)
  colnames(x_hat) <- c("lon", "lat", "theta", "speed", "timestamp","error")
  return(x_hat)
}

test_2D <- Kalman_filter(data = GPS_data_campus[GPS_data_campus$portions==389,],
                             x_pos = "lon", y_pos = "lat", theta = "theta", speed = "speed",
                             time_lapse = "time_lapse", time_stamp = "timestamp")

plot(GPS_data_campus[GPS_data_campus$portions==389,]$lon,
     GPS_data_campus[GPS_data_campus$portions==389,]$lat, type="l")
lines(test_2D$lon,test_2D$lat, col = "red")

portion10 <- GPS_data_campus[GPS_data_campus$portions==389,]

# What the R package function for kalman filter looks like

Initialization:
  if(i == 1){
    at[, i] = a0
    Pt[,, i] = P0
  }
Updating equations:
  vt[, i] = yt[, i] - ct[, i] - Zt[,,i] %*% at[, i]
Ft[,, i] = Zt[,, i] %*% Pt[,, i] %*% t(Zt[,, i]) + GGt[,, i]
Kt[,, i] = Pt[,, i] %*% t(Zt[,, i]) %*% solve(Ft[,, i])
att[, i] = at[, i] + Kt[,, i] %*% vt[, i]
Ptt[, i] = Pt[,, i] - Pt[,, i] %*% t(Zt[,, i]) %*% t(Kt[,, i])
Prediction equations:
  at[, i + 1] = dt[, i] + Tt[,, i] %*% att[, i]
Pt[,, i + 1] = Tt[,, i] %*% Ptt[,, i] %*% t(Tt[,, i]) + HHt[,, i]

plot(GPS_data[100:350,c("lon","lat")])
 var(GPS_data[100:350,c("lon","V_lon","lat","V_lat")])
#         lon         V_lon           lat         V_lat
#lon    1.327926e-09  5.040236e-11 -4.663612e-10 -1.568053e-12
#V_lon  5.040236e-11  4.247126e-11 -1.731495e-11 -9.126410e-12
#lat   -4.663612e-10 -1.731495e-11  3.968685e-10  8.948400e-12
#V_lat -1.568053e-12 -9.126410e-12  8.948400e-12  8.627073e-12
 
 
 
 lionmap <- get_map(location = c(-97.733, 30.288), zoom = 15, maptype = "hybrid")
 
 ggmap(lionmap) + geom_path(data = GPS_data_campus[GPS_data_campus$portion == 389,
                                                   c("lat","lon")]
                             , aes(x = lon, y = lat), 
                                   size = 3, color = "red", type="l")
