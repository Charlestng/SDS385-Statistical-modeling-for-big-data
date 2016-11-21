#### SDS 385 final project

library("FKF")

# Import the GPS tracking data set
GPS_data <- read.csv("C:/Users/Charles/Documents/McCombs/Statistical modeling for big data/final project/campus_gammaradiation.csv")

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
GPS_for_noise <- GPS_data[100:350,c("lon","V_lon","lat","V_lat")]
GPS_for_noise$lon <- GPS_for_noise$lon - mean(GPS_for_noise$lon)
GPS_for_noise$V_lon <- GPS_for_noise$V_lon - mean(GPS_for_noise$V_lon)
GPS_for_noise$lat <- GPS_for_noise$lat - mean(GPS_for_noise$lat)
GPS_for_noise$V_lat <- GPS_for_noise$V_lat - mean(GPS_for_noise$V_lat)

R_k <-  var(GPS_for_noise[,c("lon","V_lon","lat","V_lat")])

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
  H_k <- matrix(1, nrow = 4)
  
  # Identity matrix
  Id <- diag(x = 1, nrow = 4, ncol = 4)
  
  # Loop over the observations to estimate the filtered time serie
  for (k in 2:n_obs){
    # Predict
    # create F_k
    F_k <- t(matrix(data = c(1, data[k,time_lapse],0,0,
                           0,1,0,0,
                           0,0,1,data[k,time_lapse],
                           0,0,0,1),
                  nrow = 4,
                  ncol = 4))
    # Predicted (a priori) state estimate 
    print("F_k")
    print(F_k)
    print("x_hat")
    print(x_hat[k-1,])
    
    x_hat_k <- F_k %*% as.numeric(x_hat[k-1,])
    
    # create Q_k based on time lapse
    G_k <- matrix(c(data[k,time_lapse]^2/2,
                    data[k,time_lapse],
                    data[k,time_lapse]^2/2,
                    data[k,time_lapse]),
                  nrow = 4)
    Q_k <- G_k %*% t(G_k)
    # Predicted (a priori) estimate covariance
    P_k <- F_k %*% P_k %*% t(F_k) + Q_k
    
    # Update
    # Innovation or measurement residual
    ytilde_k = data[k, c(x_pos,x_speed,y_pos,y_speed)] - x_hat_k
    
    # Innovation (or residual) covariance
    S_k = P_k + R_k
    print(S_k)
    
    # Optimal Kalman gain	
    K_k = P_k %*% solve(S_k)
    
    # Updated (a posteriori) state estimate	
    x_hat_k = x_hat_k + K_k %*% as.numeric(ytilde_k)
    
    # Updated (a posteriori) estimate covariance	
    P_k = (Id - K_k) %*% P_k
    
    # add the new x_hat_k to the filtered time serie
    x_hat <- rbind(x_hat, t(x_hat_k))
  }
  return(x_hat)
}

# I need to convert the variables lon and lat into a distance from an origin somewhere
# because the scale of Q_k is so much bigger than R_k that S_k is not invertible


test <- Kalman_filter(data = GPS_data_campus[1:1000,], x_pos = "lon", y_pos = "lat",
                      x_speed = "V_lon", y_speed = "V_lat", time_lapse = "time_lapse")

F_k <- t(matrix(data = c(1, GPS_data_campus[2,"time_lapse"],0,0,
                         0,1,0,0,
                         0,0,1,GPS_data_campus[2,"time_lapse"],
                         0,0,0,1),
                nrow = 4,
                ncol = 4))
x_0 <-matrix(GPS_data_campus[1,c("lon","V_lon","lat","V_lat")], ncol = 4)

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
