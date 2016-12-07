#### SDS 385 final project

### DATA processing

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
