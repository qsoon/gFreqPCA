# may be redundant 

library(RColorBrewer)
library(colourvalues)
library(grDevices)
# library(SDMTools)
library(network)
library(ggraph)
library(tidygraph)
library(dplyr)
library(gasper)
library(readxl)
library(forecast)
library(ggfortify)
library(Metrics)
library(GNAR)
library(DTWBI)
library(vars)
library(geosphere)
library(xlsx)
library(scales)
library(igraph)
library(pracma)
library(R.matlab)
library(geosphere)
library(grid)
library(gridBase)
library(gridExtra)
library(expm)
library(Hmisc)



###################
#### load data ####
###################

# station Korean, English name
station.name <- read.csv("Data/seoulmetro/station_name.csv", 
                         header=TRUE, fileEncoding = "euc-kr")
colnames(station.name) <- c("code", "name_kor", "name", "line", "external_code")


# station latitude, longitude info
station.loc <- read.csv("Data/seoulmetro/station_location.csv", 
                        header=TRUE, fileEncoding = "euc-kr")
colnames(station.loc) <- c("ID", "name_kor", "line", "lon", "lat")

station.loc2 <- read.csv("Data/seoulmetro/station_location2.csv", 
                         header=TRUE, fileEncoding = "euc-kr")
station.loc2 <- station.loc2[,c(2:6)]
colnames(station.loc2) <- c("line", "station_num", "name_kor", "lon", "lat")


# distance between stations
station.distance <- read.csv("Data/seoulmetro/station_distance.csv", 
                             header=TRUE, fileEncoding = "euc-kr")

station.distance <- station.distance[, c(2:5)]
colnames(station.distance) <- c("line", "name_kor", "btwdist", "cumdist")

# 2021 hourly getting on/off info for each station
hourly.pplnum.2021 <- read.csv("Data/seoulmetro/hourly_pplnum_2021.csv", 
                               header=TRUE, fileEncoding = "euc-kr")
hourly.pplnum.2021 <- hourly.pplnum.2021[, -1]
colnames(hourly.pplnum.2021) <- c("date", "line", "station_num", "name_kor", "type",
                                  paste(rep("t",19), c(1:19), sep=""))


# We have getting on/off info only for line 1~8, target line = line 1~8
target_line <- c("01호선", "02호선", "03호선", "04호선",
                 "05호선", "06호선", "07호선", "08호선",
                 "1호선", "2호선", "3호선", "4호선",
                 "5호선", "6호선", "7호선", "8호선")

############################
#### Data preprocessing ####
############################

station.name$name_kor <- as.character(station.name$name_kor)
station.name <- station.name[station.name$name_kor != "이수",]
station.name$name_kor[which(station.name$name_kor == "4?19민주묘지")] <- "4.19민주묘지"
station.name <- station.name[station.name$line %in% target_line,]
station.name$line <- as.character(station.name$line)
station.name$line <- as.character(sapply(station.name$line,
                                         function(x) {strsplit(x, split="")[[1]][2]}))

station.loc$name_kor <- as.character(sapply(as.character(station.loc$name_kor), 
                                            function(x) {strsplit(x, split="\\(")[[1]][1]}))
station.loc <- station.loc[station.loc$line %in% target_line,]
station.loc <- as.data.frame(dplyr::select(station.loc, name_kor, lon, lat))

station.loc$name_kor[which(station.loc$name_kor=="이수")] <- "총신대입구"
# averaging location if there are several line passing the same station
station.loc <- as.data.frame(station.loc %>% group_by(name_kor) %>% 
                               summarise(lat = mean(lat), lon = mean(lon)))

# included in station.loc. lack of info 
station.loc2$name_kor <- as.character(station.loc2$name_kor)
station.loc2$name_kor[which(station.loc2$name_kor=="서울")] <- "서울역"

station.distance$name_kor <- as.character(station.distance$name_kor)
station.distance$name_kor[which(station.distance$name_kor=="이수")] <- "총신대입구"

# remove "(", ")"  in the station name
hourly.pplnum.2021$name_kor <- as.character(sapply(as.character(hourly.pplnum.2021$name_kor), 
                                                   function(x) {strsplit(x, split="\\(")[[1]][1]}))
hourly.pplnum.2021$name_kor[which(hourly.pplnum.2021$name_kor=="이수")] <- "총신대입구"
hourly.pplnum.2021 <- dplyr::select(hourly.pplnum.2021,-c("station_num"))
hourly.pplnum.2021 <- as.data.frame(hourly.pplnum.2021)

hourly.pplnum.2021_getin <- hourly.pplnum.2021[hourly.pplnum.2021$type=="승차",]
hourly.pplnum.2021_getout <- hourly.pplnum.2021[hourly.pplnum.2021$type=="하차",]

hourly.pplnum.2021_getin_line2 <- hourly.pplnum.2021_getin[hourly.pplnum.2021_getin$line=="2호선",]
hourly.pplnum.2021_getout_line2 <- hourly.pplnum.2021_getout[hourly.pplnum.2021_getout$line=="2호선",]

hourly.pplnum.2021_getin <- dplyr::select(hourly.pplnum.2021_getin,-c("line","type"))
hourly.pplnum.2021_getout <- dplyr::select(hourly.pplnum.2021_getout,-c("line","type"))
hourly.pplnum.2021_getin_line2 <- dplyr::select(hourly.pplnum.2021_getin_line2,-c("line","type"))
hourly.pplnum.2021_getout_line2 <- dplyr::select(hourly.pplnum.2021_getout_line2,-c("line","type"))


hourly.pplnum.2021_getin[,3:ncol(hourly.pplnum.2021_getin)] <- 
  sapply(hourly.pplnum.2021_getin[,3:ncol(hourly.pplnum.2021_getin)], as.numeric)

hourly.pplnum.2021_getout[,3:ncol(hourly.pplnum.2021_getout)] <- 
  sapply(hourly.pplnum.2021_getout[,3:ncol(hourly.pplnum.2021_getout)], as.numeric)

hourly.pplnum.2021_getin_line2[,3:ncol(hourly.pplnum.2021_getin_line2)] <- 
  sapply(hourly.pplnum.2021_getin_line2[,3:ncol(hourly.pplnum.2021_getin_line2)], as.numeric)

hourly.pplnum.2021_getout_line2[,3:ncol(hourly.pplnum.2021_getout_line2)] <- 
  sapply(hourly.pplnum.2021_getout_line2[,3:ncol(hourly.pplnum.2021_getout_line2)], as.numeric)

# summing the numbers of people across type, line
hourly.pplnum.2021_getin <- as.data.frame(hourly.pplnum.2021_getin %>% group_by(date, name_kor) %>% 
                                            summarise(across(everything(),sum)))

hourly.pplnum.2021_getout <- as.data.frame(hourly.pplnum.2021_getout %>% group_by(date, name_kor) %>% 
                                             summarise(across(everything(),sum)))

hourly.pplnum.2021_getin_line2 <- as.data.frame(hourly.pplnum.2021_getin_line2 %>% group_by(date, name_kor) %>% 
                                                  summarise(across(everything(),sum)))

hourly.pplnum.2021_getout_line2 <- as.data.frame(hourly.pplnum.2021_getout_line2 %>% group_by(date, name_kor) %>% 
                                                   summarise(across(everything(),sum)))


hourly.pplnum.2021_getin %>% group_by(date) %>% summarise(count = n())
hourly.pplnum.2021_getout %>% group_by(date) %>% summarise(count = n())
hourly.pplnum.2021_getin_line2 %>% group_by(date) %>% summarise(count = n())
hourly.pplnum.2021_getout_line2 %>% group_by(date) %>% summarise(count = n())

# stations with no info
station_removed <- c("신내", "강일", "하남검단산", "하남시청")

hourly.pplnum.2021_getin <- hourly.pplnum.2021_getin[!(hourly.pplnum.2021_getin$name_kor %in% station_removed), ]
hourly.pplnum.2021_getin %>% group_by(date) %>% summarise(count = n())

hourly.pplnum.2021_getout <- hourly.pplnum.2021_getout[!(hourly.pplnum.2021_getout$name_kor %in% station_removed), ]
hourly.pplnum.2021_getout %>% group_by(date) %>% summarise(count = n())

hourly.pplnum.2021_getin_line2 <- hourly.pplnum.2021_getin_line2[!(hourly.pplnum.2021_getin_line2$name_kor %in% station_removed), ]
hourly.pplnum.2021_getin_line2 %>% group_by(date) %>% summarise(count = n())

hourly.pplnum.2021_getout_line2 <- hourly.pplnum.2021_getout_line2[!(hourly.pplnum.2021_getout_line2$name_kor %in% station_removed), ]
hourly.pplnum.2021_getout_line2 %>% group_by(date) %>% summarise(count = n())


########################
#### aggregate info ####
########################

### Line 2 get in ###

# add location to getting on/off info
station.info.getin_line2 <- inner_join(hourly.pplnum.2021_getin_line2, station.loc, by='name_kor')

# add station's English name
station.info.getin_line2 <- inner_join(station.info.getin_line2, station.name[, c("name_kor", "name")][
  !duplicated(station.name[, c("name_kor", "name")]),], by='name_kor')

station.info.getin_line2$name <- as.character(station.info.getin_line2$name)

# column reordering
station.info.getin_line2 <- dplyr::select(station.info.getin_line2, 1,2, ncol(station.info.getin_line2):(ncol(station.info.getin_line2)-2),
                                          3:(ncol(station.info.getin_line2)-3))

# our target stations
target_station_line2 <- unique(station.info.getin_line2$name) # total 50 stations
target_station_kor_line2 <- unique(station.info.getin_line2$name_kor)

# fill NA values
tmp <- target_station_line2
tmp.kor <- target_station_kor_line2
tmp.loc <- sapply(station.info.getin_line2[!duplicated(station.info.getin_line2[,c("lon","lat")]),
                                           c("lon","lat")], as.numeric)

datenum <- length(unique(station.info.getin_line2$date))

tmp.mat <- as.data.frame(cbind(rep(tmp.kor, datenum), rep(tmp, datenum),
                               rep(tmp.loc[,1], datenum), rep(tmp.loc[,2], datenum)))

colnames(tmp.mat) <- c("name_kor", "name", "lon", "lat")
tmp.mat$name_kor <- as.character(tmp.mat$name_kor)
tmp.mat$name <- as.character(tmp.mat$name)
tmp.mat$lon <- as.numeric(as.character(tmp.mat$lon))
tmp.mat$lat <- as.numeric(as.character(tmp.mat$lat))
tmp.mat$date <- as.factor(as.character(rep(seq(as.Date("2021-01-01"), by = "day", length.out = datenum),
                                           each = length(target_station_line2))))

tmp2 <- left_join(tmp.mat, as.data.frame(station.info.getin_line2), 
                  by=c("date"="date","name_kor"="name_kor", 
                       "name"="name"))


station.info.getin_line2 <- dplyr::select(tmp2, 5,1,2,3,4,8:ncol(tmp2))
colnames(station.info.getin_line2)[4:5] <- c("lon", "lat")

station.info.getin_line2.timediv <- station.info.getin_line2
station.info.getin_line2.timediv$T1 <- station.info.getin_line2.timediv$t1
station.info.getin_line2.timediv$T2 <- rowSums(station.info.getin_line2.timediv[,7:10])
station.info.getin_line2.timediv$T3 <- rowSums(station.info.getin_line2.timediv[,11:16])
station.info.getin_line2.timediv$T4 <- rowSums(station.info.getin_line2.timediv[,17:20])
station.info.getin_line2.timediv$T5 <- rowSums(station.info.getin_line2.timediv[,21:24])

station.info.getin_line2.timediv <- station.info.getin_line2.timediv[, c(1,2,3,4,5,25:29)]



### Line 2 get out ###

# add location to getting on/off info
station.info.getout_line2 <- inner_join(hourly.pplnum.2021_getout_line2, station.loc, by='name_kor')

# add station's English name
station.info.getout_line2 <- inner_join(station.info.getout_line2, station.name[, c("name_kor", "name")][
  !duplicated(station.name[, c("name_kor", "name")]),], by='name_kor')

station.info.getout_line2$name <- as.character(station.info.getout_line2$name)

# column reordering
station.info.getout_line2 <- dplyr::select(station.info.getout_line2, 1,2, ncol(station.info.getout_line2):(ncol(station.info.getout_line2)-2),
                                           3:(ncol(station.info.getout_line2)-3))

# our target stations
target_station_line2 <- unique(station.info.getout_line2$name) # total 50 stations
target_station_kor_line2 <- unique(station.info.getout_line2$name_kor)

# fill NA values
tmp <- target_station_line2
tmp.kor <- target_station_kor_line2
tmp.loc <- sapply(station.info.getout_line2[!duplicated(station.info.getout_line2[,c("lon","lat")]),
                                            c("lon","lat")], as.numeric)

datenum <- length(unique(station.info.getout_line2$date))

tmp.mat <- as.data.frame(cbind(rep(tmp.kor, datenum), rep(tmp, datenum),
                               rep(tmp.loc[,1], datenum), rep(tmp.loc[,2], datenum)))

colnames(tmp.mat) <- c("name_kor", "name", "lon", "lat")
tmp.mat$name_kor <- as.character(tmp.mat$name_kor)
tmp.mat$name <- as.character(tmp.mat$name)
tmp.mat$lon <- as.numeric(as.character(tmp.mat$lon))
tmp.mat$lat <- as.numeric(as.character(tmp.mat$lat))
tmp.mat$date <- as.factor(as.character(rep(seq(as.Date("2021-01-01"), by = "day", length.out = datenum),
                                           each = length(target_station_line2))))

tmp2 <- left_join(tmp.mat, as.data.frame(station.info.getout_line2), 
                  by=c("date"="date","name_kor"="name_kor", 
                       "name"="name"))


station.info.getout_line2 <- dplyr::select(tmp2, 5,1,2,3,4,8:ncol(tmp2))
colnames(station.info.getout_line2)[4:5] <- c("lon", "lat")

station.info.getout_line2.timediv <- station.info.getout_line2
station.info.getout_line2.timediv$T1 <- station.info.getout_line2.timediv$t1
station.info.getout_line2.timediv$T2 <- rowSums(station.info.getout_line2.timediv[,7:10])
station.info.getout_line2.timediv$T3 <- rowSums(station.info.getout_line2.timediv[,11:16])
station.info.getout_line2.timediv$T4 <- rowSums(station.info.getout_line2.timediv[,17:20])
station.info.getout_line2.timediv$T5 <- rowSums(station.info.getout_line2.timediv[,21:24])

station.info.getout_line2.timediv <- station.info.getout_line2.timediv[, c(1,2,3,4,5,25:29)]


#########################################
#### edge weight matrix construction ####
#########################################
station.distance <- as.data.frame(station.distance) 
# remove "(", ")" in the name
station.distance$name_kor <- as.character(sapply(as.character(station.distance$name_kor), 
                                                 function(x) {strsplit(x, split="\\(")[[1]][1]}))
# remove blank at the end of name
station.distance$name_kor <- as.character(sapply(as.character(station.distance$name_kor), 
                                                 function(x) {strsplit(x, split=" ")[[1]][1]}))
station.distance$name_kor[which(station.distance$name_kor=="신내역")] <- "신내"
# add station's English name
station.distance_line2 <- station.distance[station.distance$name_kor %in% target_station_kor_line2,]

station.distance_line2 <- inner_join(station.distance_line2, 
                                     station.info.getin_line2[!duplicated(station.info.getin_line2[,c("name_kor","name")]),
                                                              c("name_kor", "name")], by='name_kor')

station.distance_line2[station.distance_line2$name_kor=="산성",]$btwdist <-
  station.distance_line2[station.distance_line2$name_kor=="산성",]$btwdist + 1.5 # 남위례역이 hourly data 정보에 없어서 빠지기 때문에 간격 더해줌


station.distance_line2[station.distance_line2$name_kor=="미사",]$btwdist <-
  station.distance_line2[station.distance_line2$name_kor=="미사",]$btwdist + 0.8 # 강일이 hourly data 정보에 없어서 빠지기 때문에 간격 더해줌

station.distance_line2$name <- as.character(station.distance_line2$name)

e.weight <- matrix(0,nrow=length(target_station_line2), ncol=length(target_station_line2))
colnames(e.weight) <- target_station_line2
rownames(e.weight) <- target_station_line2

for(i in 2:2){
  tmp <-station.distance_line2[station.distance_line2$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.weight["Seongsu", "Yongdap"] <- tmp$btwdist[45]
    e.weight["Yongdap", "Seongsu"] <- tmp$btwdist[45]
    for(j in (n+1):47){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
    e.weight["Sindorim", "Dorimcheon"] <- tmp$btwdist[49]
    e.weight["Dorimcheon", "Sindorim"] <- tmp$btwdist[49]
    for(j in 49:50){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.weight["Gangdong", "Dunchon-dong"] <- tmp$btwdist[47]
    e.weight["Dunchon-dong", "Gangdong"] <- tmp$btwdist[47]
    for(j in (n+1):52){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
    e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
  }
}

e.weight.old <- e.weight
e.weight.old[e.weight.old!=0] <- exp(-e.weight.old[e.weight.old!=0])
e.weight[e.weight!=0] <- exp(-(e.weight[e.weight!=0])^2 / mean(e.weight[e.weight!=0])^2)

e.weight_circle <- e.weight[(!colnames(e.weight) %in% c("Yongdap","Yongdu","Sinjeongnegeori","Sinseol-dong","Dorimcheon","Yangcheon-gu Office","Sindap")),
                            (!colnames(e.weight) %in% c("Yongdap","Yongdu","Sinjeongnegeori","Sinseol-dong","Dorimcheon","Yangcheon-gu Office","Sindap"))]

e.sp.weight <- NULL
e.color <- c() # for line color
color.cand <- c("blue", "yellowgreen", "orangered", "cyan",
                "darkorchid", "chocolate3", "darkolivegreen", "hotpink")
for(i in 2:2){
  tmp <-station.distance_line2[station.distance_line2$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.sp.weight <- rbind(e.sp.weight, c("Seongsu", "Yongdap", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Yongdap", "Seongsu", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):47){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
    e.sp.weight <- rbind(e.sp.weight, c("Sindorim", "Dorimcheon", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Dorimcheon", "Sindorim", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    
    for(j in 49:50){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.sp.weight<- rbind(e.sp.weight, c("Gangdong", "Dunchon-dong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight<- rbind(e.sp.weight, c("Dunchon-dong", "Gangdong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):52){
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
  }
}

tmp <- as.data.frame(e.sp.weight)
tmp$V1 <- as.character(tmp$V1)
tmp$V2 <- as.character(tmp$V2)
tmp2 <- exp(-as.numeric(as.character(tmp$V3)))
tmp$V3 <- exp(-(as.numeric(as.character(tmp$V3)))^2 / (mean(as.numeric(as.character(tmp$V3))))^2)

e.sp.weight <- tmp
e.sp.weight.old <- e.sp.weight
e.sp.weight.old[,3] <- tmp2


seoulmetro_line2 <- list()
# location
seoulmetro_line2$xy <- sapply(station.info.getin_line2[!duplicated(station.info.getin_line2[,c("lon","lat")]),
                                                       c("lon","lat")], as.numeric)
rownames(seoulmetro_line2$xy) <- target_station_line2

# weight matrix
seoulmetro_line2$A <- e.weight

# sparse weight matrix
seoulmetro_line2$sA <- e.sp.weight

# dist matrix
distmat <- e.weight.old
distmat[distmat!=0] <- -log(e.weight.old[e.weight.old!=0])

tmp <- e.sp.weight.old
tmp[,3] <- -log(tmp[,3])

seoulmetro_line2$dist <- distmat
seoulmetro_line2$sdist <- tmp

plot_graph(seoulmetro_line2)


L.metro_line2 <- laplacian_mat(seoulmetro_line2$A) # laplacian matrix
val1 <- eigensort(L.metro_line2)
evalues.metro_line2 <- val1$evalues
evectors.metro_line2 <- val1$evectors
# largest eigenvalue
lmax.metro_line2 <- max(evalues.metro_line2)

N.metro_line2 <- nrow(L.metro_line2)


## assign signals ##

# get in 
X.metro.getin.line2.timediv <- list()
# X.metro.getin.line2.timediv_wo_centering <- list()
date.metro_getin_line2_timediv <- as.character(unique(station.info.getin_line2.timediv$date))
 
# data.metro_getin_line2_timediv_wo_centering <- NULL  
for(j in 1:5){
  data.metro_getin_line2_timediv <- NULL
  for(i in date.metro_getin_line2_timediv){
    f <- station.info.getin_line2.timediv[station.info.getin_line2.timediv$date == i, (j+5)]
    data.metro_getin_line2_timediv <- cbind(data.metro_getin_line2_timediv, f)
  }
  
  data.metro_getin_line2_timediv <- rowMeans(data.metro_getin_line2_timediv)
  data.metro_getin_line2_timediv <- log(1+data.metro_getin_line2_timediv)
  data.metro_getin_line2_timediv <- data.metro_getin_line2_timediv - mean(data.metro_getin_line2_timediv)
  
  # X.metro.getin.line2.timediv_wo_centering[[j]] <- data.metro_getin_line2_timediv_wo_centering
  X.metro.getin.line2.timediv[[j]] <- data.metro_getin_line2_timediv
}


# get out
X.metro.getout.line2.timediv <- list()
date.metro_getout_line2_timediv <- as.character(unique(station.info.getout_line2.timediv$date))
  
for(j in 1:5){
  data.metro_getout_line2_timediv <- NULL
  for(i in date.metro_getout_line2_timediv){
    f <- station.info.getout_line2.timediv[station.info.getout_line2.timediv$date == i, (j+5)]
    data.metro_getout_line2_timediv <- cbind(data.metro_getout_line2_timediv, f)
  }
  
  data.metro_getout_line2_timediv <- rowMeans(data.metro_getout_line2_timediv)
  data.metro_getout_line2_timediv <- log(1+data.metro_getout_line2_timediv)
  data.metro_getout_line2_timediv <- data.metro_getout_line2_timediv - mean(data.metro_getout_line2_timediv)
  
  X.metro.getout.line2.timediv[[j]] <- data.metro_getout_line2_timediv
}


#############
## results ##
#############

# get in
res.metro.getin.line2.random_window <- GFPCA(X=X.metro.getin.line2.timediv, S=L.metro_line2, M=50, sigma=0.1, method="random")

M.metro.getin.line2 <- 50
# M is suitable?
(M.metro.getin.line2+1)*lmax.metro_line2 / M.metro.getin.line2^2 # sigma square
lmax.metro_line2^2 / M.metro.getin.line2^2

# scree plot?
plot(colSums(res.metro.getin.line2.random_window$tau.hat) / sum(res.metro.getin.line2.random_window$tau.hat), type="o")

res.metro.getin.line2.random_window <- GFPCA(X=X.metro.getin.line2.timediv, S=L.metro_line2, M=50, sigma=0.1, q=1, method="random")

# raw data
g1 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getin.line2.timediv[[1]], value="Passengers", ratio=0.6,
                         min=-3.47, max=2.22, mg=c(4,4,4,4), title=expression(italic(T)[1]), main.title.size = 20)
g2 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getin.line2.timediv[[2]], value="Passengers", ratio=0.6,
                         min=-3.47, max=2.22, mg=c(4,4,4,4), title=expression(italic(T)[2]), main.title.size = 20)
g3 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getin.line2.timediv[[3]], value="Passengers", ratio=0.6,
                         min=-3.47, max=2.22, mg=c(4,4,4,4), title=expression(italic(T)[3]), main.title.size = 20)
g4 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getin.line2.timediv[[4]], value="Passengers", ratio=0.6,
                         min=-3.47, max=2.22, mg=c(4,4,4,4), title=expression(italic(T)[4]), main.title.size = 20)
g5 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getin.line2.timediv[[5]], value="Passengers", ratio=0.6,
                         min=-3.47, max=2.22, mg=c(4,4,4,4), title=expression(italic(T)[5]), main.title.size = 20)
g6 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getin.line2.timediv[[5]], value="Passengers", ratio=0.6,
                         min=-3.47, max=2.22, mg=c(4,4,4,4), title=expression(italic(T)[5]), main.title.size = 20, signal=FALSE)

grid.arrange(g6,g1,g2,g3,g4,g5, nrow=2)

# eigenvectors
g7 <- plot_graph_custom3(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = evectors.metro_line2[,2], value="V", ratio=0.6,
                         min=-0.38, max=0.16, title=expression(v[2]))
g8 <- plot_graph_custom3(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = evectors.metro_line2[,4], value="V", ratio=0.6,
                         min=-0.38, max=0.16, title=expression(v[4]))
grid.arrange(g7,g8, nrow=1)



# get out
res.metro.getout.line2.random_window <- GFPCA(X=X.metro.getout.line2.timediv, S=L.metro_line2, M=50, sigma=0.1, method="random")

M.metro.getout.line2 <- 50
# M is suitable?
(M.metro.getout.line2+1)*lmax.metro_line2 / M.metro.getout.line2^2 # sigma square
lmax.metro_line2^2 / M.metro.getout.line2^2

# scree plot?
plot(colSums(res.metro.getout.line2.random_window$tau.hat) / sum(res.metro.getout.line2.random_window$tau.hat), type="o")

res.metro.getout.line2.random_window <- GFPCA(X=X.metro.getout.line2.timediv, S=L.metro_line2, M=50, sigma=0.1, q=1, method="random")

# raw data
g1 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getout.line2.timediv[[1]], value="Passengers", ratio=0.6,
                         min=-3.34, max=1.98, mg=c(4,4,4,4), title=expression(italic(T)[1]), main.title.size = 20)
g2 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getout.line2.timediv[[2]], value="Passengers", ratio=0.6,
                         min=-3.34, max=1.98, mg=c(4,4,4,4), title=expression(italic(T)[2]), main.title.size = 20)
g3 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getout.line2.timediv[[3]], value="Passengers", ratio=0.6,
                         min=-3.34, max=1.98, mg=c(4,4,4,4), title=expression(italic(T)[3]), main.title.size = 20)
g4 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getout.line2.timediv[[4]], value="Passengers", ratio=0.6,
                         min=-3.34, max=1.98, mg=c(4,4,4,4), title=expression(italic(T)[4]), main.title.size = 20)
g5 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getout.line2.timediv[[5]], value="Passengers", ratio=0.6,
                         min=-3.34, max=1.98, mg=c(4,4,4,4), title=expression(italic(T)[5]), main.title.size = 20)
g6 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = X.metro.getout.line2.timediv[[5]], value="Passengers", ratio=0.6,
                         min=-3.34, max=1.98, mg=c(4,4,4,4), title=expression(italic(T)[5]), main.title.size = 20, signal=FALSE)

grid.arrange(g6,g1,g2,g3,g4,g5, nrow=2)


par(mfrow=c(2,3), mar=c(5,5,4,2)+0.1, mgp=c(3,1,0))
plot(res.metro.getin.line2.random_window$tau.hat[,1], type="l", main="Graph spectral envelope",
     xlab="Graph frequency index", ylab=expression(tau[1]), cex.lab=2, cex.main=2, cex.axis=1.2)
abline(v=c(2,4), col="red", lty=2)

plot(-res.metro.getin.line2.random_window$H.hat[2,1,], type="l", xlab="Dimension index", ylab="Normalized amplitude", main=expression(lambda[2]^"SM"), cex.main=2, cex.lab=2, cex.axis=1.3)
points(1:5, -res.metro.getin.line2.random_window$H.hat[2,1,], pch=1, cex=1.7)
plot(res.metro.getin.line2.random_window$H.hat[4,1,], type="l", xlab="Dimension index", ylab="Normalized amplitude", main=expression(lambda[4]^"SM"), cex.main=2, cex.lab=2, cex.axis=1.3)
points(1:5, res.metro.getin.line2.random_window$H.hat[4,1,], pch=1, cex=1.7)

# par(mfrow=c(1,3), mar=c(5,5,4,2)+0.1, mgp=c(3,1,0))
plot(res.metro.getout.line2.random_window$tau.hat[,1], type="l", main="Graph spectral envelope",
     xlab="Graph frequency index", ylab=expression(tau[1]), cex.lab=2, cex.main=2, cex.axis=1.2)
abline(v=c(2,4), col="red", lty=2)

plot(res.metro.getout.line2.random_window$H.hat[2,1,], type="l", xlab="Dimension index", ylab="Normalized amplitude", main=expression(lambda[2]^"SM"), cex.main=2, cex.lab=2, cex.axis=1.3)
points(1:5, res.metro.getout.line2.random_window$H.hat[2,1,], pch=1, cex=1.7)
plot(res.metro.getout.line2.random_window$H.hat[4,1,], type="l", xlab="Dimension index", ylab="Normalized amplitude", main=expression(lambda[4]^"SM"), cex.main=2, cex.lab=2, cex.axis=1.3)
points(1:5, res.metro.getout.line2.random_window$H.hat[4,1,], pch=1, cex=1.7)
