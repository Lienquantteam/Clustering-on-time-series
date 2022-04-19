rm(list = ls())

library(dtwclust)
library(TSstudio)
library(readxl)
library(TSclust)
library(MASS) 
library(tidyverse)
library(TSstudio)
library(factoextra)
###############################################################
# USE DAILY CASES INSTEAD
###############################################################

df_covid_new = read.csv('D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/full_grouped/full_grouped - new.csv')
View(df_covid_new)

sum(is.na(df_covid_new))

library(tidyverse)
library(dplyr) # for transpose purpose
library(tidyr)

df_covid_new.T = df_covid_new %>% spread(Country.Region, New.cases)
View(df_covid_new.T)

#write.csv(df_covid_new.T,'D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/df_covid_new.T.csv')

df_covid_new.T$Date = as.POSIXct(df_covid_new.T$Date)

ts_plot(df_covid_new.T)

df_covid_new.T1 = df_covid_new %>% spread(Date, New.cases)
View(df_covid_new.T1)

rownames(df_covid_new.T1 )= df_covid_new.T1 [,1]
df_covid_new.T1 =df_covid_new.T1[-1]

write.csv(df_covid_new.T1,'D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/df_covid_new.T1.csv')

# check the stability of clustering

library(clValid)
clmethods = c("hierarchical","kmeans","pam")

internal = clValid(df_covid_new.T1, nClust = 2:5, clMethods = clmethods, validation = "internal")
summary(internal)
#Clustering Methods:
#  hierarchical kmeans pam 

#Cluster sizes:
#  2 3 4 5 

#Validation Measures:
#  2       3       4       5

#hierarchical Connectivity   2.9290  5.8579  9.7159 12.6448
#Dunn           1.0743  1.4232  0.6891  0.7203
#Silhouette     0.9541  0.9396  0.8821  0.8471
#kmeans       Connectivity   4.8579  5.8579 15.3139 26.3238
#Dunn           0.5902  1.4232  0.1396  0.1227
#Silhouette     0.9512  0.9396  0.8464  0.8198
#pam          Connectivity   2.9290 15.0246 17.9536 20.9988
#Dunn           1.0743  0.0534  0.1359  0.1359
#Silhouette     0.9541  0.8326  0.8179  0.8335

#Optimal Scores:

#             Score  Method       Clusters
#Connectivity 2.9290 hierarchical 2       
#Dunn         1.4232 hierarchical 3       
#Silhouette   0.9541 hierarchical 2       

#Hierarchical method, k=2

hc_2 = tsclust(df_covid_new.T1, type = "h", k = 2L, distance = "dtw")
cut_2 = cutree(hc_2, k=2L) 
#US, Brazil : group 2

plot(hc_2, type = "dendrogram")

#Hierarchical method, k=3
hc_3 = tsclust(df_covid_new.T1, type = "h", k = 3L, distance = "dtw")
cut_3 = cutree(hc_3, k=3L) 

#US, Brazil : group 2
#India, Russia: group 3

plot(hc_3, type = "dendrogram")

################################################################################
# ILLUSTRATE HOW TO MAKE DTW MATRIX FOR 2 TIME SERIES, USING COVID DATASET
#################################################################################
df_US_China = df_covid_new.T %>% filter(Country.Region %in% c('US', 'China'))

df_canada_netherlands = df_covid_new.T %>% select(Date,Canada, Netherlands)

df_canada_netherlands$Date =  as.POSIXct(df_canada_netherlands$Date)
ts_plot(df_canada_netherlands)
#The matching calculated above can also be portrayed on a graph. 
#Diagonal lines portray a situation where data has been matched on-to-one. 
#Vertical and horizontal lines mean that one point on one axis represents more than one point on another.
plot(dtw(df_canada_netherlands$Canada, df_canada_netherlands$Netherlands), xlab="Canada", ylab="Netherlands", 
     xaxp  = c(0,10,10), yaxp = c(0,10,10))

#The same plot can be made using type="threeway" that additionally shows the plots on x ad y axis. 
#The dtw functions has to also have keep=TRUE defined.
plot(dtw(df_canada_netherlands$Canada, df_canada_netherlands$Netherlands, keep=TRUE), 
     xlab="Canada", ylab="Netherlands", xaxp  = c(0,10,10), yaxp = c(0,10,10), type="threeway")

#Dtw step patterns
#One can also plot the connections as below using type="twoway"

plot(dtw(df_canada_netherlands$Canada, df_canada_netherlands$Netherlands, keep=TRUE), 
     xaxp  = c(0,10,10), yaxp = c(0,10,10), type="twoway", col=c('blue', 'magenta'))

#DTW by default uses symmetric2. In tabs there are calculations performed for symmetric1 and symmetric2 step patterns.

#Because an array has it's point [1,1] in the upper left corner all matrices are rotated 90 degrees to 
#the right compared to the results from dtw. Traditionally the 'beginning' of the 
#dtw matrix should be in the lower left corner

dtw(df_canada_netherlands$Canada, df_canada_netherlands$Netherlands)$stepPattern

dtw(df_canada_netherlands$Canada, df_canada_netherlands$Netherlands, step.pattern = symmetric1)$stepPattern

#Alignment plot
#Creating a zeroes matrix

dtw_matrix = matrix(rep(c(0),13225), nrow=115, ncol=115, byrow = TRUE)
dtw_matrix

dim(dtw_matrix)

#The matrix is calculated according to the step pattern formula. 
#Additionally d(0,j)=d(i,0)=d(0,0)=???. That is why it is not included in calculations - 
#will never be accounted for by min() function. I will be using the default euclidean distance between points

#First element:

a1 = df_canada_netherlands$Canada
length(a1)
a2 = df_canada_netherlands$Netherlands

dtw_matrix[1,1] = sqrt(a1[1]^2 + a2[1]^2)
dtw_matrix

for (i in 2:115){
  dtw_matrix[i,1] = sqrt((a1[i] - a2[1])^2) + dtw_matrix[i-1,1]
}
dtw_matrix

#First column:

for (i in 2:115){
  dtw_matrix[i,1] = sqrt((a1[i] - a2[1])^2) + dtw_matrix[i-1,1]
}
dtw_matrix

#First row:

for (j in 2:115){
  dtw_matrix[1,j] = sqrt((a1[1] - a2[j])^2) + dtw_matrix[1,j-1]
}
dtw_matrix


#The rest of the matrix:

for (i in 2:115){
  for (j in 2:115){
    dtw_matrix[i,j] = sqrt((a1[i] - a2[j])^2) + min(dtw_matrix[i,j-1], dtw_matrix[i-1,j], dtw_matrix[i-1,j-1] + sqrt((a1[i] - a2[j])^2))
  }
}
dtw_matrix

View(dtw_matrix)

#Find the optimal alignment:

path = c(115,115) # starting with furthest place in matrix
i = 115
j = 115
while(i>1 & j>1){
  if (j == 1) {
    j = j - 1
  } else if (i == 1) {
    i = i - 1
  } else if (dtw_matrix[i,j-1] == min(dtw_matrix[i-1, j-1], dtw_matrix[i-1, j], dtw_matrix[i, j-1])){
    j = j - 1
  } else if (dtw_matrix[i-1,j-1] == min(dtw_matrix[i-1, j-1], dtw_matrix[i-1, j], dtw_matrix[i, j-1])){
    i = i - 1
    j = j - 1
  } else {
    i = i - 1
  }
  path = rbind(path, c(i,j))
}
path = rbind(path, c(1,1))

plot(dtw(a1,a2))
points(path[,1], path[,2], type="l")

#write.csv(dtw_matrix,'D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/dtw_matrix_new.csv')
