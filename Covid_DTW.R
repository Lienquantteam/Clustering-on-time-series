library(dtwclust)
library(TSstudio)
library(readxl)
library(TSclust)
library(MASS) 
library(tidyverse)
library(TSstudio)
library(factoextra)
rm(list = ls())


df_covid = read.csv('D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/df_covid_pivot1.csv')
View(df_covid)
dim(df_covid) #115 188

df_covid$Date = as.POSIXct(df_covid$Date)
ts_plot(df_covid)

df_canada_netherlands = df_covid %>% select(Date,Canada, Netherlands)

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

write.csv()

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

###############################################################################
# find optimal k
#Silhouette method

fviz_nbclust(as.matrix(df_canada_netherlands[-1]), FUN = hcut, method = "silhouette") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Silhouette method HC")

#elbow method
fviz_nbclust(df_canada_netherlands[-1], FUN = hcut, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method HC")   

# comment: k = 2: optimal

###############################################################################
# transpose date

df_canada_netherlands_T = t(df_canada_netherlands)

#df_canada_netherlands_T  <- pivot_longer(df_canada_netherlands , cols = 2:3, names_to = "countries", values_to = "Covid_cases")

colnames(df_canada_netherlands_T) = NULL
df_canada_netherlands_T = df_canada_netherlands_T[-1]

write.csv(df_canada_netherlands_T,'D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/df_covid_pivot2.csv')

View(df_canada_netherlands_T)
###############################################################################

df_T = read.csv('D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/df_covid_pivot_T.csv')
str(df_T)
View(df_T)
rownames(df_T)= df_T[,1]
df_T = df_T[-1]

clust.pam <- tsclust(df_T, type="partitional", k=2L, distance="dtw", clustering="pam")

#clust.pam
#partitional clustering with 2 clusters
#Using dtw distance
#Using pam centroids

#Time required for analysis:
#  user  system elapsed 
#48.94    2.82   56.38 

#Cluster sizes with average intra-cluster distance:
  
#  size av_dist
#1  166  222442
#2   21 6306667

plot(clust.pam, type = "sc")
plot(clust.pam, type = "sc", clus = 1L)
plot(clust.pam, type = "series", clus = 1L)

result = t(cbind(labels, cluster = clust.pam@cluster))
result = as.data.frame(result)
View(result)

write.csv(result,'D:/Lien - Langara/Year 2/Capstone project/Lien/Lien folder/DTW/Covid 19_DTW/result.csv')

library(clValid)
clmethods = c("hierarchical","kmeans","pam")

df_T_scale = scale(df_T)
internal = clValid(df_T_scale, nClust = 2:5, clMethods = clmethods, validation = "internal")

summary(internal)


######################################################################################

#Hierarchical method

clust.hier2 <- tsclust(df_T, type = "h", k = 2L, distance = "dtw")
cutree(clust.hier2, k=2L)

#hierarchical clustering with 2 clusters
#Using dtw distance
#Using PAM (Hierarchical) centroids
#Using method average 

#Time required for analysis:
#  user  system elapsed 
#5.44    0.05    5.55 

#Cluster sizes with average intra-cluster distance:
#  size   av_dist
#1   35  8777.343
#2   80 36072.838

plot(clust.hier, type = "dendrogram")
plot(clust.hier, type = "series")
plot(clust.hier, type = "centroids")

##########################################################################################

clust.pam6 <- tsclust(df_T, type="partitional", k=6L, distance="dtw", clustering="pam")
plot(clust.pam6, type = "sc")


clust.hier6 <- tsclust(df_covid[-1], type = "h", k = 6L, distance = "dtw")
plot(clust.hier6, type = "dendrogram")
