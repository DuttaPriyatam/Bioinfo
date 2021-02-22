#R_version used:- R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
#k-mean Clustering
#For a 2D dataset, with 10 datapoints and arranged into 2 clusters


k=2       #Taking two clusters into consideration
data_points=10  #No of data points
data_num1=0     #Number of Data belonging to each cluster
data_num2=0
lower_point=1   #To obtain initial centroid values
higher_point=6

col.names=c("X", "Y")
row.names=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10") 
x= c(1.0,3.5,2.0,2.0,1.0,5.0,5.0,2.0,5.0,6.0)
y= c(1.0,3.5,3.0,2.5,3.0,5.5,6.5,3.5,5.0,6.5)

#Empty vectors to contain the information about clusters, and which cluster the data points belong to.
clust= c()   
x_cluster1=c()
y_cluster1=c()
x_cluster2=c()
y_cluster2=c()

sample=array(c(x,y), dim=c(10,2), dimnames=list(row.names,col.names))   #Making a 2D array
print(sample)   #Printing the array

#Obtaining the initial centroid values
centroid_1=array(c(x[lower_point],y[lower_point]))  
centroid_2=array(c(x[higher_point],y[higher_point]))

#Printing the initial centroid values
print("Initial centroid for cluster 1: ")
print(centroid_1)
print("Initial centroid for cluster 2: ")  
print(centroid_2)


Distance1=c()
Distance2=c()
#Initial clustering based on points which points lie closer to the centroids and taking these datapoints in separate vectors
for(i in 1:data_points)
{
    Dist_1=sqrt((((sample[i,1]-centroid_1[1])^2)+(sample[i,2]-centroid_1[2])^2))
    Dist_2=sqrt((((sample[i,1]-centroid_2[1])^2)+(sample[i,2]-centroid_2[2])^2))
    Distance1=append(Distance1, Dist_1)
    Distance2=append(Distance2, Dist_2)

    if(Dist_1<Dist_2)
    {
        clust=append(clust, "C1")
        x_cluster1=append(x_cluster1,sample[i,1])
        y_cluster1=append(y_cluster1,sample[i,2])
        data_num1=data_num1+1
    }
    else if(Dist_1>Dist_2)
    {  
        clust=append(clust, "C2")
        x_cluster2=append(x_cluster2,sample[i,1])
        y_cluster2=append(y_cluster2,sample[i,2])
        data_num2=data_num2+1
    }
}

#Printing the initial cluster
col.names_2=c("X", "Y", "Dist to cenroid(Cluster 1)", "Dist to centroid(Cluster 2)", "Cluster")
cluster=array(c(x,y,Distance1,Distance2, clust),  dim=c(10,5), dimnames=list(row.names,col.names_2))
print(cluster)

#This block will iterate the steps until an optimum cluster is formed
for(i in 1:10)
{
    #Finding the mean centroid value for both the clusters 
    for(i in 1:length(x_cluster1))
    {
        x_new1=(x_cluster1[i]+centroid_1[1])/2
        y_new1=(y_cluster1[i]+centroid_1[2])/2
        centroid_1=array(c(x_new1,y_new1)) 
    }
    for(i in 1:length(x_cluster2))    
    {
        x_new2=(x_cluster2[i]+centroid_2[1])/2
        y_new2=(y_cluster2[i]+centroid_2[2])/2
        centroid_2=array(c(x_new2,y_new2))
    }

    #New centroids after finding the mean centroid
    print("Calculated centroid for cluster 1")
    print(centroid_1)
    print("Calculated centroid for cluster 2")
    print(centroid_2)

    #Taking empty vectors to calculate the euclidean distance
    Euclidean_distance_clust1=c()
    Euclidean_distance_clust2=c()
    #Finding the Euclidean distance between the centroid and data points
    for(i in 1:data_points)
    {
        Euclid_dist_1=sqrt((((x[i]-centroid_1[1])^2)+(y[i]-centroid_1[2])^2))
        Euclid_dist_2=sqrt((((x[i]-centroid_2[1])^2)+(y[i]-centroid_2[2])^2))
        Euclidean_distance_clust1=append(Euclidean_distance_clust1, Euclid_dist_1)
        Euclidean_distance_clust2=append(Euclidean_distance_clust2, Euclid_dist_2)
    }

    #Computing which cluster it should belong to
    x_cluster1=c()
    y_cluster1=c()
    x_cluster2=c()
    y_cluster2=c()
    clust_new=c()
    for(i in 1:data_points)
    {
        if(Euclidean_distance_clust1[i]<Euclidean_distance_clust2[i])
        {
                clust_new=append(clust_new, "C1")
                x_cluster1=append(x_cluster1,x[i])
                y_cluster1=append(y_cluster1,y[i])            
        }
        else if(Euclidean_distance_clust1[i]>Euclidean_distance_clust2[i])
        {  
                clust_new=append(clust_new, "C2")
                x_cluster2=append(x_cluster2,x[i])
                y_cluster2=append(y_cluster2,y[i]) 
        }
    }

    #Making the comprehensive table to list the Euclidean distance and the cluster which the datapoint belongs to
    col.names_2=c("X", "Y", "Dist to mean(Cluster 1)", "Dist to mean(Cluster 2)", "Cluster")
    tab_cluster=array(c(x, y, Euclidean_distance_clust1, Euclidean_distance_clust2, clust_new),  dim=c(10,5), dimnames=list(row.names,col.names_2))
    print(tab_cluster)


    #Comparing cluster initially found with the one generated by the centroid
    count=0
    for(r in 1:nrow(tab_cluster))
    {
            if(clust[r]!=tab_cluster[r,5])
            {
                count=count+1
            }
    }
    #Taking the string for the new cluster assignment according to data points
    clust=c()
    for(r in 1:nrow(tab_cluster))
    {
        clust=append(clust, tab_cluster[r,5])
    }
    x=c()
    y=c()
    x=c(x_cluster1, x_cluster2) #Making new datapoints corresponding to the cluster they belong to
    y=c(y_cluster1, y_cluster2)
    

    #Condition to break the loop when the optimum cluster is formed
    if(count==0)
    {
        print("This is the optimum cluster")
        break
    }
}


