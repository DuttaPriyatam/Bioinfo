alpha=0.7	#alpha value
data=4		
cluster=2	#no of clusters

#Input vectors
x1=c(1,1,0,0)	
x2=c(0,1,0,0)
x3=c(0,0,1,0)
x4=c(0,0,1,1)
col.names=c("x1","x2","x3","x4")
row.names=c()

#Making an array of the input vector
X_vec=array(c(x1,x2,x3,x4), dim=c(4,4), dimnames=list(row.names,col.names))
print(X_vec)

#Weighted vectors
C1=c(0.16, 1.00, 0.37, 0.63)
C2=c(0.09, 0.02, 0.48, 0.03)
col.names=c("C1", "C2")

#Making array of weigthed vectors
W_vec=array(c(C1,C2), dim=c(4,2), dimnames=list(row.names,col.names))
print(W_vec)
d1=0
d2=0

#1000 iterations
for(iter in 1:1000)
{
	print(paste("Iteration ", iter))
	for(i in 1:data)
	{
		#Finding distance for each cluster
		for(n in 1:data)
		{
			d1=d1+(C1[n]-X_vec[n,i])^2
			d2=d2+(C2[n]-X_vec[n,i])^2
		}
		print(paste("Distance of C1: ",d1))
		print(paste("Distance of C2: ",d2))
		if(d1<d2)
		{
			C1=c()
		}
		if(d1>d2)
		{
			C2=c()
		}
		x=i
		for(j in 1:data)
		{
			#Updating weight vector
			if(d1<d2)
			{
				W_vec[j,1]=W_vec[j,1]*(1-alpha)+alpha*X_vec[j,x]
				C1=append(C1,W_vec[j,1])
			}
			if(d1>d2)
			{	
				W_vec[j,2]=W_vec[j,2]*(1-alpha)+alpha*X_vec[j,x]
				C2=append(C2,W_vec[j,2])
			}
		}	
		print("Updated weight vector:-" )
		print(W_vec)
		d1=0
		d2=0
	}
}
