import numpy as np
import pandas as pd
import itertools
from window_slider import Slider
import sys
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns

fasta=[]
test=[]
f1=open('500_100.contig.fa', 'r')
#Extracting sequences from the contig file
for line in f1:
        line=line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_seq=line[1:]
            if active_seq not in fasta:
                test.append(''.join(fasta))
                fasta=[]
            continue
        sequence=line
        fasta.append(sequence)
if fasta:
    test.append(''.join(fasta)) #Joining the sequences to make a list names test
test.pop(0)
f1.close()

#Making all possible combinations of tetranucleotides
nt=['A','G','C','T']
tnts_basic=np.array([''.join(i) for i in itertools.product(nt, repeat=4)])

#Empty dictonaries
dict1={}
dict3={}
dict2={}
n=0
#Opening a file which only consists of headers of the seq and putting them into a list
f2=open('try_contig.txt', 'r')
headers= f2.readlines()
#print(headers)
head=[]
badchar=['\n', ',']
for line in headers:
    line=''.join((filter(lambda i:i not in badchar, line)))
    head.append(line)
f2.close()

#Loop to run over each sequence present in test
for items in test:
    res=list(items)     #Parsing the sequence into a list
    res_arr=np.array(res)
    #Sliding window analysis
    size=4
    overlap=3
    slider=Slider(size, overlap)
    slider.fit(res_arr)
    tnts_gen=[]
    while True:
        data=slider.slide()
        tnts_data=''.join(data) #Making the tetranucleotides in the contig
        tnts_gen.append(tnts_data)
        if slider.reached_end_of_list(): break #Breaking the sliding window once the end of file is reached

    count=0
    count_arr=[]
    #Counting each element of list nt to list in genome
    for item in tnts_basic:
        count=0
        for item1 in tnts_gen:
            if item == item1:
                count=count+1
        count_arr.append(count)
    #Dictionary to hold tnts as keys and count as values
    for A,B in zip(tnts_basic, count_arr):
        dict1[A] = B

#Computing the reverse complement of each of the unique tetranucleotides
    complement={'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    tnts_rev_basic=[]
    rev_comp_arr=[]
    for i in range(128):
        rev_com="".join(complement.get(base, base) for base in reversed(tnts_basic[i]))
        tnts_rev_basic.append(tnts_basic[i])
        new_count=dict1[tnts_basic[i]]+dict1[rev_com]   #Adding values for each reverse_compliment
        rev_comp_arr.append(new_count)

#Nested dicitionary with header as the key, tnts as inner key and counts as calues
    dict2[head[n]]={}
    for key, value in zip(tnts_rev_basic, rev_comp_arr):
        dict2[head[n]][key]= value
    n=n+1

#Converting the dictionary to a pandas df
df=pd.DataFrame.from_dict(dict2)
#Transposing the dataframe
df=df.transpose()
#Writing the df as a csv file
df.to_csv("freq.csv")

#PCA for reducing the components to two dimensions
#fitting
scalar=StandardScaler()
scalar.fit(df)
scaled_data=scalar.transform(df)

pca=PCA(n_components=2)
pca.fit(scaled_data)
x_pca=pca.transform(scaled_data)

#Starting Kmeans to cluster the data into 3
mat_df=pd.DataFrame({'X':x_pca[:,0], 'Y':x_pca[:,1]}, index=zip(head))
print(mat_df)
    
#Initialize the centroids
c1=(-1, 2)
c2=(-5, -3)
c3=(2, 2)

#Calculating initial distance
def calculate_distance(centroid, X, Y):
    dist=[]
    c_x, c_y=centroid

    for x,y in list(zip(X,Y)):
        root_diff_x=(x - c_x) ** 2
        root_diff_y=(y - c_y) ** 2
        distance=np.sqrt(root_diff_x + root_diff_y)
        dist.append(distance)
    return dist

mat_df['C1_dist']=calculate_distance(c1, mat_df.X, mat_df.Y)
mat_df['C2_dist']=calculate_distance(c2, mat_df.X, mat_df.Y)
mat_df['C3_dist']=calculate_distance(c3, mat_df.X, mat_df.Y)

#Updating the column names in the dataframe
mat_df['Cluster']=mat_df[['C1_dist', 'C2_dist', 'C3_dist']].idxmin(axis=1)
mat_df['Cluster']=mat_df['Cluster'].map({'C1_dist': 'C1', 'C2_dist': 'C2', 'C3_dist': 'C3'})

#Finding new centroid points from the computed data
x_new_centroid1= mat_df[mat_df['Cluster']=='C1']['X'].mean()
y_new_centroid1= mat_df[mat_df['Cluster']=='C1']['Y'].mean()
centroid1=(x_new_centroid1,y_new_centroid1)


x_new_centroid2= mat_df[mat_df['Cluster']=='C2']['X'].mean()
y_new_centroid2= mat_df[mat_df['Cluster']=='C2']['Y'].mean()
centroid2=(x_new_centroid2,y_new_centroid2)

x_new_centroid3= mat_df[mat_df['Cluster']=='C3']['X'].mean()
y_new_centroid3= mat_df[mat_df['Cluster']=='C3']['Y'].mean()
centroid3=(x_new_centroid3,y_new_centroid3)

#Iteration to find the optimum clustering
while True:
    Clust=[]
    Clust=mat_df.to_numpy()
    mat_df['C1_dist']=calculate_distance(centroid1, mat_df.X, mat_df.Y)
    mat_df['C2_dist']=calculate_distance(centroid2, mat_df.X, mat_df.Y)
    mat_df['C3_dist']=calculate_distance(centroid3, mat_df.X, mat_df.Y)

    mat_df['Cluster']=mat_df[['C1_dist', 'C2_dist', 'C3_dist']].idxmin(axis=1)

    mat_df['Cluster']=mat_df['Cluster'].map({'C1_dist': 'C1', 'C2_dist': 'C2', 'C3_dist': 'C3'})

    Clust2=[]
    Clust2=mat_df.to_numpy()
    #Comparing two clustering lists to continue iteration or break the loop
    x=(Clust[:,5]==Clust2[:,5]).all()
    if(x==True):
        print("This is the optimum cluster")
        print(mat_df)
        mat_df.to_csv("clustered.csv")
        print("The frequency data is saved as 'freq.csv' and the clustered data is saved as 'clustered.csv'")
        sns.lmplot(data=mat_df, x='X', y='Y', hue='Cluster', fit_reg=False, legend=True, legend_out=True)
        plt.show()
        break
    else:
        x_new_centroid1= mat_df[mat_df['Cluster']=='C1']['X'].mean()
        y_new_centroid1= mat_df[mat_df['Cluster']=='C1']['Y'].mean()
        centroid1=(x_new_centroid1,y_new_centroid1)

        x_new_centroid2= mat_df[mat_df['Cluster']=='C2']['X'].mean()
        y_new_centroid2= mat_df[mat_df['Cluster']=='C2']['Y'].mean()
        centroid2=(x_new_centroid2,y_new_centroid2)

        x_new_centroid3= mat_df[mat_df['Cluster']=='C3']['X'].mean()
        y_new_centroid3= mat_df[mat_df['Cluster']=='C3']['Y'].mean()
        centroid3=(x_new_centroid3,y_new_centroid3)

    



