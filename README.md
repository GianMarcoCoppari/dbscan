# DBSCAN Algorithm - C++ Implementation

The program implements a DBSCAN algorithm from scratch in C++ language. It laso provides a fucntion for reading comma separated files to get dataset from external text files and same randomly initialized toy dataset to see the algorithm in action.

## Basic Concept
DBSCAN (*Density Based Spatial Clustering for Application with Noise*) is a clustering algorithm and an exmaple of unsupervised learning. As the name suggests the clustering is made accroding to density of points in a region of feature space: points near to each are likely to belong to the same category, whereas farer points are more likely to belong to different groups.

## The Algorithm
The algorithm labels a point $P$ with a label $l$, if there are at least $n_{\text{min}}$ neighbors in the $\epsilon$-vicinity of the point. Such a point $P$ is said to be density-connected to the others if the $\epsilon$-vicinity and in particular it will be called *core point*, in contrast to those that are density connected to other point, but do not have enough neighbours to forma group: these are *border points*. Point that are not density-connected to other points are said *isolated* or *noise points*.

Here, $\epsilon$ and $n_{\text{min}}$ are the hyperparameters of the model, i.e. those parameters that have to be provided from the user to make the algorithm work. Note that the result of the clustering will depend on the values of those hyperparameters.

The algorithm, then, goes through all the points in the dataset and classify them according to the number od points in their respective $\epsilon$-vicinity. Vicinity of points is determined by a *metric function* $d: \mathbb{R}^{n} \times \mathbb{R}^{n} \to \mathbb{R}^{+}$, which is taken to be the usual Euclidean distance in $\mathbb{R}^{n}$. Note that vicinity in $\mathbb{R}^{n}$ means similarity of feature and, moreover, a point $Q$ is in the $\epsilon$-vicinity of another point $P$ if and only if $d(Q, P) < \epsilon$.
If the $\epsilon$-vicinity of point $P$ contains at least $n_{\text{min}}$ points, then a label $l_{1}$ is attacched to point $P$. If this does not happen, i.e. $P$ has not enough points aound, it will be classified as *noise*. Note that in this case the label might change in case $P$ is in the $\epsilon$-vicinity of another point $P'$ that has a sufficient number of neighbors (i.e. more than $n_{\text{min}} - 1$ neighbors).

At the end of the procedure all points of the dataset $X$ are labelled either with a proper cluster label or as noise.

## Explanation of Code
The code I present is composed of different functions, implementing the DBSCAN Algorithm itself and some auxiliary function to make the algorithm work. The result is a set of text files containing the entry of the dataset divided by membership label.

### Pseudocode for DBSCAN Algorithm
The signature for the DBSCAN function is
`DBSCAN(dataset, epsilon, minsize)`
as it takes a dataset `dataset` and the two hyperparameters `epsilon` and `minsize`. Let `noise` and `undefined` be the labels for noise points and for point that have not been classified yet. Then a pseudocode for the algorithm might be the following

````
DBSCAN(dataset, epsilon, minsize)
    labels = EmptySet
    C = 0
    FOR entry IN X:
        IF (label != undefined) continue;
        
        neighbours = Subset(entry, dataset, query)
        N = neighbours.size()
        IF (N < minsize) label = noise
        
        label = ++C;
            
        WHILE neighbors NOT EMPTY:
           ComputeLabel(neighbors.back(), dataset, query)
            nearneaighbors = Subset(neiugbors.back(), dataset, query)
            IF (nearneighbours.size() >= minsize)
                ExpandCluster(X, neighbors.back(), neighbors, C, query, minsize)
            REMOVE neighbors.back()
                
return labels
````
In the pseudocode above `Subset` and `ExpandCluster` are two auxiliary functions that compute the $\epsilon$-vicinity of a point and the union set of two groups of neighbors, respectively. `ComputeLabel` attach the correct label to the neighbor point.

### Auxiliary Functions
The pseudocode for the two auxiliary function might be the following

````
Subset(P, X, query)
    s = EmptySet;
    FOR point IN X:
        IF (query(P, point) = TRUE):
            s.push_back(point)
            
    return s
````

````
ExpandCluster(X, Last, Neighbors, C, query, minsize)
    NeighborsToNeighbor = Subset(Last, X, query)
    
    IF (NeighorsToNeighbour.size() >= minsize)
        UnionSet(Neighbors, NeighborsToNeighbors)
````


## Computational Complexity & Improvements
The main source of complexity in the code is the function `Subset` that scans all the points of the dataset at least ones every time it is invoked. If $n$ is the number of entries in the dataset the complexity of this function is $O(n^2)$ as there are $n(n-1)$ check to compute. A great improvement might come from the implementation of a tree structure to store data rather than a plain database, such as a `Quadtree` algorithm if the dataset contains bidimensional data.

## Ouput Example
Using LaTeX we can plot reults. Here a picture of circular crowns clustering using the toy dataset produced in the code. The DBSCAN algorithm has been prepared with the folowing values of hyperparameters: $\epsilon = 0.5$ and $n_{\text{min}} = 8$. Here is the link to the image [dbscan](https://github.com/GianMarcoCoppari/dbscan/blob/master/dbscan.pdf).
The program has been compiled with the line
`g++ -std=c++2a -Wall -Wextra -Wpedantic -O3 main.cpp`
and the executable file can be execute using the command `.\a.exe`. 
