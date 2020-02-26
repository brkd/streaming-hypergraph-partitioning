# streaming-hypergraph-partitioning

A hypergraph is a pair of two sets: A set of vertices and a set of hyperedges, each of which can connect any number of vertices. 
Hypergraphs can be employed for various purposes, including the modeling of sparse datasets that describe properties of social networks.

One can imagine the hypernets of such a graph simply as the subsets of the set of vertices, visualized here with a graph of 7 vertices and 4 nets:
![hg](https://user-images.githubusercontent.com/39488383/75326044-f6594c80-588a-11ea-9cb6-9193b101317a.png)

This project aims to develop a sufficiently fast and memory-efficient method to partition a hypergraph into roughly equal parts in such a way that the number of vertices in different partitions connected by the same hyperedge is minimized, and it is also paramount that this method can run with a streaming model. 
It is understandably harder to work on a dataset when it is streaming in real-time instead of being already recorded down, since 
the data points are seen only once and memory is far more restricted.

With a way to partition large datasets in hypergraph form, one can distribute the workload onto different devices, 
possibly eliminate the needless duplication of data and improve both the yield and the pace of their work. 
The practice of effectively partitioning large graphs is crucial in real-world applications like minimizing communication cost or I/O cost and maximization of robustness in social graphs.

While partitioning a graph, the cost function we want to minimize concerns the number of "cuts". Generally, a cut in a graph causes the vertices of said graph to be divided into disjoint subsets. Similarly, if two vertices connected by an edge are placed into different partitions, thus becoming elements of disjoint sets, this edge is said to be "cut", or rather it is called a "cut-net".
In an attempt to reduce the abstractness, one can imagine the vertices to be any kind of datapoints and the edges to be a unique trait that binds tham together. Understandably, dividing such datapoints into separate partitions hurts the effectiveness of whatever task was to be carried out.

In the image given above, putting v5 and v6 into different partitions cuts e3, while placing v2 and v3 into different partitions cuts both e1 and e2. v4 and v7, as they are currently, do not pose any threats performance-wise and can be placed into whichever partition.

While partitioning is in process, multiple metrics are used for calculating the "score" of a partition, in other words the suitability of a partition for the vertex to be placed. Call the current vertex to be placed V, and the partition being evaluated P. One can keep the set of nets that already have a pin in P. The number of such nets that V is also a part of is the "connectivity" of V in P; and the bigger it is, the more suitable P is for placing V.

In addition, the sizes of the partitions are optimally close to equal. The ratio of the number of vertices to the number of partitions is the "capacity constraint". It is, in some sense, the "expected value" regarding the sizes of the partitions. For instance, dividing 10,000,000 vertices into 2048 partitions, one could expect 10000000/2048 ≈ 4883 vertices in each partition. While this is very rarely the case, the actual numbers in practice are still rather close to this number.

The "load penalty" for a partition P is the reduction to its appeal caused by its size. The ratio of the size of P to the capacity constraint is a number that clarifies the number of vertices one can still place in P without overgrowing it.    

For a vertex V, a partition P far from its capacity constraint that it has high connectivity in is ideal choice.

Other things that affect the partitioning process are values like "imbalance" or a "slack value". While a capacity constraint and a load penalty can be used for the sizes of the partitions, if no partition at the time of placement stand out, some partitions can still overgrow. An imbalance value is added to monitor the sizes of the partitions, and if need be, deny the vertex's placement to a partition completely. It defines how much imbalance the partitions can have between each other, that is to say the ratio of the sizes of the bigger partitions to the smaller partitions is capped to never go over a certain "imbalance" value.  

The variations of partitioning algorithms that are being worked on are:

*Random Partitioning

*Part-to-Net

*Net-to-Part

*Net-to-Part_i

*Bloom Filter

*Multilayered Bloom Filter

Part-to-Net and Net-to-Part approaches include keeping information on which nets have a leg
in a partition or which partitions the queried net has vertices in. That information is used
for determining which partition would raise the yield of a certain cost function the least.

Bloom filtered approaches can save time on vertex-partition association checks. There are 
a few different implementations of said algorithms present, where we try to save much more
time by grouping bloom filters into multiple layers or keeping edge-vertex partitions
together in the filter.
