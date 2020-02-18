# streaming-hypergraph-partitioning

A hypergraph is a pair of two sets: A set of vertices and a set of hyperedges, each of which can connect any number of vertices. 
Hypergraphs can be employed for various purposes, including the modeling of sparse datasets that describe properties of social networks.

This project aims to develop a sufficiently fast and memory-efficient method to partition a hypergraph into roughly equal parts in such a way that the number of vertices in different 
partitions connected by the same hyperedge is minimized, and it is also paramount that this method can run with a streaming model. 
It is understandably harder to work on a dataset when it is streaming in real-time instead of being already recorded down, since 
the data points are seen only once and memory is far more restricted.

With a way to partition large datasets in hypergraph form, one can distribute the workload onto different devices, 
possibly eliminate the needless duplication of data and improve both the yield and the pace of their work. 
The practice of effectively partitioning large graphs is crucial in real-world applications like minimizing communication cost or I/O cost and maximization of robustness in social graphs.

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
