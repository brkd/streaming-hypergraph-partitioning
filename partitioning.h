#ifndef _PARTITIONING_H
#define _PARTITIONING_H

#include "bloomclass.h"
#include<vector>
#include <string>

class Algorithms
{
private:
	//Attributes
	int* sparseMatrix;
	int* sparseMatrixIndex;
	int** partitions;
	int** partitionToNet;
	std::vector<std::vector<int>> netToPartition;
	std::vector<int> readOrder;
	Bloom<int, int>* bloomFilter;

	double capacityConstraint;
	int partitionCount;
	int vertexCount;
	int edgeCount;
	int nonzeroCount;
	 
	//Methods
	int p2nConnectivity(int, int);
	int n2pIndex(int, int*);
	int BFConnectivity(int, int);
public:
	//Constructors
	Algorithms(std::string, int);
	Algorithms(std::string, int, int, int);
	//Destructor
	~Algorithms();

	//Methods
	void LDGp2n();
	void LDGn2p();
	void LDGBF();
	void LDGMultiBF();
	int calculateCuts();
};




#endif
