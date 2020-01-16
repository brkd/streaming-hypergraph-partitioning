#ifndef _PARTITIONING_H
#define _PARTITIONING_H

#include "bloomclass.h"
#include <vector>
#include <string>

#define MAXPARTNO 4096
#define MAXIMBAL 1.20

class Algorithms
{
private:
	//Attributes
	int* sparseMatrix;
	int* sparseMatrixIndex;

	int** partitions;
	//int* partVector; //pv[i] -> where i resides

	std::vector<std::vector<int>> partitionToNet;

	std::vector<std::vector<int>> netToPartition;
	//int* netpointers;  - size: (net + 1) 
	//int* netparts; - size: pins - (netparts[netpointers[i]] to netparts[netpointers[i+1]-1 is allocated for the parts of net i. 
	//               - this array is initially assigned to -1 (net
	//int netcon; - size: net (initially 0) 

	std::vector<int> readOrder;
	Bloom<int, int>* bloomFilter;

	double capacityConstraint;
	int partitionCount;
	int vertexCount;
	int edgeCount;
	int nonzeroCount;
	 
	//Methods
	int p2nConnectivity(int, int);
	int n2pIndex(int, int*, int*, bool*, std::vector<int>);
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
