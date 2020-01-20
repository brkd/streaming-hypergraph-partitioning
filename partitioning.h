#ifndef _PARTITIONING_H
#define _PARTITIONING_H

#include "bloomclass.h"
#include <vector>
#include <string>

#define MAXPARTNO 4096
#define MAXIMBAL 1.20
#define INITVECSIZE 3

class Algorithms
{
private:
	//Attributes
	int* sparseMatrix;
	int* sparseMatrixIndex;
  int* partVec; //pv[i] -> where i resides  

	Bloom<int, int>* bloomFilter;	
  
	int partitionCount;
	int vertexCount;
	int edgeCount;
	int nonzeroCount;
  
	//Methods
  void LDGp2n(int, double);
	void LDGn2p(int, double);
	void LDGBF(int, double);
	void LDGMultiBF();
 
	//int calculateCuts();
	int p2nConnectivity(int, int, const std::vector<std::vector<int>>&);
	int n2pIndex(int, int, double, int*, int*, bool*, const std::vector<std::vector<int>*>&, const std::vector<int>&);
	int BFConnectivity(int, int);
public:
	//Constructors
	Algorithms(std::string);
	Algorithms(std::string, int, int);
 
	//Destructor
	~Algorithms();

	//Methods
	void partition(int, int, double);
};
#endif
