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
  int* partVector; //pv[i] -> where i resides  

	Bloom<int, int>* bloomFilter;	
  
	int partitionCount;
	int vertexCount;
	int edgeCount;
	int nonzeroCount;
  
	//Methods
  void LDGp2n(int, double);
	void LDGn2p();
	void LDGBF();
	void LDGMultiBF();
 
	int calculateCuts();
	int p2nConnectivity(int, int);
	int n2pIndex(int, int*, int*, bool*, std::vector<int>);
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
