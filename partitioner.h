#ifndef _PARTITIONER_H
#define _PARTITIONER_H

#include "bloomclass.h"
#include <vector>
#include <string>

#define MAXPARTNO 4096
#define MAXIMBAL 1.20
#define MAXSLACKVAL 1000
#define INITVECSIZE 3

class Partitioner
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
  void read_graph(std::string fileName);

  void LDGp2n(int, double);
  void LDGn2p(int, int, double);
  void LDGBF(int, double);
  void LDGMultiBF();
  
  //int calculateCuts();
  int p2nConnectivity(int, int, const std::vector<std::vector<int>>&);
  int n2pIndex(int, int, double, int*, int*, bool*, const std::vector<std::vector<int>*>&, const std::vector<int>&);
  int BFConnectivity(int, int);
 public:
  //Constructors
  Partitioner(std::string);
  Partitioner(std::string, int, int);
  
  //Destructor
  ~Partitioner();
  
  //Methods
  void partition(int, int, int, double);
};
#endif
