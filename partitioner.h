#ifndef _PARTITIONER_H
#define _PARTITIONER_H

#include "bloomclass.h"
#include "bf_2.hpp"
#include "bf_3.hpp"
#include <vector>
#include <string>

#define MAXPARTNO 4096
#define MAXIMBAL 1.20
#define MAXSLACKVAL 1000
#define INITVECSIZE 5

class Partitioner
{
 private:
  //Attributes
  int* sparseMatrix;
  int* sparseMatrixIndex;
  int* reverse_sparseMatrix;
  int* reverse_sparseMatrixIndex;
  int* partVec; //pv[i] -> where i resides
  double* scoreArray;
  
  int partitionCount;
  int vertexCount;
  int edgeCount;
  int nonzeroCount;
  int byteSize;
  int hashCount;
  bool symmetry;
  
  //Methods
  void read_binary_graph(std::string fileName);
  void read_mtx_and_transform_to_shpbin(std::string fileName);
  
  void LDGp2n(int, int, int, double);
  void LDGn2p(int, int, int, double);
  void LDGn2p_i(int, int, int, double);
  void LDGBF(int, int, int, double);
  void LDGBF2(int, int, int, double, int, int);
  void LDGBF3(int, int, int, double, int, int);
  void LDGMultiBF();
  
  void vertexOutput(int, int);
  int calculateCuts(int);
  int calculateCuts2(int);

  int p2nConnectivity(int, int, const std::vector<std::vector<int>>&);
  int n2pIndex(int, int, double, int*, int*, bool*, const std::vector<std::vector<int>*>&, const std::vector<int>&);
  int BFConnectivity(Bloom<int, int>*, int, int);
  int BFConnectivity2(BloomFilter* bf, int, int);
  int BFConnectivity3(BloomFilter_OT* bf, int, int);
 public:
  //Constructors
  Partitioner(std::string);
  Partitioner(std::string name, int byteSize, int hashCount):Partitioner(name)
  {
    this->byteSize = byteSize;
    this->hashCount = hashCount;
  }
  void check_and_write_binary_graph(std::string fileName);
  
  //Destructor
  ~Partitioner();
  
  //Methods
  void partition(int, int, int, int, double);
  void RandomPartition(int, int);
};
#endif
