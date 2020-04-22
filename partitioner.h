#ifndef _PARTITIONER_H
#define _PARTITIONER_H

#include "bloomclass.h"
#include "bf_2.hpp"
#include "bf_3.hpp"
#include "bf_4.hpp"
#include "mlbfclass.hpp"
#include "mlbfclass_2.hpp"
#include <vector>
#include <string>

#define MAXPARTNO 4096
#define MAXIMBAL 1.20
#define MAXSLACKVAL 1000
#define INITVECSIZE 5

struct VertexInfo
{
  int ID;
  int netCount;
  int* netInfo;
VertexInfo(int ID, int netCount): ID(ID), netCount(netCount)
  {
    netInfo = new int[netCount]; 
  }
  
};

class Partitioner
{
 private:
  //Attributes
  std::string fileName;
  
  int* sparseMatrix;
  int* sparseMatrixIndex;
  int* reverse_sparseMatrix;
  int* reverse_sparseMatrixIndex;
  int* partVec; //pv[i] -> where i resides
  int* cutArray;
  
  int partitionCount;
  size_t vertexCount;
  size_t edgeCount;
  size_t nonzeroCount;
  size_t realVertexCount;
  int byteSize;
  int hashCount;
  bool symmetry;
  int noLayers;
  
  //Methods
  void read_binary_graph(std::string fileName);
  void read_mtx_and_transform_to_shpbin(std::string fileName);
  
  void LDGp2n(int, int, int, double);
  void LDGn2p(int, int, int, double);
  void LDGn2p_i(int, int, int, double, int);
  void LDGBF(int, int, int, double);
  void LDGBF2(int, int, int, double, int, int);
  void LDGBF3(int, int, int, double, int, int);
  void LDGBF4MULTI(int, int, int, double, int, int, int);
  void LDGBF5MULTI(int, int, int, double, int, int, int);
  void LDGMultiBF();

  void LDGn2p_ref(int, int, int, double, int);
  void LDGn2p_ref2(int, int, int, double, int);
  void LDGn2p_ref3(int, int, int, double, int);
  void n2pRefine(int*, int*, int, std::vector<std::vector<int>*>, double, int, std::vector<int>);
  void n2pRefine2(int, int, double, int*, int*, int*, bool*, std::vector<std::vector<int>*>&, const std::vector<int>&);
  void n2pRefine3(int,  double, const std::vector<VertexInfo*>&,  int*, int*, bool*, std::vector<std::vector<int>*>&, const std::vector<int>&, int&);

  //void LSH(int,double);
  void MinMax(int,int,int,double);
  void LSH(int, int);
  
  
  void vertexOutput(int, int);
  int calculateCuts(int);
  int calculateCuts2(int);
  int calculateCuts3(int, int, const std::vector<int>&, const std::vector<std::vector<int>*>&, int);
  
  int p2nConnectivity(int, int, const std::vector<std::vector<int>>&);
  int n2pIndex(int, int, double, int*, int*, bool*, const std::vector<std::vector<int>*>&, const std::vector<int>&);
  int n2pIndexAndScore(int, int, double, int*, int*, bool*, const std::vector<std::vector<int>*>&, const std::vector<int>&, double&);
  int BFConnectivity(Bloom<int, int>*, int, int);
  int BFConnectivity2(BloomFilter* bf, int, int);
  int BFConnectivity3(BloomFilter_OT* bf, int, int);
  int BFConnectivityMult(mlbf* bf, int, int);
  void BFConnectivityMult2(mlbf_2* bf, bool*, int*, int, int);
 public:
  //Constructors
  Partitioner(std::string);
 Partitioner(std::string name, int byteSize, int hashCount, int noLayers):Partitioner(name)
  {
    this->byteSize = byteSize;
    this->hashCount = hashCount;
    this->noLayers = noLayers;
  }
  void check_and_write_binary_graph(std::string fileName);
  
  //Destructor
  ~Partitioner();
  
  //Methods
  void partition(int, int, int, int, double, int, int, int);
  void RandomPartition(int, int, double, int);
};
#endif
