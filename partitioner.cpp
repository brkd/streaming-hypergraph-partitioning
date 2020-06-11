#include "partitioner.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <chrono>

#include <iomanip>
#include <stdio.h>

#include <omp.h>

//#define DEBUG
//#define WATCH

/*
int readBinaryGraph(FILE* bp, etype **pxadj, vtype **padj,
                    ewtype **pewghts, vwtype **pvwghts,
                    vtype* pnov) {
  fread(pnov, sizeof(vtype), 1, bp);
  (*pxadj) = (etype*)malloc(sizeof(etype) * (*pnov + 1));
  fread(*pxadj, sizeof(etype), (size_t)(*pnov + 1), bp);
  (*padj) = (vtype*)malloc(sizeof(vtype) * (*pxadj)[*pnov]);
  fread(*padj, sizeof(vtype), (size_t)(*pxadj)[*pnov], bp);
  (*pewghts) = (ewtype*)malloc(sizeof(ewtype) * (*pxadj)[*pnov]);
  fread(*pewghts, sizeof(ewtype), (size_t)(*pxadj)[*pnov], bp);
  (*pvwghts) = (vwtype*)malloc(sizeof(vwtype) * (*pnov));
  fread(*pvwghts, sizeof(vwtype), *pnov, bp);
  return 1;
}
int writeBinaryGraph(FILE* bp, etype *xadj, vtype *adj,
                     ewtype *ewghts, vwtype *vwghts,
                     vtype nov) {
  fwrite(&nov, sizeof(vtype), (size_t)1, bp);
  fwrite(xadj, sizeof(etype), (size_t)(nov + 1), bp);
  fwrite(adj, sizeof(vtype), (size_t)(xadj[nov]), bp);
  fwrite(ewghts, sizeof(ewtype), (size_t)(xadj[nov]), bp);
  fwrite(vwghts, sizeof(vwtype), (size_t)(nov), bp);
  return 1;
}

  sprintf(bfile, "%s_bin/%s.bin", currFolder, gfile + dirindex + 1);
  printf("Binary file name: %s\n", bfile);
  bp = fopen(bfile, "rb");
  if (bp != NULL) { // read from binary 
*/

//Public methods
Partitioner::Partitioner(std::string fileName){
  //Read matrix
  this->fileName = fileName;
  std::string mtx_name = fileName + ".mtx";
  std::string bin_name = fileName + ".shpbin";
  const char* bfile = bin_name.c_str();  
  FILE* bp;
  bp = fopen(bfile, "rb");
  
  
  if(bp == NULL){
    this->read_mtx_and_transform_to_shpbin(fileName);
    this->read_binary_graph(bin_name);    
  }
  else{
    fclose(bp);
    this->read_binary_graph(bin_name);
  }
  std::cout << "Size_info:" << this->vertexCount << ":" << this->edgeCount << ":" << this->nonzeroCount << std::endl;
  this->cutArray = new int[this->vertexCount];
  this->partVec = new int[this->vertexCount];
  for(int i = 0; i < vertexCount; i++)
    {
      cutArray[i] = -1;
      partVec[i] = -1;
    }
}

  
void Partitioner::read_binary_graph(std::string fileName){
  const char* fname = fileName.c_str();
  FILE* bp;
  bp = fopen(fname, "r");
  
  size_t* nov = new size_t;
  size_t* nnet = new size_t;
  size_t* nnz = new size_t;
  
  fread(nov, sizeof(int), 1, bp);
  this->vertexCount = *nov;

  fread(nnet, sizeof(int), 1, bp);
  this->edgeCount = *nnet;
  
  fread(nnz, sizeof(int), 1, bp);
  this-> nonzeroCount = *nnz;
  
#ifdef DEBUG
  std::cout << "No vertices: " << this->vertexCount << std::endl;
  std::cout << "No nets: " << this->edgeCount << std::endl;
  std::cout << "No nonzero: " << this-> nonzeroCount << std::endl;
#endif
  this->sparseMatrixIndex = new int[*nnet];
  fread(this->sparseMatrixIndex, sizeof(int), *nnet, bp);
  this->sparseMatrix = new int[*nnz];
  fread(this->sparseMatrix, sizeof(int), *nnz, bp);

  //Reverse reads

  this->reverse_sparseMatrixIndex = new int[*nov];
  fread(this->reverse_sparseMatrixIndex, sizeof(int), *nov, bp);
  
  this->reverse_sparseMatrix = new int[*nnz];
  fread(this->reverse_sparseMatrix, sizeof(int), *nnz, bp);
  
  size_t counter = 0;
  for(int i = 0; i < this->vertexCount - 1; i++)
    {
      if(this->reverse_sparseMatrixIndex[i] > this->reverse_sparseMatrixIndex[i + 1])
	break;
      counter++;
    }
  this->realVertexCount = counter;
  /*for(int i = this->reverse_sparseMatrixIndex[0]; i < this->reverse_sparseMatrixIndex[1]; i++){
    std::cout << "i: " << i  << " " <<this->reverse_sparseMatrix[i] << std::endl;
    }*/
  //

  
#ifdef DEBUG
  std::cout << "First indexes of xadj and adj" << std::endl;
  for(int i = 0; i < this->edgeCount - 1; i++){
    std::cout << "i: " << i  << " " <<this->sparseMatrixIndex[i] << " " << sparseMatrixIndex[i + 1] - sparseMatrix[i] << std::endl;
  }
  std::cout << "###############REVERSE#################" << std::endl;
  
  std::cout << "First indexes of reverse xadj and adj" << std::endl;
  for(int i = 0; i < this->vertexCount - 1; i++){
    std::cout << "i: " << i  << " " <<this->reverse_sparseMatrixIndex[i] << " " << reverse_sparseMatrixIndex[i + 1] - this->reverse_sparseMatrixIndex[i] << std::endl;
  }
  std::cout << "###############REVERSE#################" << std::endl;

#endif
  fclose(bp);
}

bool sortbyfirst_and_sec(const pair<int,int> &a, const pair<int,int> &b) 
{ 
  if(a.first != b.first)
    return a.first < b.first;
  else
    return a.second < b.second;
} 

bool sortbysec_and_first(const pair<int,int> &a, const pair<int,int> &b) 
{ 
  if(a.second != b.second)
    return a.second < b.second;
  else
    return a.first < b.first;
}

bool sortbyfirst(const pair<int,int> &a, const pair<int,int> &b) 
{
  return a.first < b.first;
}

bool sortbysec(const pair<int,int> &a, const pair<int,int> &b) 
{ 
  return a.second < b.second;
} 

void Partitioner::read_mtx_and_transform_to_shpbin(std::string fileName){

  std::string mtx_name = fileName + ".mtx";
  std::string bin_name = fileName + ".shpbin";
  const char* bfile = bin_name.c_str();

  int* xadj;
  int* adj;

  int* reverse_xadj;
  int* reverse_adj;

  int* no_vertex;
  int* no_nets;
  
  //FIELD//
  bool _real = false;
  bool _integer = false;
  bool _complex = false;
  bool _pattern = false;
  //FIELD//
  
  //SYMMETRY//
  bool _general = false;
  bool _symmetric = false;
  //SYMMETRY//
  
  //Evaluate type
  std::ifstream fin(mtx_name);
  std::string comment;
  std::getline(fin, comment);
  std::cout << "Type: " << comment << std::endl;
  
  if(comment.find("real") != std::string::npos)
    _real = true;
  
  if(comment.find("integer") != std::string::npos)
    _integer = true;
  
  if(comment.find("complex") != std::string::npos)
    _complex = true;
  
  if(comment.find("pattern") != std::string::npos)
    _pattern = true;
  
  if(comment.find("general") != std::string::npos)
    _general = true;
  
  if(comment.find("symmetric") != std::string::npos)
    _symmetric = true;
  
#ifdef DEBUG
  std::cout << "Matrix Type: " << std::endl;
  std::cout << "Real: " << _real << std::endl;
  std::cout << "Integer: " << _integer << std::endl;
  std::cout << "Complex: " << _complex << std::endl;
  std::cout << "Pattern: " << _pattern << std::endl;
  std::cout << "General: " << _general << std::endl;
  std::cout << "Symmetric: " << _symmetric << std::endl;
#endif


  while(fin.peek() == '%')
  {    
    fin.ignore(2048, '\n');
  }  

  //int vIndex = 0, row, col, currentColumn = -1;
  float val1;
  float val2;
  int val3;

  int no_row = 0;
  int no_col = 0;
  int nnz = 0;

  int vertex, net;

  fin >> no_row >> no_col >> nnz;
  std::vector<std::pair<int,int>> intermediate;
  
  if(_general){	
    
    if(_real || _integer){
      
      while(!fin.eof())
	{      
	  fin >> vertex >> net >> val1;
	  intermediate.push_back(std::pair<int, int>(vertex-1, net-1));
	}
    }
    
    if(_integer){
      
      while(!fin.eof())
	{      
	  fin >> vertex >> net >> val3;
	  intermediate.push_back(std::pair<int, int>(vertex-1, net-1));
	}
      
    }    
    
    if(_pattern){
      
      while(!fin.eof())
	{      
	  fin >> vertex >> net;
	  intermediate.push_back(std::pair<int, int>(vertex-1, net-1));
	}
      
    }
    
    if(_complex){      
      
      while(!fin.eof())
	{      
	  fin >> vertex >> net >> val1 >> val2;
	  intermediate.push_back(std::pair<int, int>(vertex-1, net-1));
	}
      
    }
    
    std::cout << "Sorting.. " << std::endl;
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbyfirst);
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbysec);
    no_nets = &intermediate[intermediate.size()-1].second;
    adj = new int[intermediate.size()];
    xadj = new int[no_col+1];
    xadj[0] = 0;
    std::cout << *no_nets << std::endl;
    int current_net = 0;
    int xadj_cursor = 1;
    for(int i = 0; i < intermediate.size(); i++){
      adj[i] = intermediate[i].first;
      if(intermediate[i].second != current_net){
	xadj[xadj_cursor] = i;
	xadj_cursor++;
	current_net = intermediate[i].second;
      }
    }
#ifdef DEBUG
    for(int i = 0; i < 100; i ++){
      
      std::cout << "i: " << i << " ||| xadj: " << xadj[i] << " adj: " << adj[i] << std::endl;
    }

    std::cout << intermediate.size() << " ||| " << nnz+1 << std::endl;
    std::cout << intermediate[intermediate.size()-1].second << " ||| " << no_col << std::endl;
#endif


    //REVERSE ORDER
    std::cout << "Reverse sorting.. " << std::endl;
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbysec);
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbyfirst);

#ifdef DEBUG
    std::cout << "AFTER REVERSE SORTING, INTERMEDIATE: " << std::endl;
    std::cout << this->nonzeroCount;
    for(int i = 0; i < nnz; i++){
      std::cout << "i: " << i << " ||-|| " << intermediate[i].first << " " << intermediate[i].second << std::endl; 
    }

#endif
    no_vertex = &intermediate[intermediate.size()-1].first;
    reverse_adj = new int[intermediate.size()];
    reverse_xadj = new int[no_row+1];
    reverse_xadj[0] = 0;
    int current_vertex = 0;
    int r_xadj_cursor = 1;
    std::cout << "START HERE " << std::endl;
    for(int i = 0; i < intermediate.size(); i++){
      reverse_adj[i] = intermediate[i].second;

      if(intermediate[i].first != current_vertex){
	reverse_xadj[r_xadj_cursor++] = i;
	current_vertex = intermediate[i].first;
      }
    }
#ifdef DEBUG
    for(int i = 0; i < 10; i ++){
      std::cout << "i: " << i << " ||| reverse_xadj: " << reverse_xadj[i] << " adj: " << reverse_xadj[i + 1] - reverse_xadj[i] << std::endl;
    }
    
#endif


    //REVERSE ORDER
    
  }
  
  if(_symmetric){
    
    if(_real || _integer){
      
      while(!fin.eof())
	{      
	  fin >> vertex >> net >> val1;
	  intermediate.push_back(std::pair<int, int>(vertex-1, net-1));
	  intermediate.push_back(std::pair<int, int>(net-1, vertex-1));
	}
    }
    
    if(_pattern){

      while(!fin.eof())
	{      
	  fin >> vertex >> net;
	  intermediate.push_back(std::pair<int, int>(vertex-1, net-1));
	  intermediate.push_back(std::pair<int, int>(net-1, vertex-1));	  
	}
      
    }
    
    if(_complex){      
      
      while(!fin.eof())
	{      
	  fin >> vertex >> net >> val1 >> val2;
	  intermediate.push_back(std::pair<int, int>(vertex-1, net-1));
	  intermediate.push_back(std::pair<int, int>(net-1, vertex-1));
	}
      
    }
    
    std::cout << "Sorting.. " << std::endl;
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbyfirst);
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbysec);

        
    for(int i = 0; i < intermediate.size(); i++){
      if((intermediate[i].first == intermediate[i+1].first) && intermediate[i].second == intermediate[i+1].second)
	intermediate.erase(intermediate.begin()+i);
      //std::cout << "Deduplicating " << i << "/" << intermediate.size() << std::endl;
    }
    
    no_nets = &intermediate[intermediate.size()-1].second;
    adj = new int[intermediate.size()];
    xadj = new int[no_col+1];
    
    xadj[0] = 0;
    int current_net = 0;
    int xadj_cursor = 1;
    
    for(int i = 0; i < intermediate.size(); i++){
      
      if(intermediate[i].first != intermediate[i+1].first ){
	adj[i] = intermediate[i].first;
	
	if(intermediate[i].second != current_net){
	  xadj[xadj_cursor] = i;
	  xadj_cursor++;
	  current_net = intermediate[i].second;
	}
      }
    }
    
#ifdef DEBUG
    for(int i = 0; i < 100; i ++){
      
      std::cout << "i: " << i << " ||| xadj: " << xadj[i] << " adj: " << adj[i] << std::endl;
    }
    
    std::cout << intermediate.size() << " ||| " << nnz+1 << std::endl;
    std::cout << intermediate[intermediate.size()-1].second << " ||| " << no_col << std::endl;
    std::cout << intermediate[intermediate.size()].second << " ||| " << no_col << std::endl;
#endif
    
    //REVERSE ORDER
    /*std::cout << "Reverse sorting.. " << std::endl;
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbysec);
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbyfirst);*/
    
    std::cout << "here2" << std::endl;
    no_vertex = &intermediate[intermediate.size()-1].second;
    reverse_adj = new int[intermediate.size()];
    reverse_xadj = new int[no_row+1];
    
    reverse_xadj[0] = 0;
    
    int current_vertex = 0;
    int r_xadj_cursor = 1;

    for(int i = 0; i < intermediate.size(); i++){
      adj[i] = intermediate[i].second;
      
      if(intermediate[i].first != current_vertex){
	reverse_xadj[r_xadj_cursor] = i;
	r_xadj_cursor++;
	current_vertex = intermediate[i].first;
      }
    }

#ifdef DEBUG
    for(int i = 0; i < 100; i ++){
      
      std::cout << "i: " << i << " ||| reverse_xadj: " << reverse_xadj[i] << " adj: " << reverse_adj[i] << std::endl;
    }

    std::cout << intermediate.size() << " ||| " << nnz+1 << std::endl;
    std::cout << intermediate[intermediate.size()-1].second << " ||| " << no_col << std::endl;
#endif

    //REVERSE ORDER

  }
  
  if(!_general && !_symmetric)
  {
    std::cout << "I believe a problem happened during reading fields of the matrix.." << std::endl;
    exit(1);
  }    
  
  fin.close();
  
  
  int np_nz = intermediate.size();
  int* no_nz;
  no_nz = &np_nz;

  no_vertex = &no_row;
  no_nets = &no_col;
  
  const char* fname = bin_name.c_str();
  FILE* bp;
  bp = fopen(fname, "w");
  
#ifdef DEBUG
  std::cout << "No vertex: " << *no_vertex << " No nets: " << *no_nets << " No nz: " << *no_nz << std::endl;
#endif

  fwrite(no_vertex, sizeof(int), 1, bp);
  fwrite(no_nets, sizeof(int), 1, bp);
  fwrite(no_nz, sizeof(int), 1, bp);
  fwrite(xadj, sizeof(int), *no_nets-1, bp);
  fwrite(adj, sizeof(int), *no_nz, bp);
  fwrite(reverse_xadj, sizeof(int), *no_vertex-1, bp);
  fwrite(reverse_adj, sizeof(int), *no_nz, bp);
  
  fclose(bp);

  std::cout << "Wrote the bin.. " << std::endl;
  
  delete[] xadj;
  delete[] adj;

  delete[] reverse_xadj;
  delete[] reverse_adj;

}


Partitioner::~Partitioner()
{
  delete[] this->partVec;
  delete[] this->cutArray;
  delete[] this->sparseMatrix;
  delete[] this->sparseMatrixIndex;
}

void Partitioner::partition(int algorithm, int partitionCount, int slackValue, int seed, double imbal, int i, int finalRun, int refSize)
{
  std::cout << "Started partitioning" << std::endl;
  //Partition
  if(algorithm == 1)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGp2n(partitionCount, slackValue, seed, imbal);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
    }
  else if(algorithm == 2)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGn2p(partitionCount, slackValue, seed, imbal);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
      for(int i = 0; i < this->vertexCount; i++)
	partVec[i] = -1;
      if(refSize > 0)
	{
	  std::cout << "----------------------------" << std::endl;
	  start = std::chrono::high_resolution_clock::now();
	  this->LDGn2p_ref(partitionCount, slackValue, seed, imbal, refSize);
	  end = std::chrono::high_resolution_clock::now();
	  std::cout << "Duration_R1:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
	  for(int i = 0; i < this->vertexCount; i++)
	    partVec[i] = -1;

	  std::cout << "----------------------------" << std::endl;
	  start = std::chrono::high_resolution_clock::now();
	  this->LDGn2p_ref2(partitionCount, slackValue, seed, imbal, refSize);
	  end = std::chrono::high_resolution_clock::now();
	  std::cout << "Duration_R2:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
	  for(int i = 0; i < this->vertexCount; i++)
	    partVec[i] = -1;

	  std::cout << "----------------------------" << std::endl;
	  start = std::chrono::high_resolution_clock::now();
	  this->LDGn2p_ref3(partitionCount, slackValue, seed, imbal, refSize);
	  end = std::chrono::high_resolution_clock::now();
	  std::cout << "Duration_R3:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
	  for(int i = 0; i < this->vertexCount; i++)
	    partVec[i] = -1;
	}
    }
  else if(algorithm == 3)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGn2p_i(partitionCount, slackValue, seed, imbal, i);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }
  else if(algorithm == 4)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGBF(partitionCount, slackValue, seed, imbal);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }
  else if(algorithm == 5)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGBF2(partitionCount, slackValue, seed, imbal, byteSize, hashCount);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }
  else if(algorithm == 6)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGBF3(partitionCount, slackValue, seed, imbal, byteSize, hashCount);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }
  else if(algorithm == 7)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGBF4MULTI(partitionCount, slackValue, seed, imbal, byteSize, hashCount, noLayers);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }
  else if(algorithm == 8)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LDGBF5MULTI(partitionCount, slackValue, seed, imbal, byteSize, hashCount, noLayers);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }
  else if(algorithm == 9)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->MinMax(partitionCount, slackValue, seed, imbal);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }
    else if(algorithm == 10)
    {
      auto start = std::chrono::high_resolution_clock::now();
      this->LSH(partitionCount, seed);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << std::endl;
    }

  
  //std::cout << "Cuts:" << this->calculateCuts(partitionCount) << std::endl;
  //std::cout << "Cuts:" << this->calculateCuts2(partitionCount) << std::endl;
  /*if (finalRun == 1)
    this->vertexOutput(algorithm, seed);  */

  //compute cut and report  
}

//ALGO 0//
void Partitioner::RandomPartition(int partitionCount, int seed, double imbal, int slack)
{
  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
    }
  
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, partitionCount - 1);
  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  int currVertexCount = 0;
  double capacityConstraint;
  for (int i : readOrder)
    {
      if((imbal*currVertexCount) >= slack)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slack) / partitionCount;      
      int partition = distribution(generator);
      while(sizeArray[partition] >= capacityConstraint)
	partition = distribution(generator);
      sizeArray[partition] += 1;
      partVec[i] = partition;
      //calculateCuts3(partitionCount, i);
      currVertexCount++;
    }
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  //}
  
  std::cout << "Cuts:" << this->calculateCuts2(partitionCount) << std::endl;
  
  this->vertexOutput(0, seed);
  delete[] sizeArray;
}

//ALGO 1//
void Partitioner::LDGp2n(int partitionCount, int slackValue, int seed, double imbal)
{
  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      //std::cout << "SIZE ARRAY[i]:" << sizeArray[i] << std::endl;
    }
  
  //Generate random read order
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  std::vector<std::vector<int>> partitionToNet(partitionCount);  
  
  int ctr = -1;
  
  double capacityConstraint;
  int currVertexCount = 0;
  for (int i : readOrder)
    {
      if((imbal*currVertexCount) >= slackValue)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slackValue) / partitionCount;
      double maxScore = -1.0;
      int maxIndex = -1;
      for (int j = 0; j < partitionCount; j++)
	{
	  int connectivity = this->p2nConnectivity(j, i, partitionToNet);
	  //int connectivity = 0;
	  //std::cout <<  "Partition " << j << " connectivity: " << connectivity << std::endl;
	  double partOverCapacity = sizeArray[j] / capacityConstraint;
	  double penalty = 1 - partOverCapacity;
	  double score = penalty * connectivity;
	  //std::cout << "sizearray[j]: " << sizeArray[j] << " poc: " << partOverCapacity << " penalty: " << penalty << " score: " << score << std::endl;
	  if (score > maxScore)
	    {
	      maxScore = score;
	      maxIndex = j;
	    }
	  else if (score == maxScore)
	    {
	      if (sizeArray[j] < sizeArray[maxIndex])
		{
		  maxIndex = j;
		}
	    }
	}
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      //calculateCuts3(partitionCount, i);
      int maxIndexSize = partitionToNet[maxIndex].size() - 1;
      
      for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
	{
	  if(std::find (partitionToNet[maxIndex].begin(), partitionToNet[maxIndex].end(), this->sparseMatrix[k]) == partitionToNet[maxIndex].end())
	    partitionToNet[maxIndex].push_back(this->sparseMatrix[k]);
	}
      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif
      
    }
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;
  std::cout << "******PART SIZES*******" << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  //}
  
  delete[] sizeArray;
}

//ALGO 2//
void Partitioner::LDGn2p(int partitionCount, int slackValue, int seed, double imbal)
{
  int* sizeArray = new int[partitionCount];
  int* indexArray = new int[partitionCount];
  bool* markerArray = new bool[partitionCount];
  
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      indexArray[i] = -1;
      markerArray[i] = false;
    }
  
  std::vector<int> readOrder;
  
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  std::vector<std::vector<int>*> netToPartition;  
  std::vector<int> tracker(10000, -1);
  double capacityConstraint;
  int currVertexCount = 0;
  
  for (int i : readOrder) {    
    if((imbal*currVertexCount) >= slackValue)
      capacityConstraint = (imbal*currVertexCount) / partitionCount;
    else
      capacityConstraint = ((double)slackValue) / partitionCount;
    for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->reverse_sparseMatrix[k];
	if(edge >= tracker.size())
	  {
	    int currNetIndex = tracker.size() - 1;
	    for(int j = currNetIndex; j < edge; j++)
	      {
	        tracker.push_back(-1);
		if(j == edge - 1)
		  {
		    std::vector<int>* newEdge = new std::vector<int>();
		    netToPartition.push_back(newEdge);
		    int n2pSize = netToPartition.size();
		    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
		    tracker[j] = n2pSize - 1;    	      
		  }            
	      }
	  }
	if (tracker[edge] == -1)
	  {
	    std::vector<int>* newEdge = new std::vector<int>();
	    netToPartition.push_back(newEdge);
	    int n2pSize = netToPartition.size();
	    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
	    tracker[edge] = n2pSize - 1;        
	  }
      }
    int maxIndex = this->n2pIndexCorrected(i, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker); 
    //std::cout << "Vertex: " << i << " cuts: " << calculateCuts2(partitionCount) << std::endl;
    partVec[i] = maxIndex;
    sizeArray[maxIndex] += 1;
    //calculateCuts3(partitionCount, i);
    for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->reverse_sparseMatrix[k];
	if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
	  netToPartition[tracker[edge]]->push_back(maxIndex);      
      }
    
    for (int i = 0; i < partitionCount; i++) {
      indexArray[i] = -1;
      markerArray[i] = false;
    }
    currVertexCount++;
#ifdef WATCH
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
    if(((double)currVertexCount/this->vertexCount)*100 < double(100.9)){
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
    }else{
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << "???" << "%" << std::flush;;
    }
#endif
  }
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;

  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  //}

  int reg_cuts = this->calculateCuts4(partitionCount, netToPartition);
  std::cout << "Cuts:" << reg_cuts << std::endl;
  std::cout << "Final Constraint: " << capacityConstraint << std::endl;
  for(int i = 0; i < netToPartition.size(); i++)
    {
      delete netToPartition[i];
    }

  delete[] sizeArray;
  delete[] indexArray;
  delete[] markerArray;
}

void Partitioner::LDGn2p_i(int partitionCount, int slackValue, int seed, double imbal, int ii)
{
  int* sizeArray = new int[partitionCount];
  int* indexArray = new int[partitionCount];
  bool* markerArray = new bool[partitionCount];
  
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      indexArray[i] = -1;
      markerArray[i] = false;
    }
  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  std::vector<std::vector<int>*> netToPartition;  
  std::vector<int> tracker(10000, -1);
  double capacityConstraint;
  int currVertexCount = 0;
  
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0,ii - 1);
  for(int i : readOrder)
    {
      if((imbal*currVertexCount) >= slackValue)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slackValue) / partitionCount;
      for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];
	  if(edge >= tracker.size())
	    {
	      int currNetIndex = tracker.size() - 1;
	      for(int j = currNetIndex; j < edge; j++)
		{
		  tracker.push_back(-1);
		  if(j == edge - 1)
		    {
		      std::vector<int>* newEdge = new std::vector<int>();
		      netToPartition.push_back(newEdge);
		      int n2pSize = netToPartition.size();
		      netToPartition[n2pSize - 1]->reserve(ii);
		      tracker[j] = n2pSize - 1;    	      
		    }            
		}
	    }
	  
	  if (tracker[edge] == -1)
	    {
	      std::vector<int>* newEdge = new std::vector<int>();
	      netToPartition.push_back(newEdge);
	      int n2pSize = netToPartition.size();
	      netToPartition[n2pSize - 1]->reserve(ii);
	      tracker[edge] = n2pSize - 1;
	    }
	}
      
      int maxIndex = this->n2pIndex(i, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker);
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      //calculateCuts3(partitionCount, i);
      for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];      
	  if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
	    {
	      if(netToPartition[tracker[edge]]->size() != ii)
		netToPartition[tracker[edge]]->push_back(maxIndex);    
	      else
		{
		  int newIndex = distribution(generator);
		  netToPartition[tracker[edge]]->at(newIndex) = maxIndex;
		}
	    }              
	}
      for (int i = 0; i < partitionCount; i++) {
	indexArray[i] = -1;
	markerArray[i] = false;
      }
      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif
      
    }
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  //}
  for(int i = 0; i < netToPartition.size(); i++)
    {
      delete netToPartition[i];
    }
  delete[] sizeArray;
  delete[] indexArray;
  delete[] markerArray;
}


void Partitioner::LDGBF(int partitionCount, int slackValue, int seed, double imbal)
{
  Bloom<int, int>* bloomFilter = new Bloom<int, int>(this->byteSize, this->hashCount);
  
  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
    }
  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  double capacityConstraint;
  int currVertexCount = 0;
  for (int i : readOrder)
    {
      if((imbal*currVertexCount) >= slackValue)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slackValue) / partitionCount;
      
      double maxScore = -1.0;
      int maxIndex = -1;
      
      for (int j = 0; j < partitionCount; j++)
	{
	  int connectivity = this->BFConnectivity(bloomFilter, j, i);
	  double partToCapacity = sizeArray[j] / capacityConstraint;
	  double penalty = 1 - partToCapacity;
	  double score = penalty * connectivity;
	  if (score > maxScore)
	    {
	      maxScore = score;
	      maxIndex = j;
	    }
	  else if (score == maxScore)
	    {
	      if (sizeArray[j] < sizeArray[maxIndex])
		{
		  maxIndex = j;
		}
	    }
	}
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      //calculateCuts3(partitionCount, i);
      
      for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];
	  bloomFilter->insert(edge, maxIndex);
	}
      
      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif
      
      
      
    }
  
  std::cout << std::endl;
  
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  //}
  
  delete[] sizeArray;
}


void Partitioner::LDGBF2(int partitionCount, int slackValue, int seed, double imbal, int byteSize, int hashCount)
{
  int* sizeArray = new int[partitionCount];
  BloomFilter* bf = new BloomFilter[partitionCount];
  
  for(int i = 0; i < partitionCount; i++){
    bf[i] = BloomFilter(byteSize*8, hashCount, i);
  }
  
  
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
    }
  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
  {
    readOrder.push_back(i);
  }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  double capacityConstraint;
  int currVertexCount = 0;
  for (int i : readOrder)
    {
      if((imbal*currVertexCount) >= slackValue)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slackValue) / partitionCount;
      
      double maxScore = -1.0;
      int maxIndex = -1;
      for (int j = 0; j < partitionCount; j++)
	{
	  int connectivity = this->BFConnectivity2(bf, j, i);
	  double partToCapacity = sizeArray[j] / capacityConstraint;
	  double penalty = 1 - partToCapacity;
	  double score = penalty * connectivity;
	  if (score > maxScore)
	    {
	      maxScore = score;
	      maxIndex = j;
	    }
	  else if (score == maxScore)
	    {
	      if (sizeArray[j] < sizeArray[maxIndex])
		{
		  maxIndex = j;
		}
	    }
	}
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      //calculateCuts3(partitionCount, i);
      for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	{
       	  int edge = this->reverse_sparseMatrix[k];
	  bf[maxIndex].insert(edge);
	}

      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif

      
      
    }

  std::cout << std::endl;
  
  std::cout << "******PART SIZES*******" << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  //}

  delete[] bf;
  delete[] sizeArray;
}

void Partitioner::LDGBF3(int partitionCount, int slackValue, int seed, double imbal, int byteSize, int hashCount)
{
  //Bloom<int, int>* bloomFilter = new Bloom<int, int>(this->byteSize, this->hashCount);
  
  BloomFilter_OT* bf = new BloomFilter_OT(byteSize*8, hashCount, 5);
  
  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
    }
  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  double capacityConstraint;
  int currVertexCount = 0;
  for (int i : readOrder)
    {
      if((imbal*currVertexCount) >= slackValue)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slackValue) / partitionCount;
      
      double maxScore = -1.0;
      int maxIndex = -1;
      
      for (int j = 0; j < partitionCount; j++)
	{
	  int connectivity = this->BFConnectivity3(bf, j, i);
	  double partToCapacity = sizeArray[j] / capacityConstraint;
	  double penalty = 1 - partToCapacity;
	  double score = penalty * connectivity;
	  if (score > maxScore)
	    {
	      maxScore = score;
	      maxIndex = j;
	    }
	  else if (score == maxScore)
	    {
	      if (sizeArray[j] < sizeArray[maxIndex])
		{
		  maxIndex = j;
		}
	    }
	}
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      //calculateCuts3(partitionCount, i);
      for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];
	  bf->insert(edge, maxIndex);
	}
      
      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif
      
      
      
    }
  
  std::cout << std::endl;
  
  std::cout << "******PART SIZES*******" << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  //}
  
  delete[] sizeArray;
}

void Partitioner::LDGBF4MULTI(int partitionCount, int slackValue, int seed, double imbal, int byteSize, int hashCount, int noLayers)
{
  //Bloom<int, int>* bloomFilter = new Bloom<int, int>(this->byteSize, this->hashCount);
  
  //BloomFilter_OT* bf = new BloomFilter_OT(byteSize*8, hashCount, 5);
  mlbf* bf = new mlbf(noLayers, partitionCount, hashCount, byteSize);
  std::cout << "Initialized" << std::endl;

  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
    }
  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  double capacityConstraint;
  int currVertexCount = 0;
  for (int i : readOrder)
    {
      if((imbal*currVertexCount) >= slackValue)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slackValue) / partitionCount;
      
      double maxScore = -1.0;
      int maxIndex = -1;
      
      for (int j = 0; j < partitionCount; j++)
	{
	  int connectivity = this->BFConnectivityMult(bf, j, i);
	  double partToCapacity = sizeArray[j] / capacityConstraint;
	  double penalty = 1 - partToCapacity;
	  double score = penalty * connectivity;
	  if (score > maxScore)
	    {
	      maxScore = score;
	      maxIndex = j;
	    }
	  else if (score == maxScore)
	    {
	      if (sizeArray[j] < sizeArray[maxIndex])
		{
		  maxIndex = j;
		}
	    }
	}
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      //calculateCuts3(partitionCount, i);
      for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];
	  bf->insert(edge, maxIndex);
	}
      
      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif
      
      
      
    }
  
  std::cout << std::endl;
  
  std::cout << "******PART SIZES*******" << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //   std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  //}
  
  delete[] sizeArray;
}


///////ENHANCED MULTI LAYER +9 GODLY////////
void Partitioner::LDGBF5MULTI(int partitionCount, int slackValue, int seed, double imbal, int byteSize, int hashCount, int noLayers)
{
  mlbf_2* bf = new mlbf_2(noLayers, partitionCount, hashCount, byteSize);
  int no_leaf_blocks = pow(2, noLayers);
  int leaf_block_size = partitionCount/no_leaf_blocks;
  bool* existences = new bool[no_leaf_blocks];
  //
  int* connectivities = new int[no_leaf_blocks];
  //
  std::pair<int, int>* ranges = new std::pair<int,int>[no_leaf_blocks];
  for(int i = 0; i < no_leaf_blocks; i++){
    int start = i*leaf_block_size;
    int end = (i+1)*leaf_block_size;
    ranges[i].first = start;
    ranges[i].second = end;
  }
  std::cout << "Initialized" << std::endl;

  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
    }
  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  double capacityConstraint;
  int currVertexCount = 0;
  for (int i : readOrder)
    {

      for(int re = 0; re < no_leaf_blocks; re++){
	connectivities[re] = 0;
      }
      
      
      if((imbal*currVertexCount) >= slackValue)
	capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
	capacityConstraint = ((double)slackValue) / partitionCount;
      
      double maxScore = -1.0;
      int maxIndex = -1;
      
      //int connectivity = this->BFConnectivityMult2(bf, existences, connectivities, i);
      this->BFConnectivityMult2(bf, existences, connectivities, no_leaf_blocks, i);

      
      for(int j = 0; j < no_leaf_blocks; j++){
	int connectivity = connectivities[j];
	
	
	int average_size = 0;
	for(int ps = ranges[j].first; ps < ranges[j].second; ps++){
	  average_size += sizeArray[ps];
	}

	average_size /= leaf_block_size;
	
	double partToCapacity = average_size / capacityConstraint;

	double penalty = 1 - partToCapacity;
	double score = penalty * connectivity;
	
	
      if (score > maxScore)
	{
	  maxScore = score;
	  maxIndex = j;
	}
      else if (score == maxScore)
	{
	  if (sizeArray[j] < sizeArray[maxIndex])
	    {
	      maxIndex = j;
	    }
	}
      }
      int part = rand() % leaf_block_size + (maxIndex*leaf_block_size);
      partVec[i] = part;
      sizeArray[part] += 1;
      //calculateCuts3(partitionCount, i);
      for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];
	  bf->insert(edge, part);
	}
      
      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif
      
      
      
    }
  
  std::cout << std::endl;
  
  std::cout << "******PART SIZES*******" << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  //}
  
  delete[] sizeArray;
  delete[] existences;
  delete[] connectivities;
}


void Partitioner::LDGMultiBF()
{
  /*	int* sizeArray = new int[this->partitionCount];
	for (int i = 0; i < this->partitionCount; i++)
	{
		sizeArray[i] = 0;
	}
	int partitions = this->partitionCount;
	int layers = 0;
	while (partitions != 0)
	{
		layers++;
		partitions = partitions >> 2;
	}
	int filterCount = pow(2, layers + 1) - 1;
	std::vector<Bloom*> filters;
	int bytes = 1000000;
	for (int i = 1; i < filterCount; i++)
	{
		filters[i] = new Bloom(bytes, this->hashCount);
	}
	double maxScore = -1.0;
	int maxIndex = -1;
	for (int i : this->readOrder)
	{
		int* connectivities = new int[this->partitionCount];
		for (int p = 0; p < this->partitionCount; p++)
		{
			connectivities[p] = 0;
		}
		for (int j = this->sparseMatrixIndex[i]; j < this->sparseMatrixIndex[i + 1]; j++)
		{
			for (int k = 1; k < filterCount; k++)
			{
				if (filters[k]->contains(j, 0))
				{
					if (k * 2 >= filterCount)
					{
						connectivities[k % this->partitionCount]++;//Final layer, partition found
					}
					else
					{
						k = (k * 2) + 1; //Next filter on the path
					}
				}
			}
		}
		delete[] connectivities;
	}
	
	//Might use a stack in order not to skip partitions
	//Insertion?
	*/
}

void Partitioner::LDGn2p_ref(int partitionCount, int slackValue, int seed, double imbal, int refSize)
{
  int* sizeArray = new int[partitionCount];
  int* indexArray = new int[partitionCount];
  bool* markerArray = new bool[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      indexArray[i] = -1;
      markerArray[i] = false;
    }
  int refCount = refSize * 250000;
  int* refTracker = new int[refCount];
  for (int i = 0; i < refCount; i++)
    {
      refTracker[i] = -1;
    }

  std::vector<int> readOrder;
  
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  std::vector<std::vector<int>*> netToPartition;  
  std::vector<int> tracker(10000, -1);
  double capacityConstraint;
  int currVertexCount = 0;
  int vertexCountToRefine = 0;
  for (int i : readOrder) {
    
    if((imbal*currVertexCount) >= slackValue)
      capacityConstraint = (imbal*currVertexCount) / partitionCount;
    else
      capacityConstraint = ((double)slackValue) / partitionCount;
    for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->reverse_sparseMatrix[k];
	if(edge >= tracker.size())
	  {
	    int currNetIndex = tracker.size() - 1;
	    for(int j = currNetIndex; j < edge; j++)
	      {
	        tracker.push_back(-1);
		if(j == edge - 1)
		  {
		    std::vector<int>* newEdge = new std::vector<int>();
		    netToPartition.push_back(newEdge);
		    int n2pSize = netToPartition.size();
		    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
		    tracker[j] = n2pSize - 1;    	      
		  }            
	      }
	  }
	
	if (tracker[edge] == -1)
	  {
	    std::vector<int>* newEdge = new std::vector<int>();
	    netToPartition.push_back(newEdge);
	    int n2pSize = netToPartition.size();
	    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
	    tracker[edge] = n2pSize - 1;        
	  }
      }
    int maxIndex = this->n2pIndexCorrected(i, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker); 
    //std::cout << "Vertex: " << i << " cuts: " << calculateCuts2(partitionCount) << std::endl;
    partVec[i] = maxIndex;
    sizeArray[maxIndex] += 1;
    //calculateCuts3(partitionCount, i);
    for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->reverse_sparseMatrix[k];      
	if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
	  netToPartition[tracker[edge]]->push_back(maxIndex);      
      }
    
    for (int i = 0; i < partitionCount; i++) {
      indexArray[i] = -1;
      markerArray[i] = false;
    }
    currVertexCount++;
    refTracker[vertexCountToRefine++] = i;
    if(vertexCountToRefine == refCount)
      {
	std::cout << "Starting refinement round..." << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	//this->n2pRefine(sizeArray, refTracker, refCount, netToPartition, capacityConstraint, partitionCount, tracker);
	this->n2pRefine2(partitionCount, refCount, capacityConstraint, refTracker, sizeArray, indexArray, markerArray, netToPartition, tracker);
	auto end  = std::chrono::high_resolution_clock::now();
	std::cout << "Refinement finished round..." << std::endl;
	std::cout << "Round duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
	for(int i = 0; i < refCount; i++)
	  {
	    refTracker[i] = -1;
	  }
	vertexCountToRefine = 0;
      }

#ifdef WATCH
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
    if(((double)currVertexCount/this->vertexCount)*100 < double(100.9)){
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
    }else{
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << "???" << "%" << std::flush;;
    }
#endif
  }
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;

  int r1_cuts = this->calculateCuts4(partitionCount, netToPartition);
  std::cout << "Cuts_R1:" << r1_cuts << std::endl;
  std::cout << "Final Constraint: " << capacityConstraint << std::endl;


  /*for(int i = 0; i < netToPartition.size(); i++)
    {
      delete netToPartition[i];
    }
    delete[] sizeArray;
    delete[] indexArray;
    delete[] markerArray;*/
}

void Partitioner::LDGn2p_ref2(int partitionCount, int slackValue, int seed, double imbal, int refSize)
{
  int* sizeArray = new int[partitionCount];
  int* indexArray = new int[partitionCount];
  bool* markerArray = new bool[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      indexArray[i] = -1;
      markerArray[i] = false;
    }
  int refCount = refSize * 250000;
  double refScore = -1.0;
  int* refTracker = new int[refCount];
  for (int i = 0; i < refCount; i++)
    {
      refTracker[i] = -1;
    }

  std::vector<int> readOrder;
  
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  std::vector<std::vector<int>*> netToPartition;  
  std::vector<int> tracker(10000, -1);
  double capacityConstraint;
  int currVertexCount = 0;
  int vertexCountToRefine = 0;
  for (int i : readOrder) {
    
    if((imbal*currVertexCount) >= slackValue)
      capacityConstraint = (imbal*currVertexCount) / partitionCount;
    else
      capacityConstraint = ((double)slackValue) / partitionCount;
    for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->reverse_sparseMatrix[k];
	if(edge >= tracker.size())
	  {
	    int currNetIndex = tracker.size() - 1;
	    for(int j = currNetIndex; j < edge; j++)
	      {
	        tracker.push_back(-1);
		if(j == edge - 1)
		  {
		    std::vector<int>* newEdge = new std::vector<int>();
		    netToPartition.push_back(newEdge);
		    int n2pSize = netToPartition.size();
		    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
		    tracker[j] = n2pSize - 1;    	      
		  }            
	      }
	  }
	if (tracker[edge] == -1)
	  {
	    std::vector<int>* newEdge = new std::vector<int>();
	    netToPartition.push_back(newEdge);
	    int n2pSize = netToPartition.size();
	    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
	    tracker[edge] = n2pSize - 1;        
	  }
      }
    int maxIndex = this->n2pIndexAndScoreCorrected(i, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker, refScore); 
    //std::cout << "Vertex: " << i << " cuts: " << calculateCuts2(partitionCount) << std::endl;
    partVec[i] = maxIndex;
    sizeArray[maxIndex] += 1;
    //calculateCuts3(partitionCount, i);
    for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->reverse_sparseMatrix[k];      
	if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
	  netToPartition[tracker[edge]]->push_back(maxIndex);      
      }
    
    for (int i = 0; i < partitionCount; i++) {
      indexArray[i] = -1;
      markerArray[i] = false;
    }
    currVertexCount++;
    if(refScore == 0.0 && !(this->reverse_sparseMatrixIndex[i + 1] <= this->reverse_sparseMatrixIndex[i]))
      {
	refTracker[vertexCountToRefine++] = i;
	refScore = -1.0;
      }

    if(vertexCountToRefine == refCount)
      {
	std::cout << "Starting refinement round..." << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	//this->n2pRefine(sizeArray, refTracker, refCount, netToPartition, capacityConstraint, partitionCount, tracker);
	this->n2pRefine2(partitionCount, refCount, capacityConstraint, refTracker, sizeArray, indexArray, markerArray, netToPartition, tracker);
	auto end  = std::chrono::high_resolution_clock::now();
	std::cout << "Refinement finished round..." << std::endl;
	std::cout << "Round duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
	for(int i = 0; i < refCount; i++)
	  {
	    refTracker[i] = -1;
	  }
	vertexCountToRefine = 0;
      }

#ifdef WATCH
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
    if(((double)currVertexCount/this->vertexCount)*100 < double(100.9)){
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
    }else{
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << "???" << "%" << std::flush;;
    }
#endif
  }

  for(int i = 0; i < vertexCountToRefine; i++)
    {
      std::cout << "Starting refinement round..." << std::endl;
      auto start = std::chrono::high_resolution_clock::now();
      //this->n2pRefine(sizeArray, refTracker, refCount, netToPartition, capacityConstraint, partitionCount, tracker);
      this->n2pRefine2(partitionCount, refCount, capacityConstraint, refTracker, sizeArray, indexArray, markerArray, netToPartition, tracker);
      auto end  = std::chrono::high_resolution_clock::now();
      std::cout << "Refinement finished round..." << std::endl;
      std::cout << "Round duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
      for(int i = 0; i < refCount; i++)
	{
	  refTracker[i] = -1;
	}
      vertexCountToRefine = 0;
    }
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;

  int r2_cuts = this->calculateCuts4(partitionCount, netToPartition);
  std::cout << "Cuts_R2:" << r2_cuts << std::endl;
  std::cout << "Final Constraint: " << capacityConstraint << std::endl;
  /*for(int i = 0; i < netToPartition.size(); i++)

    {
      delete netToPartition[i];
    }
  delete[] sizeArray;
  delete[] indexArray;
  delete[] markerArray;*/
}

void Partitioner::LDGn2p_ref3(int partitionCount, int slackValue, int seed, double imbal, int refSize)
{
  int* sizeArray = new int[partitionCount];
  int* indexArray = new int[partitionCount];
  bool* markerArray = new bool[partitionCount];
  std::vector<VertexInfo*> putAsideVector;
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      indexArray[i] = -1;
      markerArray[i] = false;
    }
  int refCount = refSize * 250000;
  double refScore = -1.0;

  std::vector<int> readOrder;
  
  /*for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
      }*/
  for (int i = 0; i < this->vertexCount; i++)
    readOrder.push_back(i);
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  std::vector<std::vector<int>*> netToPartition;  
  std::vector<int> tracker(10000, -1);
  double capacityConstraint;
  int currVertexCount = 0;
  int current_pa_size = 0;
  for (int i : readOrder) {
    
    if((imbal*currVertexCount) >= slackValue)
      capacityConstraint = (imbal*currVertexCount) / partitionCount;
    else
      capacityConstraint = ((double)slackValue) / partitionCount;
    for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->reverse_sparseMatrix[k];
	if(edge >= tracker.size())
	  {
	    int currNetIndex = tracker.size() - 1;
	    for(int j = currNetIndex; j < edge; j++)
	      {
	        tracker.push_back(-1);
		if(j == edge - 1)
		  {
		    std::vector<int>* newEdge = new std::vector<int>();
		    netToPartition.push_back(newEdge);
		    int n2pSize = netToPartition.size();
		    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
		    tracker[j] = n2pSize - 1;    	      
		  }            
	      }
	  }
	if (tracker[edge] == -1)
	  {
	    std::vector<int>* newEdge = new std::vector<int>();
	    netToPartition.push_back(newEdge);
	    int n2pSize = netToPartition.size();
	    netToPartition[n2pSize - 1]->reserve(INITVECSIZE);
	    tracker[edge] = n2pSize - 1;        
	  }
      }
    int maxIndex = this->n2pIndexAndScoreCorrected(i, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker, refScore); 

    if(refScore <= 0.0 && !(this->reverse_sparseMatrixIndex[i + 1] <= this->reverse_sparseMatrixIndex[i]))
      {
	VertexInfo* new_pav = new VertexInfo(i, this->reverse_sparseMatrixIndex[i + 1] - this->reverse_sparseMatrixIndex[i]);
	
	for(int offset = 0; offset < new_pav->netCount; offset++)
	  {
	    new_pav->netInfo[offset] = this->reverse_sparseMatrix[this->reverse_sparseMatrixIndex[i] + offset];
	  }
	putAsideVector.push_back(new_pav);
	current_pa_size += 1;
	if(current_pa_size >= refCount)
	  {
	    std::cout << "Starting refinement round..." << std::endl;
	    auto start = std::chrono::high_resolution_clock::now();
	    //refine burada
	    this->n2pRefine3(partitionCount, imbal, putAsideVector, sizeArray, indexArray, markerArray, netToPartition, tracker, currVertexCount);
	    auto end = std::chrono::high_resolution_clock::now();
	    std::cout << "Refinement round finished..." << std::endl;
	    std::cout << "Round duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
	    for(int i = 0; i < putAsideVector.size(); i++)
	      {
		delete putAsideVector[i];
	      }
	    putAsideVector.clear();

	    current_pa_size = 0;
	  }	
      }
    else
      {
	partVec[i] = maxIndex;
	sizeArray[maxIndex] += 1;
	for (int k = this->reverse_sparseMatrixIndex[i]; k < this->reverse_sparseMatrixIndex[i + 1]; k++)
	  {
	    int edge = this->reverse_sparseMatrix[k];      
	    if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
	      netToPartition[tracker[edge]]->push_back(maxIndex);      
	  }
	
	currVertexCount++;
      }
    for (int i = 0; i < partitionCount; i++) {
      indexArray[i] = -1;
      markerArray[i] = false;
    }

    //calculateCuts3(partitionCount, i);

#ifdef WATCH
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
    if(((double)currVertexCount/this->vertexCount)*100 < double(100.9)){
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
    }else{
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << "???" << "%" << std::flush;;
    }
#endif
  }
  
  if(putAsideVector.size() != 0)
    {
      std::cout << "Starting refinement round..." << std::endl;
      auto start = std::chrono::high_resolution_clock::now();
      //refine burada
      this->n2pRefine3(partitionCount, imbal, putAsideVector, sizeArray, indexArray, markerArray, netToPartition, tracker, currVertexCount);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Refinement round finished..." << std::endl;
      std::cout << "Round duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
      for(int i = 0; i < putAsideVector.size(); i++)
	{
	  delete putAsideVector[i];
	}
      putAsideVector.clear();
      current_pa_size = 0;
    }	
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;
  
  //for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  //}

  int r3_cuts = this->calculateCuts4(partitionCount, netToPartition);
  capacityConstraint = (imbal*currVertexCount) / partitionCount;
  std::cout << "Cuts_R3:" << r3_cuts << std::endl;
  std::cout << "Final Constraint: " << capacityConstraint << std::endl;
  /*for(int i = 0; i < netToPartition.size(); i++)
    {
      delete netToPartition[i];
      }*/
  //delete[] sizeArray;
  //delete[] indexArray;
  //delete[] markerArray;
}

void Partitioner::n2pRefine(int* sizeArray, int* refTracker, int refCount, std::vector<std::vector<int>*> netToPartition, double capacityConstraint, int partitionCount, std::vector<int> tracker)
{
  for(int i = 0; i < refCount; i++)
    {
      int vertex = refTracker[i];
      int cut = this->calculateCuts3(partitionCount, vertex, tracker, netToPartition, -1);
      if(cut == 0)
	continue;
      int new_part = this->partVec[vertex];
      int new_cut = cut;
      for(int j = reverse_sparseMatrixIndex[vertex]; j < reverse_sparseMatrixIndex[vertex + 1]; j++)
	{
	  int net = reverse_sparseMatrix[j];
	  for(int k = 0; k < netToPartition[tracker[net]]->size(); k++)
	    {
	      int candidate_part = netToPartition[tracker[net]]->at(k);

	      if(sizeArray[candidate_part] >= capacityConstraint)
		continue;
	      int cand_cut = this->calculateCuts3(partitionCount, vertex, tracker, netToPartition, candidate_part);
	      if(cand_cut < cut)
		{
		  new_part = candidate_part;
		  new_cut = cand_cut;
		}
	    }
	}
      this->partVec[vertex] = new_part;
    }
}

void Partitioner::n2pRefine2(int partitionCount, int refCount, double capacityConstraint, int* refTracker, int* sizeArray, int* indexArray, bool* markerArray, std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker)
{
  for(int i = 0; i < refCount; i++)
    {
      int vertex = refTracker[i];
      int maxIndex = this->n2pIndexCorrected(vertex, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker);

      if(maxIndex != this->partVec[vertex] && sizeArray[maxIndex] < capacityConstraint)
	{
	  int old_part = this->partVec[vertex];
	  sizeArray[old_part] -= 1;
	  
	  this->partVec[vertex] = maxIndex;
	  sizeArray[maxIndex] += 1;
	  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
	    {
	      int edge = this->reverse_sparseMatrix[k];
	      bool dc_part = true;
	      for(int j = this->sparseMatrixIndex[edge]; j < this->sparseMatrixIndex[edge + 1]; j++)
		{
		  int v = this->sparseMatrix[j];
		  if(this->partVec[v] == old_part)
		    {
		      dc_part = false;
		      break;
		    }
		  if(dc_part)
		    {
		      std::vector<int>::iterator position = std::find(netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), old_part);
		      if(position != netToPartition[tracker[edge]]->end())
			netToPartition[tracker[edge]]->erase(position);
		    }
		}
	      if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
		netToPartition[tracker[edge]]->push_back(maxIndex);      
	    }
	}

      for (int i = 0; i < partitionCount; i++) {
	indexArray[i] = -1;
	markerArray[i] = false;
      }      
    }
}

void Partitioner::n2pRefine3(int partitionCount, double imbal, const std::vector<VertexInfo*>& putAsideVector, int* sizeArray, int* indexArray, bool* markerArray, std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker, int& currVertexCount)
{
  for(int i = 0; i < putAsideVector.size(); i++)
    {
      double capacityConstraint = (imbal*currVertexCount) / partitionCount;
      int vertex = putAsideVector[i]->ID;
      int maxIndex = this->n2pIndexCorrected(vertex, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker); 

      this->partVec[vertex] = maxIndex;
      sizeArray[maxIndex] += 1;
      for(int offset = 0; offset < putAsideVector[i]->netCount; offset++)
	{
	  int edge = putAsideVector[i]->netInfo[offset];
	  if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
	    netToPartition[tracker[edge]]->push_back(maxIndex);
	}
      for (int i = 0; i < partitionCount; i++) {
	indexArray[i] = -1;
	markerArray[i] = false;
      }  
      currVertexCount++;	
    }
}
void Partitioner::MinMax(int partitionCount, int SLACK, int seed, double imbal)
{
  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      //std::cout << "SIZE ARRAY[i]:" << sizeArray[i] << std::endl;
    }
  
  //Generate random read order
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());

  std::vector<std::vector<int>> partitionToNet(partitionCount); 

  int currVertexCount = 0;
  double capacityConstraint;
  for (int i : readOrder)
    {
       if((imbal*currVertexCount) >= SLACK)
        capacityConstraint = (imbal*currVertexCount) / partitionCount;
      else
         capacityConstraint = ((double)SLACK) / partitionCount;

      double maxScore = -1.0;
      int maxIndex = -1;
      
      int minLoad = sizeArray[0];
      for(int a = 1; a < partitionCount; a++)
      {
        if(sizeArray[a] < minLoad)
          minLoad = sizeArray[a];
      }

for (int j = 0; j < partitionCount; j++)
  {
    if(sizeArray[j] <= minLoad + SLACK)
      {
        int connectivity = this->p2nConnectivity(j, i, partitionToNet);
        double partOverCapacity = sizeArray[j] / capacityConstraint;
        double penalty = 1 - partOverCapacity;
        double score = penalty * connectivity;

        if (score > maxScore)
          {
            maxScore = score;
            maxIndex = j;
          }
        else if (score == maxScore && sizeArray[j] < sizeArray[maxIndex])
          {
            maxIndex = j;
          }
      }
  }
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      //calculateCuts3(partitionCount, i);
      int maxIndexSize = partitionToNet[maxIndex].size() - 1;


  for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
  {
    if(std::find (partitionToNet[maxIndex].begin(), partitionToNet[maxIndex].end(), this->sparseMatrix[k]) == partitionToNet[maxIndex].end())
      partitionToNet[maxIndex].push_back(this->sparseMatrix[k]);
  }

      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif
      
    }
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;
  std::cout << "******PART SIZES*******" << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  //}
  
  delete[] sizeArray;
}

void Partitioner::LSH(int partitionCount, int seed)
{
  int* sizeArray = new int[partitionCount];
  for (int i = 0; i < partitionCount; i++)
    {
      sizeArray[i] = 0;
      //std::cout << "SIZE ARRAY[i]:" << sizeArray[i] << std::endl;
    }
  
  //Generate random read order
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());
  int currVertexCount = 0;
  for (int i : readOrder)
  {
    int hash1min = 4099;
    int hash2min = 2053;
    int hash3min = 1031;
    int minnet1 = 0;
    int minnet2 = 0;
    int minnet3 = 0;
    
    for(int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
  	 {  
  	 	int hash1 = (16 * this->sparseMatrix[k] + 25)%4099;
  	 	int hash2 = (27 * this->sparseMatrix[k] + 49)%2053;
  	 	int hash3 = (121 * this->sparseMatrix[k] + 169)%1031;

  	 	if(hash1min > hash1)
		{
			hash1min = hash1;
			minnet1 = this->sparseMatrix[k];
		}

		if(hash2min > hash2)
		{
			hash2min = hash2;
			minnet2 = this->sparseMatrix[k];
		}

		if(hash3min > hash3)
		{
			hash3min = hash3;
			minnet3 = this->sparseMatrix[k];
		}
  	 }

  	 int part = (196 * minnet1 + 289 * minnet2 + 361 * minnet3 + 529)%partitionCount;
     
     if(part%2 == 1 && rand()%4 == 0)
     {
       part--;
     }
     
  	 partVec[i] = part;
  	 sizeArray[part] += 1;
  	 //calculateCuts3(partitionCount, i);

  	 currVertexCount++;
      
	#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
	#endif

  }
  
  std::cout << std::endl;
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << "-1" << std::endl;
  std::cout << "******PART SIZES*******" << std::endl;
  
  //  for(int i = 0; i < partitionCount; i++){
  //  std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  //}
  
  delete[] sizeArray;
}

void Partitioner::vertexOutput(int algorithm, int seed)
{
  string textName;
  if(algorithm == 0)
    textName = "RANDOMvertex.txt";
  else if(algorithm == 1)
    textName = "P2Nvertex.txt";
  else if(algorithm == 2)
    textName = "N2Pvertex.txt";
  else if(algorithm == 3)
    textName = "N2P_Kvertex.txt";
  else if(algorithm == 4)
    textName = "BFvertex.txt";
  else if(algorithm == 5)
    textName = "BF2vertex.txt";
  else if(algorithm == 6)
    textName = "BF3vertex.txt";
  else if(algorithm == 7)
    textName = "BF4MULTIvertex.txt";
  std::ofstream outfile;
  outfile.open(textName);  
  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
    {
      readOrder.push_back(i);
    }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());

  outfile << "%OUTPUT FOR RUN ID " << seed << " MATRIX ID " << this->fileName << "\n";
  int cur_vertex = 0;
  for (int i : readOrder)
    {
      outfile << cur_vertex++ << "," << this->partVec[i] << "," << this->cutArray[i] << "\n";
    }
  outfile.close();
}

int Partitioner::calculateCuts(int partitionCount)
{
  int cuts = 0;
  std::vector<std::vector<int>> netsInParts(partitionCount);
  for(int i = 0; i < this->vertexCount; i++)
    {
      int part = this->partVec[i];
      for(int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
	{
	  if(std::find(netsInParts[part].begin(), netsInParts[part].end(), this->sparseMatrix[k]) == netsInParts[part].end())
	    {
	      cuts++;
	      netsInParts[part].push_back(this->sparseMatrix[k]);
	    }
	}
    }
  return cuts;
}

//Private methods

int Partitioner::p2nConnectivity(int partitionID, int vertex, const std::vector<std::vector<int>>& partitionToNet)
{
  int connectivityCount = 0;		
  for(int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++)
    {
      if(std::find(partitionToNet[partitionID].begin(), partitionToNet[partitionID].end(), this->sparseMatrix[k]) != partitionToNet[partitionID].end())
	connectivityCount += 1;
    }
  return connectivityCount;
}


///Manual connectivity
/*
  int Partitioner::p2nConnectivity(int partitionID, int vertex, const std::vector<std::vector<int>>& partitionToNet){
  
  int connectivityCount = 1;
  
  for(int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++){//For all edges vertex connected
  int net = this->sparseMatrix[k];
  int vec_size = partitionToNet[partitionID].size();
  for(int n = 0; n < vec_size; n++){ //Check if partition have that net
  -	 if(partitionToNet[partitionID][n] == net)
  connectivityCount++;
  }
  }
  return connectivityCount;
  
*/

int Partitioner::n2pIndex(int vertex, int partitionCount, double capacityConstraint, int* sizeArray, int* indexArray, bool* markerArrayLocal, const std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker)
{   
  for(int i = 0; i < partitionCount; i++)
    {
      indexArray[i] = -1;
      markerArrayLocal[i] = false;
    }
  std::vector<int> encounterArray;
  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
    {
      int edge = this->reverse_sparseMatrix[k];
      int edgeIndex = tracker[edge];
      for (int i = 0; i < netToPartition[edgeIndex]->size(); i++)
	{
	  int part = netToPartition[edgeIndex]->at(i);
	  if(part == -1)
	    continue;
	  if (markerArrayLocal[part])
	    { 
 	      encounterArray[indexArray[part]] += 1;
	    }		    
	  else
	    {
	      encounterArray.push_back(1);
	      indexArray[part] = encounterArray.size() - 1;
	      markerArrayLocal[part] = true;
	    }
	}
    }

  double maxScore = -1.0;
  int maxIndex = -1;
  for (int i = 0; i < partitionCount; i++)
    {
      int connectivity;      
      if (indexArray[i] == -1)
	connectivity = 0;
      else
	connectivity = encounterArray[indexArray[i]];
      double partOverCapacity = sizeArray[i] / capacityConstraint;     
      double penalty = 1 - partOverCapacity;
      double score = penalty * connectivity;
      if (score > maxScore)
	{
	  maxScore = score;
	  maxIndex = i;
	}      
      else if (score == maxScore)
	{
	  
	  if (sizeArray[i] < sizeArray[maxIndex])
	    {
	      maxIndex = i;
	    } 	      
	}
      
    }
  for(int i = 0; i < partitionCount; i++)
    {
      indexArray[i] = -1;
      markerArrayLocal[i] = false;
    }
  return maxIndex;
}

int Partitioner::n2pIndexCorrected(int vertex, int partitionCount, double capacityConstraint, int* sizeArray, int* indexArray, bool* markerArrayLocal, const std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker)
{   
  double maxScore = -1.0;
  int maxIndex = -1;
  int connectivity;
  int degree = this->reverse_sparseMatrixIndex[vertex + 1] - this->reverse_sparseMatrixIndex[vertex];
  if(degree <= 0)
    {
      maxIndex = partitionCount - 1;
      for (int i = partitionCount - 2; i >= 0; i--)
	{
	  if(sizeArray[i] + 1 < capacityConstraint && sizeArray[i] < sizeArray[maxIndex])
	    maxIndex = i;
	}      
    }

  else
    {
      std::vector<int> encounterArray;
      std::vector<int> partArray;
      for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];
	  int edgeIndex = tracker[edge];
	  for (int i = 0; i < netToPartition[edgeIndex]->size(); i++)
	    {
	      int part = netToPartition[edgeIndex]->at(i);
	      if(part == -1)
		continue;
	      if (markerArrayLocal[part])
		{ 
		  encounterArray[indexArray[part]] += 1;
		}		    
	      else
		{
		  encounterArray.push_back(1);
		  indexArray[part] = encounterArray.size() - 1;
		  markerArrayLocal[part] = true;
		  partArray.push_back(part);
		}
	    }
	}
      if(partArray.size() == 0)
	{
	  maxIndex = partitionCount - 1;
	  for (int i = partitionCount - 2; i >= 0; i--)
	    {
	      if(sizeArray[i] + 1 < capacityConstraint && sizeArray[i] < sizeArray[maxIndex])
		maxIndex = i;
	    }      
	}
      else
	{
	  for (int part : partArray)
	    {
	      int size = sizeArray[part];
	      connectivity = encounterArray[indexArray[part]];
	      double partOverCapacity = size / capacityConstraint;
	      double penalty = 1 - partOverCapacity;
	      double score = penalty * connectivity;
	      if (score > maxScore)
		{
		  maxScore = score;
		  maxIndex = part;
		}      
	      else if (score == maxScore)
		{
		  
		  if (size < sizeArray[maxIndex])
		    {
		      maxIndex = part;
		    } 	      
		}
	    }
	}
    }
  
  return maxIndex;
}

int Partitioner::n2pIndexAndScoreCorrected(int vertex, int partitionCount, double capacityConstraint, int* sizeArray, int* indexArray, bool* markerArrayLocal, const std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker, double& scoreToReturn)
{   
  double maxScore = -1.0;
  int maxIndex = -1;
  int connectivity;      
  int degree = this->reverse_sparseMatrixIndex[vertex + 1] - this->reverse_sparseMatrixIndex[vertex];
  if(degree <= 0)
    {
      maxIndex = partitionCount - 1;
      maxScore = 1.0;
      for (int i = partitionCount - 2; i >= 0; i--)
	{
	  if(sizeArray[i] + 1 < capacityConstraint && sizeArray[i] < sizeArray[maxIndex])
	    maxIndex = i;
	}      
    }

  else
    {
      std::vector<int> encounterArray;
      std::vector<int> partArray;
      for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
	{
	  int edge = this->reverse_sparseMatrix[k];
	  int edgeIndex = tracker[edge];
	  for (int i = 0; i < netToPartition[edgeIndex]->size(); i++)
	    {
	      int part = netToPartition[edgeIndex]->at(i);
	      if(part == -1)
		continue;
	      if (markerArrayLocal[part])
		{ 
		  encounterArray[indexArray[part]] += 1;
		}		    
	      else
		{
		  encounterArray.push_back(1);
		  indexArray[part] = encounterArray.size() - 1;
		  markerArrayLocal[part] = true;
		  partArray.push_back(part);
		}
	    }
	}
      if(partArray.size() == 0)
	{
	  maxIndex = partitionCount - 1;
	  maxScore = 0.0;
	  for (int i = partitionCount - 2; i >= 0; i--)
	    {
	      if(sizeArray[i] + 1 < capacityConstraint && sizeArray[i] < sizeArray[maxIndex])
		maxIndex = i;
	    }      
	}
      
      else
	{
	  for (int part : partArray)
	    {
	      connectivity = encounterArray[indexArray[part]];
	      double partOverCapacity = sizeArray[part] / capacityConstraint;
	      double penalty = 1 - partOverCapacity;
	      double score = penalty * connectivity;
	      if (score > maxScore)
		{
		  maxScore = score;
		  maxIndex = part;
		}      
	      else if (score == maxScore)
		{
		  
		  if (sizeArray[part] < sizeArray[maxIndex])
		    {
		      maxIndex = part;
		    } 	      
		}
	    }
	}
    }
  
  for(int i = 0; i < partitionCount; i++)
    {
      indexArray[i] = -1;
      markerArrayLocal[i] = false;
    }
  
  scoreToReturn = maxScore;
  return maxIndex;
}

int Partitioner::n2pIndexAndScore(int vertex, int partitionCount, double capacityConstraint, int* sizeArray, int* indexArray, bool* markerArrayLocal, const std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker, double& scoreToReturn)
{   
  for(int i = 0; i < partitionCount; i++)
    {
      indexArray[i] = -1;
      markerArrayLocal[i] = false;
    }
  std::vector<int> encounterArray;
  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
    {
      int edge = this->reverse_sparseMatrix[k];
      int edgeIndex = tracker[edge];
      for (int i = 0; i < netToPartition[edgeIndex]->size(); i++)
	{
	  int part = netToPartition[edgeIndex]->at(i);
	  if(part == -1)
	    continue;
	  if (markerArrayLocal[part])
	    { 
 	      encounterArray[indexArray[part]] += 1;
	    }		    
	  else
	    {
	      encounterArray.push_back(1);
	      indexArray[part] = encounterArray.size() - 1;
	      markerArrayLocal[part] = true;
	    }
	}
    }
  
  double maxScore = -1.0;
  int maxIndex = -1;
  for (int i = 0; i < partitionCount; i++)
    {
      int connectivity;      
      if (indexArray[i] == -1)
	connectivity = 0;
      else
	connectivity = encounterArray[indexArray[i]];
      double partOverCapacity = sizeArray[i] / capacityConstraint;     
      double penalty = 1 - partOverCapacity;
      double score = penalty * connectivity;
      if (score > maxScore)
	{
	  maxScore = score;
	  maxIndex = i;
	}      
      else if (score == maxScore)
	{
	  
	  if (sizeArray[i] < sizeArray[maxIndex])
	    {
	      maxIndex = i;
	    } 	      
	}
      
    }

  for(int i = 0; i < partitionCount; i++)
    {
      indexArray[i] = -1;
      markerArrayLocal[i] = false;
    }

  scoreToReturn = maxScore;
  return maxIndex;
}

int Partitioner::BFConnectivity(Bloom<int, int>* bloomFilter, int partitionID, int vertex)
{
  int connectivityCount = 0;
  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
    {
      int edge = this->reverse_sparseMatrix[k];
      if (bloomFilter->contains(edge, partitionID)){
	connectivityCount++;
      }
    }
  
  return connectivityCount;
}

int Partitioner::BFConnectivity2(BloomFilter* bf, int partitionID, int vertex)
{
  int connectivityCount = 0;
  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
    {
      int val = this->reverse_sparseMatrix[k];
      if (bf[partitionID].query(val))
	connectivityCount++;
    }
  
  return connectivityCount;
}

int Partitioner::BFConnectivity3(BloomFilter_OT* bloomFilter, int partitionID, int vertex)
{
  int connectivityCount = 0;
  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
    {
      int edge = this->reverse_sparseMatrix[k];
      if (bloomFilter->query(edge, partitionID))
	connectivityCount++;
    }
  
  return connectivityCount;
}  

int Partitioner::BFConnectivityMult(mlbf* bloomFilter, int partitionID, int vertex)
{
  int connectivityCount = 0;
  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
    {
      int edge = this->reverse_sparseMatrix[k];
      if (bloomFilter->query(edge, partitionID))
	connectivityCount++;
    }
  
  return connectivityCount;
}

void Partitioner::BFConnectivityMult2(mlbf_2* bloomFilter, bool* existences, int* connectivities, int no_leaf_blocks, int vertex)
{
  for (int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++){

    int edge = reverse_sparseMatrix[k];
    
    for(int i = 0; i < no_leaf_blocks; i++){
      existences[i] = false;
    }
    
    bloomFilter->query(edge, existences);
    
    for(int i = 0; i < no_leaf_blocks; i++){
      if(existences[i])
	connectivities[i] += 1;
    }
    
  }  
}
int Partitioner::calculateCuts2(int partitionCount)
{
  int cuts = 0;
  bool* arr = new bool[partitionCount];
  
  for(int b = 0; b < partitionCount; b++){
    arr[b] = false;
  }
  
  for(int i = 0; i < this->edgeCount; i++)
    {
      int net = i;
      int cut = 0;
      for(int k = this->sparseMatrixIndex[net]; k < this->sparseMatrixIndex[net + 1]; k++)
	{
	  int vertex = sparseMatrix[k];
	  int part = this->partVec[vertex];
	  
	  arr[part] = true;	  
	}
      for(int b = 0; b < partitionCount; b++){
	if(arr[b])
	  {
	    cut++;
	  }
      }
      if(cut != 0)
	cuts += cut-1;
      for(int b = 0; b < partitionCount; b++){
	arr[b] = false;
      }
    }
  
  //delete[] arr;
  return cuts;
}

int Partitioner::calculateCuts3(int partitionCount, int vertex, const std::vector<int>& tracker, const std::vector<std::vector<int>*>& netToPartition, int pti)
{
  int cuts = 0;

  for(int k = this->reverse_sparseMatrixIndex[vertex]; k < this->reverse_sparseMatrixIndex[vertex + 1]; k++)
    {
      int net = reverse_sparseMatrix[k];
      int cut = netToPartition[tracker[net]]->size() - 1;
      if(cut != 0)
	cuts += (cut - 1);      
    }  
  this->cutArray[vertex] = cuts;
  return cuts;
}

int Partitioner::calculateCuts4(int partitionCount, const std::vector<std::vector<int>*>& netToPartition)
{
  int cuts = 0;
  for(int i = 0; i < netToPartition.size(); i++)
    {
      int cut = netToPartition[i]->size();
      if (cut != 0)
	{
	  cuts += netToPartition[i]->size() - 1;
	}
    }
  return cuts;
}
