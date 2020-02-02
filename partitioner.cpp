#include "partitioner.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <chrono>
//
#include <iomanip>
#include <stdio.h>

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
  this->scoreArray = new double[this->vertexCount];
}

  
  void Partitioner::read_binary_graph(std::string fileName){
  const char* fname = fileName.c_str();
  FILE* bp;
  bp = fopen(fname, "r");
  
  int* nov = new int;
  int* nnet = new int;
  int* nnz = new int;
  
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

  this->partVec = new int[this->vertexCount];
  this->bloomFilter = nullptr;
  
#ifdef DEBUG
  std::cout << "First indexes of xadj and adj" << std::endl;
  for(int i = 0; i < 150; i++){
    std::cout << "i: " << i  << " " <<this->sparseMatrixIndex[i] << " " << sparseMatrix[i] << std::endl;
  }
#endif
  
  fclose(bp);
  }

bool sortbyfirst(const pair<int,int> &a, const pair<int,int> &b) 
{ 
  return (a.first < b.first); 
} 

bool sortbysec(const pair<int,int> &a, const pair<int,int> &b) 
{ 
  return (a.second < b.second); 
} 

void Partitioner::read_mtx_and_transform_to_shpbin(std::string fileName){

  std::string mtx_name = fileName + ".mtx";
  std::string bin_name = fileName + ".shpbin";
  const char* bfile = bin_name.c_str();

  int* xadj;
  int* adj;
  
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
    
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbyfirst);
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbysec);
    
    adj = new int[intermediate.size()];
    xadj = new int[no_col+1];
    
    xadj[0] = 0;
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
      
      std::cout << "i: " << i << " ||| xadj: " << xadj[i] << " edge: " << adj[i] << std::endl;
    }

    std::cout << intermediate.size() << " ||| " << nnz+1 << std::endl;
    std::cout << intermediate[intermediate.size()-1].second << " ||| " << no_col << std::endl;
#endif
        
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
    
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbyfirst);
    std::stable_sort(intermediate.begin(), intermediate.end(), sortbysec);

    //need vector::erase and remove replicatons BEFORE continuing to adj and xadj because we need the exact size of purged intermediate to write binary files 
    
    for(int i = 0; i < intermediate.size(); i++){
      if((intermediate[i].first == intermediate[i+1].first) && intermediate[i].second == intermediate[i+1].second)
	intermediate.erase(intermediate.begin()+i);
      //std::cout << "Deduplicating " << i << "/" << intermediate.size() << std::endl;
    }
    
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
    
  }
  
  if(!_general && !_symmetric)
  {
    std::cout << "I believe a problem happened during reading fields of the matrix.." << std::endl;
    exit(1);
  }    
  
  fin.close();
  
  
  int* no_vertex;
  no_vertex = &no_row;
  int* no_nets;
  no_nets = &intermediate[intermediate.size()-1].second;
  int np_nz = intermediate.size();
  int* no_nz;
  no_nz = &np_nz;
  
  const char* fname = bin_name.c_str();
  FILE* bp;
  bp = fopen(fname, "w");
  
#ifdef DEBUG
  std::cout << "No vertex: " << *no_vertex << " No nets: " << *no_nets << " No nz: " << *no_nz << std::endl;
#endif

  fwrite(no_vertex, sizeof(int), 1, bp);
  fwrite(no_nets, sizeof(int), 1, bp);
  fwrite(no_nz, sizeof(int), 1, bp);
  fwrite(xadj, sizeof(int), *no_nets, bp);
  fwrite(adj, sizeof(int), *no_nz, bp);
  
  fclose(bp);

  delete[] xadj;
  delete[] adj;
}


Partitioner::~Partitioner()
{
  delete[] this->partVec;
  delete[] this->scoreArray;
  delete[] this->sparseMatrix;
  delete[] this->sparseMatrixIndex;
  
  if(this->bloomFilter)
    delete this->bloomFilter;
}

void Partitioner::partition(int algorithm, int partitionCount, int slackValue, int seed, double imbal)
{ 
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
    std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
  }
  else if(algorithm == 3)
  {
    auto start = std::chrono::high_resolution_clock::now();
    this->LDGn2p_i(partitionCount, slackValue, seed, imbal);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
  }
  else if(algorithm == 4)
  {
    auto start = std::chrono::high_resolution_clock::now();
    this->LDGBF(partitionCount, slackValue, seed, imbal);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Duration:" << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
  }
  
  std::cout << "Cuts:" << this->calculateCuts(partitionCount) << std::endl;
  this->vertexOutput(algorithm, seed);
  //compute cut and report  
}

//ALGO 0//
void Partitioner::RandomPartition(int partitionCount, int seed)
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
  
  for (int i : readOrder)
  {
    int partition = distribution(generator);
    sizeArray[partition] += 1;
    partVec[i] = partition;
  }
  
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  }

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
    scoreArray[i] = maxScore;
    partVec[i] = maxIndex;
    sizeArray[maxIndex] += 1;
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
  
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  }

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
  //THIS LOOKS LIKE ITS CORRECT BUT ITS NOT//  
  /*for(int i = 0; i < sparseMatrixIndex[this->vertexCount]; i++){
    //std::cout << "i: " << i <<" sparseMatrixIndex[vertexCount]: " << sparseMatrixIndex[this->vertexCount] << " vertexCount: " << this->vertexCount  << " ";
    //std::cout << "spsM[" <<i <<"]: "<< sparseMatrix[i] << std::endl;
  }*/

  /*
  for(int i = 0; i < this->vertexCount; i++){
    std::cout << "sparseMatrixIndex[" << i << "]: " << this->sparseMatrixIndex[i] << std::endl;
  }
  */

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

    for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
      {
	//std::cout << "k: " << k << std::endl;
	//std::cout << "index[i]: " << sparseMatrixIndex[i] << std::endl;
	//std::cout << "index[vertexCount]: " << this->sparseMatrixIndex[this->vertexCount]<< " k: " << k << std::endl;
	int edge = this->sparseMatrix[k];
	//std::cout << "edge: " << edge << std::endl;
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
    
    int maxIndex = this->n2pIndex(i, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker); 
    partVec[i] = maxIndex;
    sizeArray[maxIndex] += 1;
    
    for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
      {
	int edge = this->sparseMatrix[k];      
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
  std::cout << "******PART SIZES*******" << std::endl;
    
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  }
  
  delete[] sizeArray;
  delete[] indexArray;
  delete[] markerArray;
}

void Partitioner::LDGn2p_i(int partitionCount, int slackValue, int seed, double imbal)
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
  std::uniform_int_distribution<int> distribution(0,INITVECSIZE - 1);
  for(int i : readOrder)
  {
    if((imbal*currVertexCount) >= slackValue)
      capacityConstraint = (imbal*currVertexCount) / partitionCount;
    else
      capacityConstraint = ((double)slackValue) / partitionCount;
    for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
    {
      int edge = this->sparseMatrix[k];
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
    int maxIndex = this->n2pIndex(i, partitionCount, capacityConstraint, sizeArray, indexArray, markerArray, netToPartition, tracker);
    partVec[i] = maxIndex;
    sizeArray[maxIndex] += 1;
    
    for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
    {
      int edge = this->sparseMatrix[k];      
      if(std::find (netToPartition[tracker[edge]]->begin(), netToPartition[tracker[edge]]->end(), maxIndex) == netToPartition[tracker[edge]]->end())
      {
        if(netToPartition[tracker[edge]]->size() != INITVECSIZE)
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
  std::cout << "******PART SIZES*******" << std::endl;
    
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size:" << sizeArray[i] << std::endl;
  }
  
  delete[] sizeArray;
  delete[] indexArray;
  delete[] markerArray;
}

void Partitioner::LDGBF(int partitionCount, int slackValue, int seed, double imbal)
{
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
	  int connectivity = this->BFConnectivity(j, i);
	  //std::cout <<"partition " << j <<  " Conn: " << this->BFConnectivity(j, i) << std::endl;
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
      scoreArray[i] = maxScore;
      sizeArray[maxIndex] += 1;
      
      for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
	{
	  int edge = this->sparseMatrix[k];
	  //if (!(this->bloomFilter->contains(edge, maxIndex)))
	    this->bloomFilter->insert(edge, maxIndex);
	}
      
      currVertexCount++;
      
#ifdef WATCH
      std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << std::flush;
      std::cout << std::fixed << std::setprecision(2) << "Progress: " << ((double)currVertexCount/this->vertexCount)*100 << "%" << std::flush;;
#endif

      
      
    }

  std::cout << std::endl;
  
  std::cout << "******PART SIZES*******" << std::endl;
  
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  }
  
  delete[] sizeArray;
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

void Partitioner::vertexOutput(int algorithm, int seed)
{
  string textName;

  if(algorithm == 1)
      textName = "P2Nvertex.txt";
  else if(algorithm == 2)
      textName = "N2Pvertex.txt";
  else if(algorithm == 3)
      textName = "N2P_Kvertex.txt";
  else if(algorithm == 3)
      textName = "BFvertex.txt";
    
  std::ofstream outfile;
  outfile.open(textName, std::ios_base::app);

  std::vector<int> readOrder;
  for (int i = 0; i < this->vertexCount; i++)
  {
    readOrder.push_back(i);
  }
  std::srand(seed);
  std::random_shuffle(readOrder.begin(), readOrder.end());

  for (int i : readOrder)
  {
        outfile << i << "," << partVec[i] << "," << scoreArray[i] << "\n";
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
	 if(partitionToNet[partitionID][n] == net)
	   connectivityCount++;
       }
     }
   return connectivityCount;
 }
  */

int Partitioner::n2pIndex(int vertex, int partitionCount, double capacityConstraint, int* sizeArray, int* indexArray, bool* markerArray, const std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker)
{    
  std::vector<int> encounterArray;
	for (int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++)
	{
	  int edge = this->sparseMatrix[k];
	  int edgeIndex = tracker[edge];
    
	  for (int i = 0; i < netToPartition[edgeIndex]->size(); i++)
	  {
      	      int part = netToPartition[edgeIndex]->at(i);
              if (markerArray[part])
              {       
	        encounterArray[indexArray[part]] += 1;        
              }		    
	      else
	      {
	        encounterArray.push_back(1);
	        indexArray[part] = encounterArray.size() - 1;
	        markerArray[part] = true;
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

    this->scoreArray[vertex] = maxScore;	
    return maxIndex;
}

int Partitioner::BFConnectivity(int partitionID, int vertex)
{
  int connectivityCount = 0;
  for (int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++)
    {
      int edge = this->sparseMatrix[k];
      if ((this->bloomFilter)->query(edge, partitionID))
	connectivityCount++;
	}
  
  return connectivityCount;
}
