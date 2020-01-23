#include "partitioner.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

#define DEBUG


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
  std::ifstream fin(fileName);
  std::string comment;
  std::getline(fin, comment);
  
  
  while(fin.peek() == '%')
  {    
    fin.ignore(2048, '\n');
  }  
  
  //Getting net, pin, non-zero counts
  fin >> this->edgeCount >> this->vertexCount >> this->nonzeroCount;
  
  std::cout << "Row count: " << this->edgeCount << " Column count: " << this->vertexCount << " Non-zero count: " << this->nonzeroCount << std::endl;
  
  //Init partition matrix
  this->partVec = new int[this->vertexCount];  
  this->bloomFilter = nullptr;
  
  //Init sparse matrix representation
  this->sparseMatrixIndex = new int[this->vertexCount + 1];
  sparseMatrixIndex[0] = 0;
  if(comment.find("symmetric") == std::string::npos)
    this->sparseMatrix = new int[this->nonzeroCount + 1];
  else
    this->sparseMatrix = new int[this->nonzeroCount*2 + 1];
  
  fin.close();
  
  this->read_graph(fileName);
  
}
  

void Partitioner::read_graph(std::string fileName){
    
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
  std::ifstream fin(fileName);
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
  std::getline(fin, comment);//We already acquired that values
  
  
  int vIndex = 0, row, col, currentColumn = -1;
  float val1;
  float val2;
  int val3;
  if(_general){	
    
    if(_real){
      
      for(int i = 0; i < this->nonzeroCount + 1; i++)
	    {      
	      fin >> row >> col >> val1;
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
	        this->sparseMatrixIndex[vIndex] = i;
	        currentColumn = col;
	        vIndex++;
	      }
	    }
    }
    
    if(_integer){
      
      for(int i = 0; i < this->nonzeroCount + 1; i++)
	    {      
	      fin >> row >> col >> val2;
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
	        this->sparseMatrixIndex[vIndex] = i;
	        currentColumn = col;
	        vIndex++;
	      }
	    }  
    }
    
    
    if(_pattern){
    
      for(int i = 0; i < this->nonzeroCount; i++)
	    {            
	      fin >> row >> col;
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
	        this->sparseMatrixIndex[vIndex] = i;
	        currentColumn = col;
	        vIndex++;
	      }
	    }
    }
    
    if(_complex){      
      
      for(int i = 0; i < this->nonzeroCount + 1; i++)
	    {      
	      fin >> row >> col >> val1 >> val2;
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
	        this->sparseMatrixIndex[vIndex] = i;
	        currentColumn = col;
	        vIndex++;
	      }
	    }
    }
  }
  
  if(_symmetric){
    
    if(_real){
     
      for(int i = 0; i < this->nonzeroCount + 1; i++)
	    {      
	      fin >> row >> col >> val1;
	      //ALSO NEED TO ADD (j,i) as well as (i,j)
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
		      this->sparseMatrixIndex[vIndex] = i;
		      currentColumn = col;
		      vIndex++;
	      }
	    }
    }
    
    if(_integer){
      
      for(int i = 0; i < this->nonzeroCount + 1; i++)
	    {      
	      fin >> row >> col >> val2;
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
	        this->sparseMatrixIndex[vIndex] = i;
	        currentColumn = col;
	        vIndex++;
	      }
	    }  
    }
    
    if(_pattern){
      for(int i = 0; i < this->nonzeroCount; i++)
	    {            
	      fin >> row >> col;
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
		      this->sparseMatrixIndex[vIndex] = i;
		      currentColumn = col;
		      vIndex++;
	      }
	    }
    }
    
    if(_complex){
      
      for(int i = 0; i < this->nonzeroCount + 1; i++)
	    {      
	      fin >> row >> col >> val1 >> val2;
	      this->sparseMatrix[i] = row - 1;
	      if (col != currentColumn)
	      {
	        this->sparseMatrixIndex[vIndex] = i;
	        currentColumn = col;
	        vIndex++;
	      }
	    }
    }
  }  
  
  if(!_general && !_symmetric)
  {
    std::cout << "I believe a problem happened during reading fields of the matrix" << std::endl;
    exit(1);
  }
    
  
  this->sparseMatrix[this->nonzeroCount] = this->sparseMatrix[this->nonzeroCount - 1] + 1;
  this->sparseMatrixIndex[this->vertexCount] = this->nonzeroCount + 1;
  std::cout << "Matrix integration: DONE!" << std::endl;	
}


Partitioner::Partitioner(std::string fileName, int byteSize, int hashCount)
{
  //Read matrix
	std::ifstream fin(fileName);
	while (fin.peek() == '%') fin.ignore(2048, '\n');
	fin >> this->edgeCount >> this->vertexCount >> this->nonzeroCount;
	std::cout << "Edge count: " << this->edgeCount << " Vertex count: " << this->vertexCount << " Nonzero count: " << this->nonzeroCount << std::endl;

	//Init partition matrix
	this->partVec = new int[this->vertexCount];

	//Init bloom filter
	this->bloomFilter = new Bloom<int, int>(byteSize, hashCount);

	//Init sparse matrix representation
	this->sparseMatrix = new int[this->nonzeroCount + 1];
	this->sparseMatrixIndex = new int[this->vertexCount + 1];
	sparseMatrix[0] = 0;
	
	fin.close();
	this->read_graph(fileName);

}

Partitioner::~Partitioner()
{
  delete[] this->partVec;
  delete[] this->sparseMatrix;
  delete[] this->sparseMatrixIndex;
  
  if(this->bloomFilter)
    delete this->bloomFilter;
}

void Partitioner::partition(int algorithm, int partitionCount, int slackValue, double imbal)
{ 
  //Partition
  if(algorithm == 1)
  {
    this->LDGp2n(partitionCount, slackValue, imbal);
  }
  else if(algorithm == 2)
  {
    this->LDGn2p(partitionCount, slackValue, imbal);  
  }
  else if(algorithm == 3)
  {
    this->LDGBF(partitionCount, slackValue, imbal);
  }

  //compute cut and report  
}


//ALGO 1//
void Partitioner::LDGp2n(int partitionCount, int slackValue, double imbal)
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
      capacityConstraint = slackValue / partitionCount;
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
    int maxIndexSize = partitionToNet[maxIndex].size() - 1;
    
    for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
	  {
	    if(std::find (partitionToNet[maxIndex].begin(), partitionToNet[maxIndex].end(), this->sparseMatrix[k]) == partitionToNet[maxIndex].end())
         partitionToNet[maxIndex].push_back(this->sparseMatrix[k]);
	  }
    currVertexCount++;
  }
  
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;
  std::cout << "******PART SIZES*******" << std::endl;
  
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  }

  delete[] sizeArray;
}

//ALGO 2//
void Partitioner::LDGn2p(int partitionCount, int slackValue, double imbal)
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
  std::random_shuffle(readOrder.begin(), readOrder.end());
  
  std::vector<std::vector<int>*> netToPartition;  
  std::vector<int> tracker(10000, -1);
  double capacityConstraint;
  int currVertexCount = 0;
  for (int i : readOrder) {    
    if((imbal*currVertexCount) >= slackValue)
      capacityConstraint = (imbal*currVertexCount) / partitionCount;
    else
      capacityConstraint = slackValue / partitionCount;
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
        netToPartition[tracker[edge]]->push_back(maxIndex);      
    }
 
    for (int i = 0; i < partitionCount; i++) {
      indexArray[i] = -1;
      markerArray[i] = false;
    }
    
    currVertexCount++;
  }
  
  std::cout << "MAX ALLOWED PART COUNT: " << MAXPARTNO << " - PART COUNT: " << partitionCount <<  std::endl;
  std::cout << "MAX ALLOWED IMBALANCE RATIO: " << MAXIMBAL << " - IMBALANCE RATIO: " << imbal << std::endl;
  std::cout << "******PART SIZES*******" << std::endl;
    
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  }
  
  delete[] sizeArray;
  delete[] indexArray;
  delete[] markerArray;
}

void Partitioner::LDGBF(int partitionCount, int slackValue, double imbal)
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
  std::random_shuffle(readOrder.begin(), readOrder.end());
	
  double capacityConstraint;
  int currVertexCount = 0;
	for (int i : readOrder)
	{
    if((imbal*currVertexCount) >= slackValue)
      capacityConstraint = (imbal*currVertexCount) / partitionCount;
    else
      capacityConstraint = slackValue / partitionCount;
    double maxScore = -1.0;
	  int maxIndex = -1;
		for (int j = 0; j < partitionCount; j++)
		{
			int connectivity = this->BFConnectivity(j, i);
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
		for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
		{
			int edge = this->sparseMatrix[k];
			if (!(this->bloomFilter->contains(k, maxIndex)))
				this->bloomFilter->insert(k, maxIndex);
		}
    currVertexCount++;
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

/*int Partitioner::calculateCuts()
{
	int cuts = 0;
	for (std::vector<int> edge : this->netToPartition)
	{
		int partCount = 0;
		for (int i : edge)
		{
			partCount++;
		}		
		cuts += partCount - 1;
	}
	return cuts;
}*/

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
	
	return maxIndex;
}

int Partitioner::BFConnectivity(int partitionID, int vertex)
{
	int connectivityCount = 0;
	for (int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++)
	{
		int edge = this->sparseMatrix[k];
		if ((this->bloomFilter)->contains(edge, partitionID))
			connectivityCount++;
	}

	return connectivityCount;
}
