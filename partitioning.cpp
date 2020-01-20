#include "partitioning.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>


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

std::default_random_engine generator;
std::uniform_int_distribution<int> distribution(0,16);

//Public methods
Algorithms::Algorithms(std::string fileName) {
  //Read matrix
  std::ifstream fin(fileName);
  while(fin.peek() == '%') fin.ignore(2048, '\n');
  fin >> this->edgeCount >> this->vertexCount >> this->nonzeroCount;
  std::cout << "Row count: " << this->edgeCount << " Column count: " << this->vertexCount << " Non-zero count: " << this->nonzeroCount << std::endl;

  //Init partition matrix
  this->partVec = new int[this->vertexCount];  
  this->bloomFilter = nullptr;
  
  //Init sparse matrix representation
  this->sparseMatrix = new int[this->nonzeroCount + 1];
  this->sparseMatrixIndex = new int[this->vertexCount + 1];
  sparseMatrix[0] = 0;
  int vIndex = 0, row, col, currentColumn = -1;
  double value;	
  for(int i = 0; i < this->nonzeroCount; i++)
  {
      fin >> row >> col >> value;
      this->sparseMatrix[i] = row;
      if (col != currentColumn)
    	{
    	  this->sparseMatrixIndex[vIndex] = i;
    	  currentColumn = col;
    	  vIndex++;
    	}
  }
  
  /*
  for(int i = 0; i < this->vertexCount + 1; i++)
  {
    std::cout << this->sparseMatrixIndex[i] << std::endl;
  }
  */

  this->sparseMatrix[this->nonzeroCount] = this->sparseMatrix[this->nonzeroCount - 1] + 1;
  this->sparseMatrixIndex[this->vertexCount] = this->nonzeroCount + 1;
  std::cout << "Matrix integration: DONE!" << std::endl;	
}

Algorithms::Algorithms(std::string fileName, int byteSize, int hashCount)
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
	int vIndex = 0, row, col, currentColumn = -1;
	double value;
	for (int i = 0; i < this->nonzeroCount; i++)
	{
		fin >> row >> col >> value;
		this->sparseMatrix[i] = row;
		if (col != currentColumn)
		{
			this->sparseMatrixIndex[vIndex] = col;
			currentColumn = col;
			vIndex++;
		}
	}
	this->sparseMatrix[this->nonzeroCount] = this->sparseMatrix[this->nonzeroCount - 1] + 1;
	this->sparseMatrixIndex[this->vertexCount] = this->nonzeroCount + 1;
	std::cout << "Matrix integration: DONE!" << std::endl;
}

Algorithms::~Algorithms()
{
  delete[] this->partVec;
  delete[] this->sparseMatrix;
  delete[] this->sparseMatrixIndex;
  
  if(this->bloomFilter)
    delete this->bloomFilter;
}

void Algorithms::partition(int algorithm, int partitionCount, double imbal)
{ 
  //Partition
  if(algorithm == 1)
  {
    LDGp2n(partitionCount, imbal);
  }
  else if(algorithm == 2)
  {
    LDGn2p(partitionCount, imbal);  
  }
  else if(algorithm == 3)
  {
    LDGBF(partitionCount, imbal);
  }

  //compute cut and report  
}


//ALGO 1//
void Algorithms::LDGp2n(int partitionCount, double imbal)
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
  
  double capacityConstraint = (imbal*this->vertexCount) / partitionCount;
  for (int i : readOrder)
    {
      std::cout << "Partitioning vertex: " << ++ctr << "/" << this->vertexCount << std::endl;
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
      std::cout << "Vertex " << i << " assigned to partition " << maxIndex << std::endl;
      partVec[i] = maxIndex;
      sizeArray[maxIndex] += 1;
      int maxIndexSize = partitionToNet[maxIndex].size() - 1;


      for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
	{
	  if(std::find (partitionToNet[maxIndex].begin(), partitionToNet[maxIndex].end(), this->sparseMatrix[k]) == partitionToNet[maxIndex].end())
       partitionToNet[maxIndex].push_back(this->sparseMatrix[k]);
	}
    }
  
  std::cout << "******PART SIZES*******" << std::endl;
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  }

  delete sizeArray;
}

//ALGO 2//
void Algorithms::LDGn2p(int partitionCount, double imbal)
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
  double capacityConstraint = (imbal*this->vertexCount) / partitionCount;
  
  for (int i : readOrder) {
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
      markerArray[i] = 0;
    }
  }
  
  for(int i = 0; i < partitionCount; i++){
    std::cout << "part " << i << " size: " << sizeArray[i] << std::endl;
  }
  
  delete[] sizeArray;
  delete[] indexArray;
  delete[] markerArray;
}

void Algorithms::LDGBF(int partitionCount, double imbal)
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

	double maxScore = -1.0;
	int maxIndex = -1;
  double capacityConstraint = (imbal*this->vertexCount) / partitionCount;
	for (int i : readOrder)
	{
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

	}
	
	delete[] sizeArray;
}

void Algorithms::LDGMultiBF()
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

/*int Algorithms::calculateCuts()
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
 
int Algorithms::p2nConnectivity(int partitionID, int vertex, const std::vector<std::vector<int>>& partitionToNet)
{
  int connectivityCount = 0;		
  for(int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++)
    {
      if(std::find(partitionToNet[partitionID].begin(), partitionToNet[partitionID].end(), this->sparseMatrix[k]) != partitionToNet[partitionID].end())
	connectivityCount += 1;
	}
  
  if(connectivityCount > 3)
    std::cout << "cc: " << connectivityCount << std::endl;
  return connectivityCount;
}
 

  ///Manual connectivity
  /*
 int Algorithms::p2nConnectivity(int partitionID, int vertex, const std::vector<std::vector<int>>& partitionToNet){
   
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

int Algorithms::n2pIndex(int vertex, int partitionCount, double capacityConstraint, int* sizeArray, int* indexArray, bool* markerArray, const std::vector<std::vector<int>*>& netToPartition, const std::vector<int>& tracker)
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

int Algorithms::BFConnectivity(int partitionID, int vertex)
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