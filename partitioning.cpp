#include "partitioning.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>

//Public methods
Algorithms::Algorithms(std::string fileName, int partitionCount)
{
	//Read matrix
	std::ifstream fin(fileName);
	while(fin.peek() == '%') fin.ignore(2048, '\n');
	fin >> this->edgeCount >> this->vertexCount >> this->nonzeroCount;
	std::cout << "Edge count: " << this->edgeCount << " Vertex count: " << this->vertexCount << " Nonzero count: " << this->nonzeroCount << std::endl;

	//Init variables
	this->partitionCount = partitionCount;
	this->capacityConstraint = this->vertexCount / this->partitionCount;
	int capacity = this->vertexCount / this->partitionCount;

	//Generate random read order
	for (int i = 0; i < this->vertexCount; i++)
	{
		this->readOrder.push_back(i);
	}
	std::random_shuffle(this->readOrder.begin(), this->readOrder.end());

	//Init partition matrix
	this->partitions = new int*[this->partitionCount];
	for (int i = 0; i < this->partitionCount; i++)
	{
		this->partitions[i] = new int[capacity];
	}

	//Init partitionToNet structure
	this->partitionToNet = new bool*[this->partitionCount];
	for (int i = 0; i < this->partitionCount; i++)
	{
		this->partitionToNet[i] = new bool[this->edgeCount]; //Kapasite fazla 
	}

	//Init netToPartition structure
	for (int i = 0; i < this->edgeCount; i++)
	{
		std::vector<int> edge;
		this->netToPartition.push_back(edge);
	}
	std::cout << "Initialization: DONE!" << std::endl;

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
			this->sparseMatrixIndex[vIndex] = col;
			currentColumn = col;
			vIndex++;
		}
	}
	this->sparseMatrix[this->nonzeroCount] = this->sparseMatrix[this->nonzeroCount - 1] + 1;
	this->sparseMatrixIndex[this->vertexCount] = this->nonzeroCount + 1;
	std::cout << "Matrix integration: DONE!" << std::endl;	
}

Algorithms::Algorithms(std::string fileName, int partitionCount, int byteSize, int hashCount)
{
	//Read matrix
	std::ifstream fin(fileName);
	while (fin.peek() == '%') fin.ignore(2048, '\n');
	fin >> this->edgeCount >> this->vertexCount >> this->nonzeroCount;
	std::cout << "Edge count: " << this->edgeCount << " Vertex count: " << this->vertexCount << " Nonzero count: " << this->nonzeroCount << std::endl;

	//Init variables
	this->partitionCount = partitionCount;
	this->capacityConstraint = this->vertexCount / this->partitionCount;
	int capacity = this->vertexCount / this->partitionCount;

	//Generate random read order
	for (int i = 0; i < this->vertexCount; i++)
	{
		this->readOrder.push_back(i);
	}
	std::random_shuffle(this->readOrder.begin(), this->readOrder.end());

	//Init partition matrix
	this->partitions = new int*[this->partitionCount];
	for (int i = 0; i < this->partitionCount; i++)
	{
		this->partitions[i] = new int[capacity];
	}

	//Init bloom filter
	this->bloomFilter = new Bloom<int, int>(byteSize, hashCount);
	std::cout << "Initialization: DONE!" << std::endl;

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
	for (int i = 0; i < this->partitionCount; i++)
	{
		delete[] this->partitions[i];
		delete[] this->partitionToNet[i];
	}
	delete[] this->partitions;
	delete[] this->partitionToNet;

	delete[] this->sparseMatrix;
	delete[] this->sparseMatrixIndex;
}

void Algorithms::LDGp2n()
{
	int* sizeArray = new int[this->partitionCount];
	for (int i = 0; i < this->partitionCount; i++)
	{
		sizeArray[i] = 0;
	}	
	
	for (int i : this->readOrder)
	{
		double maxScore = -1.0;
		int maxIndex = -1;
		for (int j = 0; j < this->partitionCount; j++)
		{
			int connectivity = this->p2nConnectivity(j, i);
			double partToCapacity = sizeArray[j] / this->capacityConstraint;
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
		partitions[maxIndex][sizeArray[maxIndex]] = i;
		sizeArray[maxIndex] += 1;
		for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
		{
			this->partitionToNet[maxIndex][this->sparseMatrix[k]] = true;
		}
	}
	
	delete[] sizeArray;
}

void Algorithms::LDGn2p()
{
	int* sizeArray = new int[this->partitionCount];
	for (int i = 0; i < this->partitionCount; i++)
	{
		sizeArray[i] = 0;
	}

	for (int i : this->readOrder)
	{
		int maxIndex = this->n2pIndex(i, sizeArray);
		partitions[maxIndex][sizeArray[maxIndex]] = i;
		sizeArray[maxIndex] += 1;
		for (int k = this->sparseMatrixIndex[i]; k < this->sparseMatrixIndex[i + 1]; k++)
		{
			int edge = this->sparseMatrix[k];
			if (std::find(netToPartition[edge].begin(), netToPartition[edge].end(), maxIndex) == netToPartition[edge].end())
			{
				netToPartition[edge].push_back(maxIndex);
			}			
		}
	}
	
	delete[] sizeArray;
}

void Algorithms::LDGBF()
{
	int* sizeArray = new int[this->partitionCount];
	for (int i = 0; i < this->partitionCount; i++)
	{
		sizeArray[i] = 0;
	}

	double maxScore = -1.0;
	int maxIndex = -1;
	for (int i : this->readOrder)
	{
		for (int j = 0; j < this->partitionCount; j++)
		{
			int connectivity = this->BFConnectivity(j, i);
			double partToCapacity = sizeArray[j] / this->capacityConstraint;
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

int Algorithms::calculateCuts()
{
	int cuts = 0;
	for (std::vector<int> edge : this->netToPartition)
	{
		int vertexCount = 0;
		for (int i : edge)
		{
			vertexCount++;
		}		
		cuts += (vertexCount * (vertexCount - 1)) / 2;
	}

	return cuts;
}

//Private methods
int Algorithms::p2nConnectivity(int partitionID, int vertex)
{
	int connectivityCount = 0;		
	for(int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++)
	{
		int edge = this->sparseMatrix[k];
		if (this->partitionToNet[partitionID][edge])
			connectivityCount++;
	}

	return connectivityCount;
}

int Algorithms::n2pIndex(int vertex, int* sizeArray)
{
	int* connectivities = new int[this->partitionCount]();
	for (int k = this->sparseMatrixIndex[vertex]; k < this->sparseMatrixIndex[vertex + 1]; k++)
	{
		int edge = this->sparseMatrix[k];
		for (int i = 0; i < this->netToPartition[edge].size(); i++)
		{
			int part = this->netToPartition[edge][i];
			connectivities[part] += 1;
		}
	}
	double maxScore = -1.0;
	int maxIndex = -1;

	for (int i = 0; i < this->partitionCount; i++)
	{
		int connectivity = connectivities[i];
		double partToCapacity = sizeArray[i] / this->capacityConstraint;
		double penalty = 1 - partToCapacity;
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
	delete[] connectivities;
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










