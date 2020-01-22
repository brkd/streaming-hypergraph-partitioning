#include "partitioning.cpp"
#include <chrono>

int main(int argc, char** argv) {
  std::string fileName;
  int partitionCount, randomizationCount;
  int algorithm = atoi(argv[1]);
  if (algorithm != 1 && algorithm != 2 && algorithm != 3) {
    std::cout << "wrong algorithm" << endl;
    exit(1);
  }
  partitionCount = atoi(argv[2]);
  if(partitionCount > MAXPARTNO)
  {
    std::cout << "Maximum allowed number of partitions is " << MAXPARTNO << ", setting partition count to " << MAXPARTNO << "." << std::endl;
    partitionCount = MAXPARTNO;
  }  
  double imbal = atof(argv[3]);
  if(imbal > MAXIMBAL)
  {
    std::cout << "Maximum allowed tolerance ratio is " << MAXIMBAL << ", setting tolerance ratio to " << MAXIMBAL << "." << std::endl;
    imbal = MAXIMBAL;
  }
    
  int slackValue = atoi(argv[4]);  
  fileName = argv[5];		
  randomizationCount = atoi(argv[6]);

  int byteSize;
  int hashCount;
  Algorithms* Partitioner;
  if(algorithm == 3)
  { 
    if(argc > 8)
    {
      byteSize = atoi(argv[7]);
      hashCount = atoi(argv[8]);
      Partitioner = new Algorithms(fileName, byteSize, hashCount);
    }
    else
    {
      std::cout << "Missing info for BF " << std::endl;
    }
  }
  else
    Partitioner = new Algorithms(fileName);
  
  for(int i = 0; i < randomizationCount; i++)
  {
    auto start = std::chrono::high_resolution_clock::now();
    Partitioner->partition(algorithm, partitionCount, slackValue, imbal);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Run: "<< i + 1 << " - Duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
  }
  
  delete Partitioner;
  return 0;	
}
