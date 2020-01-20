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

  double imbal = atof(argv[3]);  
  fileName = argv[4];		
  randomizationCount = atoi(argv[5]);

  int byteSize;
  int hashCount;
  Algorithms* Partitioner;
  if(algorithm == 3)
  { 
    if(argc > 7)
    {
      byteSize = atoi(argv[6]);
      hashCount = atoi(argv[7]);
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
    Partitioner->partition(algorithm, partitionCount, imbal);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Run: "<< i + 1 << " - Duration: " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
  }
  
  return 0;	
}
