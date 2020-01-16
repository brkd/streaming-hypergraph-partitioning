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
  if(state == 3) { 
    if(argc > 7) {
      byteSize = atoi(argv[6]);
      hashCount = atoi(argv[7]);
    } else {
      std::cout << "Missing info for BF " << std::endl;
    }
  }

  if(algorithm == 1) {
      Algorithms n2pPartitioner(fileName);
      for (size_t i = 0; i < randomizationCount; i++) {
	auto start = std::chrono::high_resolution_clock::now();
	n2pPartitioner.LDGn2p(partitionCount, imbal);
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "N2P with " << fileName << ", " << partitionCount << ": " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
	int cuts = n2pPartitioner.calculateCuts();
	std::cout << "Cuts with N2P:" << cuts << std::endl;
	std::cout << std::endl;
      }
  } else if(algorithm == 2) {
    Algorithms p2nPartitioner(fileName);
    for (size_t i = 0; i < randomizationCount; i++) {
      auto start = std::chrono::high_resolution_clock::now();
      p2nPartitioner.LDGp2n(partitionCount, imbal);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "P2N with " << fileName << ", " << partitionCount << ": " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;					
      std::cout << std::endl;
    }
  } else if(algorithm == 3) {
    Algorithms bfPartitioner(fileName, byteSize,hashCount);
    for (size_t i = 0; i < randomizationCount; i++) {
      auto start = std::chrono::high_resolution_clock::now();
      bfPartitioner.LDGBF(partitionCount, imbal);
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Regular BF with " << fileName << ", " << partitionCount << ": " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
      std::cout << std::endl;
    }
  }
  return 0;	
}
