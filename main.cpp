#include "partitioner.cpp"

int main(int argc, char** argv) {

  if (argv[1] == "h")
    {
      std::cout << "--------------------\n";
      std::cout << "arg1: Algorithm\n";
      std::cout << "0 - Random Partition\n";
      std::cout << "1 - P2N\n";
      std::cout << "2 - N2P\n";
      std::cout << "3 - N2Pi\n";
      std::cout << "4 - BF\n";
      std::cout << "5 - BF2\n";
      std::cout << "6 - BF3\n";
      std::cout << "7 - Multilayer BF\n";
      std::cout << "--------------------\n";
      std::cout << "arg2: Partition Count/n";
      std::cout << "arg3: Imbalance\n";
      std::cout << "arg4: Slack\n";
      std::cout << "arg5: Filename\n";
      std::cout << "arg6: Randomization Count\n";
      std::cout << "--------------------\n";
      std::cout << "If Alg is BF:\n";
      std::cout << "arg7: Byte size\n";
      std::cout << "arg8: Hash Count\n";
      std::cout << "--------------------\n";
      std::cout << "If Alg is N2Pi:\n";
      std::cout << "arg7: i\n";
    }
  
  else
    {
      
      int algorithm = atoi(argv[1]);
      
      if (algorithm != 0 && algorithm != 1 && algorithm != 2 && algorithm != 3 && algorithm != 4 && algorithm != 5 && algorithm != 6 && algorithm != 7) {
	std::cout << "wrong algorithm" << endl;
	exit(1);
      }
      
      int partitionCount = atoi(argv[2]);
      if (partitionCount > MAXPARTNO)
	{
	  std::cout << "Maximum allowed number of partitions is " << MAXPARTNO << ", setting partition count to " << MAXPARTNO << "." << std::endl;
	  partitionCount = MAXPARTNO;
	}
      
      
      std::string fileName;
      
      int byteSize;
      int hashCount;
      int randomizationCount;
      int num_layer;
      Partitioner* partitioner;
      
      if (algorithm == 0)
	{
	  fileName = argv[3];
	  randomizationCount = atoi(argv[4]);
	  Partitioner randomPartitioner(fileName);
	  
	  for (int i = 0; i < randomizationCount; i++)
	    {
	      randomPartitioner.RandomPartition(partitionCount, i);
	    }
	}
      else
	{
	  double imbal = atof(argv[3]);
	  if (imbal > MAXIMBAL)
	    {
	      std::cout << "Maximum allowed tolerance ratio is " << MAXIMBAL << ", setting tolerance ratio to " << MAXIMBAL << "." << std::endl;
	      imbal = MAXIMBAL;
	    }
	  
	  int slackValue = atoi(argv[4]);
	  fileName = argv[5];
	  randomizationCount = atoi(argv[6]);
	  
	  int byteSize;
	  int hashCount;
	  int num_layer;
	  Partitioner* partitioner;
	  
	  if (algorithm == 4 || algorithm == 5 || algorithm == 6 || algorithm == 7)
	    {
	      if (argc > 8)
		{
		  byteSize = atoi(argv[7]);
		  hashCount = atoi(argv[8]);
		  if (argc > 9)
		    num_layer = atoi(argv[9]);
		  else
		    num_layer = -1;
		  partitioner = new Partitioner(fileName, byteSize, hashCount, num_layer);
		}
	      
	      else
		{
		  std::cout << "Missing info for BF " << std::endl;
		  return 1;
		}
	    }
	  else
	    partitioner = new Partitioner(fileName);
	  
	  int i = 0; //doesn't do anything unless alg = 3
	  
	  if(algorithm == 3)
	    {
	      i = atoi(argv[7]);
	    }
	  
	  for (int i = 0; i < randomizationCount; i++)
	    {
	      partitioner->partition(algorithm, partitionCount, slackValue, i + 1, imbal, i);
	    }
	  
	  delete partitioner;
	}
    }
  
  return 0;
}
