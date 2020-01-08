#include "partitioning.cpp"
#include <chrono>

int main()
{
	std::ofstream output;
	output.open("output.txt", std::ios_base::app);
	while (true)
	{
		std::string fileName;
		int partitionCount, randomizationCount;
		int state;
		std::cout << "[1]: N2P - [2]: P2N - [3]: BF - Other: Exit" << std::endl;
		std::cout << "Choice: "; std::cin >> state;
		if (state != 1 && state != 2 && state != 3)
		{
			break;
		}
		std::cout << "Enter partition count: "; std::cin >> partitionCount;
		std::cout << "Enter matrix name: "; std::cin >> fileName;
		std::cout << "Enter randomization count: "; std::cin >> randomizationCount;
		switch (state)
		{			
			case 1:
			{
				for (size_t i = 0; i < randomizationCount; i++)
				{
					Algorithms n2pPartitioner(fileName, partitionCount);
					auto start = std::chrono::high_resolution_clock::now();
					n2pPartitioner.LDGn2p();
					auto end = std::chrono::high_resolution_clock::now();
					output << "N2P with " << fileName << ", " << partitionCount << ": " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
					int cuts = n2pPartitioner.calculateCuts();
					output << "Cuts with N2P:" << cuts << std::endl;
					output << std::endl;
				}
				break;
			}				
			case 2:
			{
				for (size_t i = 0; i < randomizationCount; i++)
				{
					Algorithms p2nPartitioner(fileName, partitionCount);
					auto start = std::chrono::high_resolution_clock::now();
					p2nPartitioner.LDGp2n();
					auto end = std::chrono::high_resolution_clock::now();
					output << "P2N with " << fileName << ", " << partitionCount << ": " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;					
					output << std::endl;
				}
				break;
			}
			case 3:
			{
				int byteSize, hashCount;
				std::cout << "Enter BF size in bytes: "; std::cin >> byteSize;
				std::cout << "Enter hash function count for BF: "; std::cin >> hashCount;
				Algorithms bfPartitioner(fileName, partitionCount, byteSize, hashCount);
				auto start = std::chrono::high_resolution_clock::now();
				bfPartitioner.LDGBF();
				auto end = std::chrono::high_resolution_clock::now();
				output << "Regular BF with " << fileName << ", " << partitionCount << ": " << std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count() << "s" << std::endl;
				output << std::endl;
				break;
			}
			
			default:
				break;				
		}
	}
	output.close();

	return 0;	
}