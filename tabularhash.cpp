#include "tabularhash.h"
#include <chrono>
#include <functional>
#include <random>
#include <iostream>
#include <stdio.h>
#include <tchar.h>
#include <limits.h>

TabularHash::TabularHash()
{
	std::mt19937::result_type seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	auto rand_32bit = std::bind(std::uniform_int_distribution<unsigned int>(0, UINT_MAX), std::mt19937(seed));

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			Table[i][j] = rand_32bit();
		}
	}
}

TabularHash::TabularHash(unsigned int seed)
{
	auto rand_32bit = std::bind(std::uniform_int_distribution<unsigned int>(0, UINT_MAX), std::mt19937(seed));

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			Table[i][j] = rand_32bit();
		}
	}
}

unsigned int TabularHash::Hash(int key)
{
	unsigned int result = 0;
	unsigned int hashArray[4];

	unsigned int mask = bitMask(8);

	unsigned int first = key & bitMask(8);
	unsigned int second = (key >> 8)  & mask;
	unsigned int third = (key >> 16) & mask;
	unsigned int fourth = (key >> 24) & mask;

	hashArray[0] = Table[0][first];
	hashArray[1] = Table[1][second];
	hashArray[2] = Table[2][third];
	hashArray[3] = Table[3][fourth];

	for (int i = 0; i < 4; i++)
	{
		result ^= hashArray[i];
	}

	return result;
}

unsigned int TabularHash::bitMask(unsigned int n)
{
	unsigned int mask = 0;
	mask = ((1 << n) - 1);

	return mask;
}
