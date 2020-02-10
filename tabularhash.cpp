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

unsigned int TabularHash::Hash(int key1, int key2)
{
  unsigned int result1 = 0;
  unsigned int result2 = 0;	
  
  unsigned int hashArrayFirst[4];
  unsigned int hashArraySecond[4];
  
  unsigned int mask = bitMask(8);
  
  unsigned int first1 = key1 & mask;
  unsigned int second1 = (key1 >> 8)  & mask;
  unsigned int third1 = (key1 >> 16) & mask;
  unsigned int fourth1 = (key1 >> 24) & mask;
  
  hashArrayFirst[0] = Table[0][first1];
  hashArrayFirst[1] = Table[1][second1];
  hashArrayFirst[2] = Table[2][third1];
  hashArrayFirst[3] = Table[3][fourth1];
  
  unsigned int first2 = key2 & mask;
  unsigned int second2 = (key2 >> 8)  & mask;
  unsigned int third2 = (key2 >> 16) & mask;
  unsigned int fourth2 = (key2 >> 24) & mask;
  
  hashArraySecond[0] = Table[0][first2];
  hashArraySecond[1] = Table[1][second2];
  hashArraySecond[2] = Table[2][third2];
  hashArraySecond[3] = Table[3][fourth2];
  
  for (int i = 0; i < 4; i++)
    {
      result1 ^= hashArrayFirst[i];
      result2 ^= hashArraySecond[i];	  
    }
    
  
  return result1 ^ result2;
}

unsigned int TabularHash::bitMask(unsigned int n)
{
  unsigned int mask = 0;
  mask = ((1 << n) - 1);
  
  return mask;
}
