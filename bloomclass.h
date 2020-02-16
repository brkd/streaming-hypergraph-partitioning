#ifndef _BLOOM_H
#define _BLOOM_H
#include <iostream>
#include <math.h>
#include <cmath>
#include <random>
#include <vector>
#include <functional>

using namespace std;

template <class T1, class T2>
  class Bloom {
 private:
  int modPrime;
  uint64_t size;
  uint64_t hash_count;
  vector<unsigned char> bloomFilter;
  unsigned int basicHash(T1 item1, T2 item2, unsigned int seed);
  //unsigned int tabularHash(T1 item1, T2 item2);

 public:
  Bloom(uint64_t byteSize, uint64_t hashCount);
  void insert(T1 item1, T2 item2);
  bool contains(T1 item1, T2 item2);	
};


template<class T1, class T2>
  unsigned int Bloom<T1, T2>::basicHash(T1 item1, T2 item2, unsigned int seed) {
  return ((23*item1) + (47*item2) + (101*item1)*item2 + seed) % modPrime;
  //return ((431*item1) + (439*item2) + (521*item1)*item2 + seed) % modPrime;
}

template <class T1, class T2>
  Bloom<T1, T2>::Bloom(uint64_t byteSize, uint64_t hashCount) {
  hash_count = hashCount;
  size = byteSize;
  modPrime = 2954646407;
  
  for (int i = 0; i < size; i++) {
    bloomFilter.push_back(0);
  }
}

template <class T1, class T2>
  void Bloom<T1, T2>::insert(T1 item1, T2 item2) {
  uint32_t cur_bit;
  uint32_t cur_char;
  uint8_t one = 1;
  
  for (unsigned int i = 0; i < hash_count; i++)
    {
      cur_bit = basicHash(item1, item2, i) % size;
      cur_char = cur_bit / 8;
      cur_bit = cur_bit % 8;
      
      bloomFilter[cur_char] = bloomFilter[cur_char] | (one << cur_bit);
    }
}

template <class T1, class T2>
  bool Bloom<T1, T2>::contains(T1 item1, T2 item2) {
  uint32_t cur_bit;
  uint32_t cur_char;
  uint8_t one = 1;
  
  for (unsigned int i = 0; i < hash_count; i++)
    {
      cur_bit = basicHash(item1, item2, i) % size;
      cur_char = cur_bit / 8;
      cur_bit = cur_bit % 8;
      
      if ((bloomFilter[cur_char] & (one << cur_bit)) == 0)
	return false;
    }
  return true;
}
#endif
