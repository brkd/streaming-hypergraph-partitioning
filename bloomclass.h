#pragma once
#include <iostream>
#include <math.h>
#include <cmath>
#include <random>
#include <vector>
#include <functional>

//#include "tabularhash.h"
using namespace std;

template <class T1, class T2>
class Bloom {
private:
	int modPrime;
	uint64_t size;
	uint64_t hash_count;
	vector<unsigned char> bit_array;
	unsigned int basicHash(T1 item1, T2 item2, unsigned int seed);
	//TabularHash* tabular_hash;

public:
	Bloom(uint64_t byteSize, uint64_t hashCount);
	void insert(T1 item1, T2 item2);
	bool contains(T1 item1, T2 item2);
	/*void calculate_size(uint64_t n, float fp_prob);
	void calculate_hash_count(uint64_t expected_amount);*/
};


template<class T1, class T2>
inline unsigned int Bloom<T1, T2>::basicHash(T1 item1, T2 item2, unsigned int seed)
{
  //auto rand_32bit = bind(std::uniform_int_distribution<unsigned int>(0, modPrime), std::mt19937(seed));

  //return (rand_32bit()*item1 + rand_32bit()*item2 + rand_32bit()*item1*item2) % modPrime;

  return item1*item2+seed%6;
}

template <class T1, class T2>
Bloom<T1, T2>::Bloom(uint64_t byteSize, uint64_t hashCount) {
	hash_count = hashCount;
	size = byteSize;
	modPrime = 2147483647;

	/*calculate_size(n, fp_prob);
	calculate_hash_count(n);*/

	for (int i = 0; i < size; i++) {
		bit_array.push_back(0);
	}
}

template <class T1, class T2>
void Bloom<T1, T2>::insert(T1 item1, T2 item2) {
	uint32_t cur_bit;
	uint32_t cur_char;
	uint8_t one = 1;

	for (unsigned int i = 0; i < hash_count; i++)
	{
		//tabular_hash = new TabularHash(i);

		cur_bit = basicHash(item1, item2, i) % size;
		cur_char = cur_bit / 8;
		cur_bit = cur_bit % 8;

		bit_array[cur_char] = bit_array[cur_char] | (one << cur_bit);

		//delete(tabular_hash);
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

		if ((bit_array[cur_char] & (one << cur_bit)) == 0)
			return false;
	}
	return true;
}
//
//template <class T1, class T2>
//void Bloom<T1, T2>::calculate_size(uint64_t n, float pb) {
//	size = -1 * (n * log(pb)) / pow(log(2), 2);
//}

//template <class T1, class T2>
//void Bloom<T1, T2>::calculate_hash_count(uint64_t n) {
//	hash_count = (size / n) * log(2);
//}
