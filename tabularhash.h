#ifndef TABULAR_HASH
#define TABULAR_HASH

class TabularHash {
public:
	TabularHash();
	TabularHash(unsigned int seed);
	unsigned int Hash(int key1, int key2);
private:	
	unsigned int Table[4][256];
	unsigned int bitMask(unsigned int n);
};
#endif // !TABULAR_HASH
