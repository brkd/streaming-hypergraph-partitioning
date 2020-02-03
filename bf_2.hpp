
#include<bitset>

//int primes[4] = {1129, 2063, 3217, 4093};
//int primes[4] = {11290, 20630, 320017, 390093};
int primes[16] = {643873, 98507, 181277, 25367, 325231, 340933, 416401, 519371, 647033, 735107, 837461, 917239, 990469, 1060177, 1136411, 1269173};

struct ax_b_hash{
private:
  int a;
  int b;
  int p;
  
public:
  
  ax_b_hash()
  {}
  
  ax_b_hash(int num, int bf_id){
    a = 712 + bf_id + num*322;
    b = 1909 + bf_id + num*4043;
    p = primes[num];
  }
  
  int hash(int val){
    
    int sum_val = std::hash<int>{}(val);
    //int sum_val = val1 ^ val2;
    
    int hash_val = ((a*sum_val + b)%p);//%4096;
    
    /*if(hash_val < 0){
      hash_val -= hash_val;
      hash_val /= 2;
      }*/
    
    //std::cout << "Hash val: " << hash_val << std::endl;
    return hash_val;
  }
  
  
};

#define NOBITS 8192
#define NOHASHES 3
struct BloomFilter {
private:
  //int* bits;
  uint32_t no_bits;
  int bf_id;
  //uint32_t no_ints;
  std::bitset<NOBITS> bits;
  ax_b_hash* hashes[NOHASHES];
  
public:
  int hash(int val, int k) {
    int index = hashes[k]->hash(val);
    return index;
  }
  
  void insert(int val){
    for(int k = 0; k < NOHASHES; k++){
      int index = hash(val, k);
      bits[index] = 1;
    }
  }

  bool query(int val){
    for(int k = 0; k < NOHASHES; k++){
      int index = hash(val,k);
      if(bits[index] == 0)
	return false;
    }
    return true;
  }
  
  BloomFilter(int no_bits, int bf_id) : no_bits(no_bits), bf_id(bf_id){
    //no_ints = ceil(no_bits / (sizeof(int) * 8.0)); 
    //bits = new int[no_ints];
    //memset(bits, 0, sizeof(int) * no_ints);
    
    for(int i = 0; i < NOHASHES; i++){
      hashes[i] = new ax_b_hash(i, bf_id);
    }
  }

  BloomFilter()
  {}
  
  ~BloomFilter() {
    //delete [] bits;
    //delete [] hashes;
  }
};
