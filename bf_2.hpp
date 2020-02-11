#include<bitset>
#include<cstring>


////////int primes[4] = {1129247, 3063522, 5752225, 8191999};
int primes[29] = {11293, 301363, 54752, 8191, 23063, 531752, 845191,14555129, 30454563, 575211, 81191, 13063, 75211, 81191,132129, 3354063, 545752, 448191, 434063, 55752, 8542191,112900, 3063, 22522, 4319121, 13063, 57521, 8191};

struct ax_b_hash{
private:
  int a;
  int b;
  int p;
  int size;

public:
  
  ax_b_hash()
  {}
  
  ax_b_hash(int num, int bf_id, int csize){
    a = 712 + bf_id + num*322;
    b = 1909 + bf_id + num*4043;
    p = primes[num];
    size = csize;
  }
  
  int hash(int val){

    int hash_val = (a*val + b)%p%(size/4);//;%(size/4));//%size;//%4096;
  
    if(hash_val < 0){
      hash_val *= -1;
      hash_val = hash_val%p%(size/4);//%(size/4);
    }
    
    return hash_val;
  }
  
  
};

void set_bit(uint32_t* &ints, int big_loc, int small_loc){
  ints[big_loc] |= (1U << (small_loc-1));
}

void clear_bit(uint32_t* &ints, int big_loc, int small_loc){
  ints[big_loc] &= ~(1U << (small_loc-1));
}

bool check_bit(uint32_t* &ints, int big_loc, int small_loc){
  return (ints[big_loc] & (1 << (small_loc - 1)));
}

void toggle_bit(uint32_t* &ints, int big_loc, int small_loc){
  if(!check_bit(ints, big_loc, small_loc))
    ints[big_loc] ^= (1U << (small_loc-1));
}

#define NOBITS 8192
#define NOHASHES 3
struct BloomFilter {
private:
  uint32_t* ints;
  uint32_t no_bits;
  int bf_id;
  int hash_count;
  uint32_t no_ints;
  std::bitset<NOBITS> bits;
  ax_b_hash** hashes;//[NOHASHES];
  
public:
  int hash(int val, int k) {
    int index = hashes[k]->hash(val);
    return index;
  }
  
  void insert(int val){
    for(int k = 0; k < hash_count; k++){
      int index = hash(val, k);
      //bits[index] = 1;
      int big_loc = index/32;
      int small_loc = index - (big_loc*32);
      set_bit(ints, big_loc, small_loc);
    }
  }
  
  bool query(int val){
    for(int k = 0; k < hash_count; k++){
      int index = hash(val,k);
      int big_loc = index/32;
      int small_loc = index - (big_loc*32);
      //if(bits[index] == 0)
      check_bit(ints, big_loc, small_loc);
      if(!check_bit(ints, big_loc, small_loc))
	return false;
    }
    return true;
  }
  
  BloomFilter(int no_bits, int hashCount, int bf_id) : no_bits(no_bits), bf_id(bf_id){
    no_ints = no_bits/32;
    hash_count = hashCount;
    ints = new uint32_t[no_ints];
    memset(ints, 0, sizeof(uint32_t) * no_ints);
    
    hashes = new ax_b_hash*[hashCount];
    for(int i = 0; i < hashCount; i++){
      hashes[i] = new ax_b_hash(i, bf_id, no_bits);
    }
  }
  
  BloomFilter()
  {}
  
  ~BloomFilter() {
    //delete [] ints;
    //delete [] hashes;
  }
};
