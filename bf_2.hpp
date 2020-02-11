#include<bitset>
#include<cstring>

//int primes[4] = {1129, 2063, 3217, 4093};
//int primes[4] = {1129247, 2063522, 67959190, 80932000};
int primes[4] = {1129247, 3063522, 5752225, 8191999};
//int primes[4] = {11290, 20630, 320017, 390093};
//int primes[16] = {643873, 98507, 181277, 2536, 325231, 340933, 416401, 519371, 647033, 735107, 837461, 917239, 990469, 1060177, 1136411, 1269173};

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
    int hash_val = ((a*val + b)%p);//%4096;
    
    return hash_val;
  }
  
  
};

void set_bit(uint32_t* &ints, int big_loc, int small_loc){
  //std::cout << "I'm setting loc: " << loc << std::endl;
  //uint32_t m = 1;
  ints[big_loc] |= (1U << (small_loc-1));
}

void clear_bit(uint32_t* &ints, int big_loc, int small_loc){
  ints[big_loc] &= ~(1U << (small_loc-1));
}

bool check_bit(uint32_t* &ints, int big_loc, int small_loc){
  //return !((ints[big_loc] >> small_loc) & 1);
  return (ints[big_loc] & (1 << (small_loc - 1)));
}

void toggle_bit(uint32_t* &ints, int big_loc, int small_loc){
  if(!check_bit(ints, big_loc, small_loc))
    ints[big_loc] ^= (1U << (small_loc-1));
}

#define NOBITS 8192000
#define NOHASHES 3
struct BloomFilter {
private:
  uint32_t* ints;
  uint32_t no_bits;
  int bf_id;
  uint32_t no_ints;
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
      //bits[index] = 1;
      int big_loc = index/32;
      int small_loc = index - (big_loc*32);
      //std::cout << "Before setting bit, big_loc: " << big_loc << " small_loc: " << small_loc << " index: " << index << std::endl;
      toggle_bit(ints, big_loc, small_loc);
    }
  }
  
  bool query(int val){
    for(int k = 0; k < NOHASHES; k++){
      int index = hash(val,k);
      int big_loc = index/32;
      int small_loc = index - (big_loc*32);
      //std::cout << "Querying, big_loc: " << big_loc << " small_loc: " << small_loc << " index: " << index << std::endl;
      bool ret_val = check_bit(ints, big_loc, small_loc);
      //std::cout << "check_bit returned: " << ret_val << " index: " << index <<std::endl;
      //if(bits[index] == 0)
      if(!check_bit(ints, big_loc, small_loc))
	return false;
    }
    //std::cout << "val: " << val << " returning true " << std::endl;
    return true;
  }
  
  BloomFilter(int no_bits, int bf_id) : no_bits(no_bits), bf_id(bf_id){
    //no_ints = ceil(no_bits / (sizeof(int) * 8.0)) + 1;
    no_ints = NOBITS/32;
    //std::cout << "No ints: " << no_ints << std::endl;
    //exit(1);
    ints = new uint32_t[no_ints];
    memset(ints, 0, sizeof(uint32_t) * no_ints);
    
    for(int i = 0; i < NOHASHES; i++){
      hashes[i] = new ax_b_hash(i, bf_id);
    }
  }
  
  BloomFilter()
  {}
  
  ~BloomFilter() {
    //delete [] ints;
    //delete [] hashes;
  }
};
