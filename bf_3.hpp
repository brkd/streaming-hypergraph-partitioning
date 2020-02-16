#include<bitset>
#include<cstring>


////////int primes[4] = {1129247, 3063522, 5752225, 8191999};
//int OTprimes[29] = {11293, 301363, 54752, 8191, 23063, 531752, 845191,14555129, 30454563, 575211, 81191, 13063, 75211, 81191,132129, 3354063, 545752, 448191, 434063, 55752, 8542191,112900, 3063, 22522, 4319121, 13063, 57521, 8191};

//PRIME//

bool isPrime(int n)  
{  
  // Corner cases  
  if (n <= 1)  return false;  
  if (n <= 3)  return true;  
    
  // This is checked so that we can skip   
  // middle five numbers in below loop  
  if (n%2 == 0 || n%3 == 0) return false;  
    
  for (int i=5; i*i<=n; i=i+6)  
    if (n%i == 0 || n%(i+2) == 0)  
      return false;  
    
  return true;  
}  
  
// Function to return the smallest 
// prime number greater than N 
int nextPrime(int N) 
{ 
  
  // Base case 
  if (N <= 1) 
    return 2; 
  
  int prime = N; 
  bool found = false; 
  
  // Loop continuously until isPrime returns 
  // true for a number greater than n 
  while (!found) { 
    prime++; 
  
    if (isPrime(prime)) 
      found = true; 
  } 
  
  return prime; 
} 

//PRIME//



struct OTax_b_hash{
private:
  int a;
  int b;
  int p;
  int size;

public:
  
  OTax_b_hash()
  {}
  
  OTax_b_hash(int num, int bf_id, int csize, int hash_count){
    a = 712 + bf_id + num*322;
    b = 1909 + bf_id + num*4043;

    //std::cout << "num: " << num << " size: " << csize <<std::endl;
    
    size = csize;
    
    for(int i = 0; i < hash_count; i++){
      if(num == i){
	p = nextPrime(ceil(2*size/(i+3)));
	//p = nextPrime(size/(i+1));
	//p+=23;
      }
    }
    
    if(num == hash_count-1){
      this->p = nextPrime(csize*0.85);
      while(this->p < csize-(csize*0.05))
	{
	  this->p = nextPrime(this->p);
	}
    }
    //this->p = csize-1;
    
    //std::cout << "hash_no: " << num << " p is: " << p << std::endl;
    //std::cout << "a: " << a << " b: " << b << std::endl;
    //delete[] OTPrimes;
  }
  
  int OThash(int val1, int val2){
   

    int pre_sum_val = (val1 * val2);
    //xint sum_val = std::hash<int>{}(pre_sum_val);
    int hash_val = (a*pre_sum_val + b)%p;//;%(size/4));//%size;//%4096;
  
    //std::cout << "p: " << p << std::endl;

    if(hash_val < 0){
      hash_val *= -1;
      hash_val = hash_val%p;//%(size/4);
    }
    
    //std::cout << "val1: " << val1 << " val2: " << val2 <<" Returning: " << hash_val << std::endl;
    return hash_val;
  }
  
  
};

/*
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
*/

struct BloomFilter_OT {
private:
  uint32_t* ints;
  uint32_t no_bits;
  int bf_id;
  int hash_count;
  uint32_t no_ints;
  std::bitset<NOBITS> bits;
  OTax_b_hash** hashes;//[NOHASHES];
  
public:
  int hash(int val1, int val2, int k) {
    int index = hashes[k]->OThash(val1, val2);
    return index;
  }
  
  void insert(int val1, int val2){
    for(int k = 0; k < hash_count; k++){
      int index = hash(val1, val2, k);
      //bits[index] = 1;
      int big_loc = index/32;
      int small_loc = index - (big_loc*32);
      set_bit(ints, big_loc, small_loc);
    }
  }
  
  bool query(int val1, int val2){
    for(int k = 0; k < hash_count; k++){
      int index = hash(val1, val2, k);
      int big_loc = index/32;
      int small_loc = index - (big_loc*32);
      //if(bits[index] == 0)
      check_bit(ints, big_loc, small_loc);
      if(!check_bit(ints, big_loc, small_loc))
	return false;
    }
    return true;
  }
  
  BloomFilter_OT(int no_bits, int hashCount, int bf_id) : no_bits(no_bits), bf_id(bf_id){
    no_ints = no_bits/32;
    hash_count = hashCount;
    ints = new uint32_t[no_ints];
    memset(ints, 0, sizeof(uint32_t) * no_ints);
    
    hashes = new OTax_b_hash*[hashCount];
    for(int i = 0; i < hashCount; i++){
      hashes[i] = new OTax_b_hash(i, bf_id, no_bits, hashCount);
    }
  }
  
  BloomFilter_OT()
  {}
  
  ~BloomFilter_OT() {
    delete [] ints;
    //delete [] hashes;
  }
};
