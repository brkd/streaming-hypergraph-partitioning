#include <tgmath.h>

class mlbf{
  
private:
  int num_layers;
  int num_partitions;
  int byte_size;
  int layer_size;
  int num_filters;
  BloomFilter_OT** bf_heap;
  
public:
  mlbf(int NUM_LAYERS, int NUM_PARTITIONS, int HASH_COUNT){
    byte_size = 8000000;
    num_layers = NUM_LAYERS;
    layer_size = byte_size / num_layers;
    
    num_filters = 0;
    
    for(int i = 1; i <= num_layers; i++){
      
      num_filters += pow(2, i);
      
    }
    
    bf_heap = new BloomFilter_OT*[num_filters+1];
    
    int filter_ctr = 1;
    int filter_in_layer = 1;

    for(int i = 1; i <= num_layers; i++){
      std::cout << "Layer size: " << layer_size << std::endl;
      filter_in_layer *= 2;
      
      for(int f = 0; f < filter_in_layer; f++){
      bf_heap[filter_ctr] = new BloomFilter_OT(layer_size/pow(2,i), HASH_COUNT, i);
      std::cout << "Filter no: " << filter_ctr <<" with size: " << layer_size/pow(2,i) << std::endl;
      filter_ctr++;
      }
    }
    
    exit(1);

  }
};
