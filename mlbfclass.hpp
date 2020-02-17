#include <tgmath.h>

class mlbf{
  
private:
  int num_layers;
  int num_partitions;
  int byte_size;
  int layer_size;
  int num_filters;
  int hash_count;
  BloomFilter_OT** bf_heap;
  
public:
  mlbf(int NUM_LAYERS, int NUM_PARTITIONS, int HASH_COUNT, int NO_BYTES){
    hash_count = HASH_COUNT;
    byte_size = NO_BYTES;
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
      //std::cout << "Layer size: " << layer_size << std::endl;
      filter_in_layer *= 2;
      
      for(int f = 0; f < filter_in_layer; f++){
      bf_heap[filter_ctr] = new BloomFilter_OT(layer_size/pow(2,i), HASH_COUNT, i);
      //std::cout << "Filter no: " << filter_ctr <<" with size: " << layer_size/pow(2,i) << std::endl;
      filter_ctr++;
      }
    }
  }
  
    ~mlbf(){
      delete[] bf_heap;
    }
  
  void layer_series(int layer, int cords[]){
    
    int first = 1;
    int last = 0;

    for(int i = 0; i < layer-1; i++){
      first += pow(2, i);
    }

    last = first + pow(2,layer);
    
    //int cords[2] = {first, last};
    cords[0] = first;
    cords[1] = last;
    std::cout << "In series: " << first << " " << last << std::endl;

    return cords;
  }
  
  bool recursive_insert(int layer, int child, int* range, int edge, int part){
    
    if(layer == num_layers){
      bf_heap[child]->insert(edge, part);
      return true;
    }
    
    if(child == 0){
      int filter_num = 0;
      int child = -1;
      
      int start = 0;
      int end = num_partitions-1;
      
      if(part < (start+end)/2){
	filter_num = 1;
	start = 0;
	end = floor((num_partitions-1)/2);
	if(part < (start+end)/2){
	  child = 2*filter_num + 1;
	  end = (start+end)/2;
	}
	else{
	  child = 2*filter_num + 2;
	  start = (start+end)/2;
	}
	int range[2] = {start, end};
      }
      else{ 
	filter_num = 2;
	start = floor((num_partitions-1)/2);
	end = num_partitions-1;
	if(part < (start+end)/2){
	  child = 2*filter_num + 1;
	  end = (start+end)/2;
	}
	else{
	  child = 2*filter_num + 2;
	  start = (start+end)/2;
	}
	int range[2] = {start, end};
      }
      
      if(recursive_insert(layer+1, child, range, edge, part)){
	bf_heap[filter_num]->insert(edge, part);
	return true;
      }
    }

      int start = range[0];
      int end = range[1];
      int r_child = -1;
            
      
      if(part < (start+end)/2){
	r_child = child * 2 + 1;
	end = (start+end)/2;
	int range[2] = {start, end};
      }
      else{
	r_child = child * 2 + 2;
	start = (start+end)/2;
	int range[2] = {start, end};
      }

      if(recursive_insert(layer+1, r_child, range, edge, part)){
	bf_heap[child]->insert(edge, part);
	return true;
      }   
    
  }

  void insert(int edge, int part){
    //std::cout << "Inserting" << std::endl;
    int range[2] = {-1,-1};
    recursive_insert(1, 0, range, edge, part); //Start from layer 1
  }

  bool recursive_query(int edge, int part, int child, int layer, int* range){
    
    if(layer == num_layers)
      return bf_heap[child]->query(edge, part);

    if(bf_heap[child]->query(edge,part)){
      int r_child = -1;
      int r_range[2] = {-1,-1};
      int start = range[0];
      int end = range[1];
      
      if(part < (start+end)/2){
	r_child = (child*2)+1;
	r_range[0] = start;
	r_range[1] = (start+end)/2;
	return recursive_query(edge, part, r_child, layer+1, r_range);
      }
      else{
	r_child = (child*2)+2;
	r_range[0] = (start+end)/2;
	r_range[1] = end;
	return recursive_query(edge, part, r_child, layer+1, r_range);
      }
      
      
    }
    else{ 
      return false;
    }

  }

  bool query(int edge, int part){
    
    int location_in_heap = -1;
    
    if(part < num_filters/2){
      location_in_heap = 1;
      int range[2] = {0, num_filters/2};
      return recursive_query(edge, part, location_in_heap, 1, range);
    }
    else{
      location_in_heap = 2;
      int range[2] = {num_filters/2, num_filters};
      return recursive_query(edge, part, location_in_heap, 1, range);
    }
    
    
  }

   
};
