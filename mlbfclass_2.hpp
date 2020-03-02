#include <tgmath.h>

class mlbf_2{
  
private:
  int num_layers;
  int num_partitions;
  int byte_size;
  int layer_size;
  int num_filters;
  int hash_count;
  BloomFilter_OT_OV** bf_heap;
  
public:
  mlbf_2(int NUM_LAYERS, int NUM_PARTITIONS, int HASH_COUNT, int NO_BYTES){
    hash_count = HASH_COUNT;
    byte_size = NO_BYTES;
    num_layers = NUM_LAYERS;
    num_partitions = NUM_PARTITIONS;
    layer_size = byte_size / num_layers;
    
    num_filters = 0;
    
    for(int i = 1; i <= num_layers; i++){
      
      num_filters += pow(2, i);
      
    }
    
    bf_heap = new BloomFilter_OT_OV*[num_filters+1];
    
    int filter_ctr = 1;
    int filter_in_layer = 1;

    for(int i = 1; i <= num_layers; i++){
      //std::cout << "Layer size: " << layer_size << std::endl;
      filter_in_layer *= 2;
      
      for(int f = 0; f < filter_in_layer; f++){
	bf_heap[filter_ctr] = new BloomFilter_OT_OV(layer_size/pow(2,i)*8, HASH_COUNT, i);
	//std::cout << "Filter no: " << filter_ctr <<" with size: " << layer_size/pow(2,i) << std::endl;
	filter_ctr++;
      }
    }
  }
  
    ~mlbf_2(){
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
    //std::cout << "In series: " << first << " " << last << std::endl;

    return cords;
  }
  
  bool recursive_insert(int layer, int child, int* range, int edge, int part){

    //std::cout << "Part: " << part << " edge: " << edge << " child: " << child << std::endl;
    
    if(child == 0){
      int filter_num = -1;
      int r_child = -1;
      
      int start = 0;
      int end = num_partitions-1;
      
      int range[2] = {-1, -1};
      
      if(part <= (start+end)/2){
	filter_num = 1;
	start = 0;
	end = (num_partitions-1)/2;
	//r_child = 2*filter_num + 1;  
      	range[0] = start;
	range[1] = end;
	if(part <= (start+end)/2){
	  r_child = 2*filter_num + 1;  
	  range[1] = (start+end)/2;
	}
	else{
	  r_child = 2*filter_num + 2;
	  range[0] = (start+end)/2 + 1;
	}
      }
      else{ 
	filter_num = 2;
	start = ((num_partitions-1)/2)+1;
	end = num_partitions-1;
	range[0] = start;
	range[1] = end;
	//r_child = 2*filter_num + 2;
	if(part <= (start+end)/2){
	  r_child = 2*filter_num + 1;  
	  range[1] = (start+end)/2;
	}
	else{
	  r_child = 2*filter_num + 2;
	  range[0] = (start+end)/2 + 1;
	}
      }
      
      
      if(num_layers == 1){
	//std::cout << "Part: " << part << "-->Inserting filter num: " << filter_num << std::endl;
	bf_heap[filter_num]->insert(edge);
	return true;
      }
      
      else if(recursive_insert(layer+1, r_child, range, edge, part)){
	//std::cout << "Inserting filter num: " << filter_num << std::endl;
	bf_heap[filter_num]->insert(edge);
	return true;
      }
    }
    
    if(layer == num_layers){
      bf_heap[child]->insert(edge);
      //std::cout << "Part: " << part << " start: " << range[0] << " end: " << range[1] << std::endl;
      //std::cout << "Insert: Part " << part << " chose filter no: " << child <<std::endl;
      return true;
    }
    
    int start = range[0];
    int end = range[1];
    int r_child = -1;
    int r_range[2] = {-1, -1};
    
    //std::cout << "Part: " << part << " start: " << start << " end: " << end << std::endl;
    
    
    if(part <= (start+end)/2){
      r_child = child * 2 + 1;
      end = (start+end)/2;
      r_range[0] = start;
      r_range[1] = end;
  }
    else{
      r_child = child * 2 + 2;
      start = (start+end)/2 + 1;
      r_range[0] = start;
      r_range[1] = end;
    }
    
    if(recursive_insert(layer+1, r_child, r_range, edge, part)){
      bf_heap[child]->insert(edge);
      return true;
    }   
    
  }
  
  void insert(int edge, int part){
    //std::cout << "Inserting" << std::endl;
    int range[2] = {-1,-1};
    bool status = recursive_insert(1, 0, range, edge, part); //Start from layer 1
    //std::cout << "Status: " << status << std::endl;
  }
  /*
  bool recursive_query(int edge, int part, int child, int layer, int* range){
    
    if(layer == num_layers){
      //std::cout << "Query: Part " << part << " chose filter no: " << child <<std::endl;
      return bf_heap[child]->query(edge, part);
    }
    
    if(bf_heap[child]->query(edge,part)){
      int r_child = -1;
      int r_range[2] = {-1,-1};
      int start = range[0];
      int end = range[1];
      
      //std::cout << "Part: " << part << " start: " << start << " end: " << end << std::endl;

      if(part <= (start+end)/2){
	r_child = (child*2)+1;
	r_range[0] = start;
	r_range[1] = (start+end)/2;
	return recursive_query(edge, part, r_child, layer+1, r_range);
      }
      else{
	r_child = (child*2)+2;
	r_range[0] = (start+end)/2 + 1;
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
    
    if(part <= (num_partitions-1)/2){
      location_in_heap = 1;
      int range[2] = {0, ((num_partitions-1)/2)};
      return recursive_query(edge, part, location_in_heap, 1, range);
    }
    else{
      location_in_heap = 2;
      int range[2] = {(num_partitions-1)/2+1, num_partitions - 1};
      return recursive_query(edge, part, location_in_heap, 1, range);
    }
    
    
  }
  */


  bool recursive_insert(int child, int edge, int layer){
    


  }
  
  void insert(int edge){
    recursive_insert(1, edge, 1);
    recursive_insert(2, edge, 1);
  }

  void recursive_query(int child, int edge, bool* existences, int  layer){
    //Hardcoded for 2 children per node
    
    if(bf_heap[child]->query(edge)){
      
      
      if(layer == num_layers){
	int leaf_offset = 0;
	for(int i = 0; i < layer; i++){
	  leaf_offset += pow(2,i);
	}
	
	int leaf_size = num_partitions/pow(2,layer);
	int order = child-leaf_offset;
	existences[order] = true;
	
	/*
	  int start = order*leaf_size;
	  int end = (order+1)*leaf_size;
	  
	  for(int r = start; r < end; r++){
	  existences[r] = true;
	  } 
	*/
	return;
      }
      
      for(int i = 1; i <=2; i++){
	recursive_query((child*2)+i, edge, existences, layer+1);
      }
    }
    return;
  }

  void query(int edge, bool* existences){
    recursive_query(1, edge, existences, 1);
    recursive_query(2, edge, existences, 1);
  }
   
};
