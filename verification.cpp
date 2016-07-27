#include <vector>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <fstream>
#define _USE_MATH_DEFINES



int main(){

  bool output_array[8]{0,0,0,0,0,0,0,0};

  for (auto g : output_array){std::cout << g;}
  std::cout << std::endl;

  for(int z = 0; z != 8; ++z){
  for (int n = 0; n != 8; ++n){
    output_array[n] = (z%2)*(n%2);
  }
  for (auto g : output_array){std::cout << g;}
  std::cout << std::endl;

}
  return 0;
}
