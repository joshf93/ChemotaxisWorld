#include <vector>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <fstream>
#define _USE_MATH_DEFINES



inline int is_set(const uint32_t &num, const int &idx){
  return ((num & (1 << idx)) >> idx);
}

int main(){

  for (uint32_t j = 0; j != 1024; ++j){
  for(int i = 0; i != 32; ++i){
    std::cout << is_set(j, i);
  }
  std::cout << std::endl;
}
  std::cout << std::endl;

  return 0;
}
