#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main()
{

  int list[5];
  for(unsigned int i = 0; i < 5; ++i)
    list[i] = (i+1)*(i+1);
  int * a = &list[3] - 2;
  std::cout << *a << std::endl;
}
