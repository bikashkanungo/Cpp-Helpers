#include <cstdlib>

#include <iostream>

#include <string>

#include <vector>


int main()
{

  std::vector<std::vector<double> > c(10);
  for(unsigned int i = 0; i < c.size(); ++i)
    c[i] = std::vector<double>(2*i+1,i);

  for(unsigned int i = 0; i < c.size(); ++i)
    c[i].push_back(i+1);

  for(unsigned int i = 0; i < c.size(); ++i)
  {
    for(unsigned int j = 0; j < c[i].size(); ++j)
      std::cout << c[i][j] << " ";
    std::cout << std::endl;
  }
}

