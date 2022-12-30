// Example program
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

void f(double ** x,
       std::vector<double> & vec)
{

  x = (double*)malloc(vec.size()*sizeof(double));
  for(int i = 0; i < vec.size(); ++i)
    x[i] = vec[i];
  
}
int main()
{

  std::vector<double> vec(10);
  for(int i = 0; i < vec.size(); ++i)
    vec[i] = i;
  double * x;
  f(&x, vec);
  for(int i = 0; i < vec.size(); ++i)
    std::cout << x[i] << std::endl;
    
}
