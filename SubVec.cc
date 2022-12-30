#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>

void print(int a)
{

  std::cout << a << std::endl;

}

int main()
{

  std::vector<int> a(10);

  for(unsigned int i = 0; i < 10; ++i)
    a[i] = 10 - i;

  std::vector<int>::const_iterator iter1 = a.begin() + 3;
  std::vector<int>::const_iterator iter2 = a.begin() + 8;

  int length = std::distance(iter1, iter2);

  std::vector<int> b(length);

  std::copy(iter1, iter2, b.begin());

  for(unsigned int i = 0; i < b.size(); ++i)
    std::cout << b[i] << " ";

  std::cout << std::endl;

  print(0.0);

  std::cout << fabs(-9.9) << std::endl;
  
  std::complex<double> x(1,2);

  std::cout << abs(x) << std::endl;

  double * c;
  double y  = 10.8989;
  c = &y;
  
  std::cout << *c << std::endl;

  double & d = *c;
  d = 1.434;

  std::cout << d << " " << *c << std::endl; 


  for(unsigned int  i = 0; i < 5; ++i)
  {

    int z = i;
    std::cout << z << std::endl;
   
    for(unsigned int j = 0; j < 5; ++j)
    {

      int z = j;
      std::cout << z << " ";

    }

    std::cout << std::endl;

  }
  
 
}
     
