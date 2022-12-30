#include <cstdlib>

#include <iostream>

#include <vector>

#include <complex>

std::complex<double> 
testFunc1()
{

  return std::complex<double> (0.0);

}

template <typename T>
void
testFunc2()
{

  T a = testFunc1();

  std::cout << a << std::endl; 

}

int main()
{

  testFunc2<std::complex<double> >();
  
  return 0;

}

