#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <complex>
#include <valarray>
#include <cassert>
#include <time.h>

int main()
{

  std::vector<std::complex<double> > v(10);
  for(unsigned int i = 0; i < v.size(); ++i)
  {
    v[i].real((rand()+0.0)/RAND_MAX);
    v[i].imag((rand()+0.0)/RAND_MAX);
  }

  double * re = &reinterpret_cast<double*>( v[0] )[0];
  double * im = &reinterpret_cast<double*>( v[0] )[1];

  std::cout << "\nPrinting vector as complex: " << std::endl;
  for(unsigned int i = 0; i < v.size(); ++i)
    std::cout << v[i] << "\t";
  
  std::cout << std::endl;
  std::cout << "\nPrinting vector as real and imag parts: " << std::endl;
  for(unsigned int i = 0; i < v.size(); ++i)
    std::cout << "(" << re[i] << "," << im[i] << ")\t";
}
