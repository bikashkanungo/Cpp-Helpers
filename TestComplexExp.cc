
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

typedef double Point[3]; 

int main()
{

  const std::complex<double> iota(0,1);
  double theta = rand()*1.0/RAND_MAX;

  //std::complex<double> exponential = exp(iota*M_PI);
  //std::complex<double> trig;
  //trig.real() = cos(theta);
  //trig.imag() = sin(theta);

  //std::cout << "Theta: " << theta << std::endl;
  //std::cout << exponential << " " << trig << " " << abs(exponential-trig) << std::endl;
  //const double pi = acos(-1);
     const double pi = std::acos(-1);
   const std::complex<double> i(0, 1);
 
   std::cout << std::fixed << " exp(i*pi) = " << std::exp(iota * theta) << '\n';
  return 0;

}

