#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#define TOL 1e-12

double f1(const double x)
{
  if(x <= 0.0)
    return 0.0;
  
  else
    return exp(-1.0/x);

}

double f1Der(const double x)
{

  return f1(x)/(x*x);

}

double f2(const double x)
{

  return (f1(x)/(f1(x)+f1(1-x)));

}


double f2Der(const double x)
{

  if(fabs(x-0.0) < TOL || fabs(1-x) < TOL)
    return 0.0;
	
  else
  return ((f1Der(x)*f1(1-x) + f1(x)*f1Der(1-x))/(pow(f1(x) + f1(1-x), 2.0)));

}

double Y(const double x, const double r)
{

  return  (1 - .5*(x-r)/r);

}

double YDer(const double x, const double r)
{

  return (-.5/r);

}


double getMollifierValue(const double x, const double r)
{

   const double y = Y(x,r);
   return pow(f2(y), 1.0);

}

double getMollifierDerivative(const double x, const double r)
{

  const double y = Y(x,r);
  return f2Der(y)*YDer(x,r);
   

}


int main()
{

   //std::cout << getMollifierDerivative(3.0,1.0) << std::endl;

int polyOrder = 3;
  for(int i = polyOrder; i >= 0; i--)
  {
      for(int j = polyOrder; j >= 0; j--)
      {
          for(int k = polyOrder; k >= 0; k--)
          {
              if(i+j+k == polyOrder)
              {
                  std::vector<int> exponents(3);
                  exponents[0] = i;
                  exponents[1] = j;
                  exponents[2] = k;
                  std::cout  << i << " " << j << " " << k << std::endl;
              }
        }
    }
    
  }

   return 0;

}

