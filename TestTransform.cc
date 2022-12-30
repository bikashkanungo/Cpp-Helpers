#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>
int main()
{

  std::vector<double> v1(100000000);
  std::fill(v1.begin(), v1.end(), 1.0);
  std::vector<double> v2;

  std::clock_t timeBegin = clock();
  //v2.reserve(v1.size());
  //v2.insert(v2.end(),v1.begin(), v1.end());

  v2.resize(v1.size());
  std::copy(v1.begin(), v1.end(), v2.begin());
  std::clock_t timeEnd = clock();

  double elapsedTime = double(timeEnd - timeBegin)/CLOCKS_PER_SEC/60.0;

  std::cout << "Time taken for insert: " << elapsedTime << std::endl;

  return 0;

  const int length = 10, numVecs = 5;
  std::vector<double> a(length);
  std::vector<double> b(length);
  std::vector<double> c(length);

  for(unsigned int i = 0; i < length; ++i)
  {

    a[i] = 1.0*rand()/RAND_MAX;
    b[i] = 1.0*rand()/RAND_MAX;
    c[i] = a[i]*b[i];
  
  }

  std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::multiplies<double>());

  std::cout << std::endl;
  for(unsigned int i = 0; i < length; ++i)
    std::cout << a[i] << " " << c[i] << " " << a[i] - c[i] << std::endl;

  std::vector<double> X(length*numVecs);
  std::vector<double> Y(length*numVecs);
  std::vector<double> Z(length*numVecs);

  for(unsigned int  i = 0; i < length; ++i)
  {

    for(unsigned int j = 0; j < numVecs; ++j)
    {

      X[i*numVecs + j] = 1.0*rand()/RAND_MAX;
      Y[i*numVecs + j] = a[i]*X[i*numVecs + j];

    }

  }

  std::vector<double>::iterator iter1 = X.begin();
  std::vector<double>::iterator iter3 = Z.begin();
  for(unsigned int i = 0; i < length; ++i)
  {

    std::vector<double>::iterator iter2 = iter1 + numVecs;
    std::transform(iter1, iter2, iter3, std::bind1st(std::multiplies<double>(),a[i]));
    
    iter1 = iter2;

    iter3 += numVecs;

  }

  for(unsigned int  i = 0; i < length; ++i)
  {

    for(unsigned int j = 0; j < numVecs; ++j)
    {
 
       std::cout << "("<< Y[i*numVecs + j] << "," << Z[i*numVecs + j] << "," << Y[i*numVecs + j] - Z[i*numVecs + j] << ")\t";

    }

    std::cout << std::endl;

  }
     

  std::vector<double> d(length);
  std::vector<double> e(length);

  for(unsigned int i = 0; i < length; ++i)
    d[i] = a[i] - b[i];

  std::transform(a.begin(), a.end(), b.begin(), e.begin(), std::minus<double>());
   
  std::cout << std::endl;
  for(unsigned int i = 0; i < length; ++i)
    std::cout << d[i] << " " << e[i] << " " << d[i] - e[i] << std::endl; 

  return 0;

} 
