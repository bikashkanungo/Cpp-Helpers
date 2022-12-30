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
#include <deque>

double cosSumAngles(std::deque<double> cosVal,
		    std::deque<double> sinVal);

double sinSumAngles(std::deque<double> cosVal,
		    std::deque<double> sinVal);


double sinSumAngles(std::deque<double> cosVal,
		    std::deque<double> sinVal)

{

  if(cosVal.size() == 2)
    return (sinVal[0]*cosVal[1] + cosVal[0]*sinVal[1]);

  else
  {

    double cosVal0 = cosVal[0];
    double sinVal0 = sinVal[0];
    cosVal.pop_front();
    sinVal.pop_front();
    return (sinVal0*cosSumAngles(cosVal,sinVal) + cosVal0*sinSumAngles(cosVal, sinVal));

  }

}

double cosSumAngles(std::deque<double> cosVal,
		    std::deque<double> sinVal)

{

  if(cosVal.size() == 2)
    return (cosVal[0]*cosVal[1] - sinVal[0]*sinVal[1]);

  else
  {

    double cosVal0 = cosVal[0];
    double sinVal0 = sinVal[0];
    cosVal.pop_front();
    sinVal.pop_front();
    return (cosVal0*cosSumAngles(cosVal,sinVal) - sinVal0*sinSumAngles(cosVal, sinVal));

  }

}

int
main()
{

  std::vector<std::deque<double> > cosVals(100, std::deque<double>(3));
  std::vector<std::deque<double> > sinVals(100, std::deque<double>(3));

  for(unsigned int i = 0; i < 100; ++i)
  {

    std::vector<double> thetas(3);
    double sumTheta = 0.0;
    for(unsigned int j = 0; j < 3; ++j)
    {

      thetas[j] = rand()*2*M_PI/RAND_MAX;
      if(rand()%2 == 0)
        thetas[j] = -thetas[j];
      cosVals[i][j] = cos(thetas[j]);
      sinVals[i][j] = sin(thetas[j]);
      sumTheta += thetas[j];

    }

    double cosSum = cosSumAngles(cosVals[i], sinVals[i]);
    double sinSum = sinSumAngles(cosVals[i], sinVals[i]);
 
    double cosSum1 = cos(sumTheta);
    double sinSum1 = sin(sumTheta);

    std::cout << "Cos: " << cosSum << " " << cosSum1 << " Sin: " << sinSum << " " << sinSum1 << std::endl;

  }

  std::vector<double> alphas(100);
  for(unsigned int i = 0; i < 100; ++i)
    alphas[i] = rand()*2*M_PI/RAND_MAX;

  std::vector<double> cosAlphas(100);

  std::transform(alphas.begin(), alphas.end(), cosAlphas.begin(), cos);
  for(unsigned int i = 0; i < 100; ++i)
    std::cout << cosAlphas[i] << " " << cos(alphas[i]) << std::endl;

}
  
