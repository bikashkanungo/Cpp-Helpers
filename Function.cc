#include <iostream>
#include <cmath>
#include <vector>

namespace mathConstants
{

  constexpr double e = std::exp(1.0);
}

template<typename T, typename Q>
class Function {

  public:
	virtual Q operator()(const T & t) = 0;
	virtual std::vector<Q> operator()(const std::vector<T> & t) = 0;

};

typedef Function<double,double> Scalar1DFunction;

class LogModX : public Scalar1DFunction
{

  public:
	LogModX(const double base = mathConstants::e):
	  d_base(base) {}

	double operator()(const double & x) override
	{

	  return std::log(fabs(x))/std::log(d_base);
	}
	
	std::vector<double> operator()(const std::vector<double> & x) override
	{

	  std::vector<double> returnValue(x.size(),0.0);
	  const double logBase = std::log(d_base);
	  for(unsigned int i = 0; i < x.size(); ++i)
		returnValue[i] = std::log(fabs(x[i]))/logBase;
	  return returnValue; 
	}

  private:
	double d_base;

};


int main()
{

  Function<double,double> * f = new LogModX();
  std::vector<double> x(3,4.0);
  std::cout << (*f)(-2.0) << std::endl;
  std::vector<double> y = (*f)(x);
  for(unsigned int i = 0; i < y.size(); ++i)
	std::cout << y[i] << " ";
  std::cout << std::endl;
}
