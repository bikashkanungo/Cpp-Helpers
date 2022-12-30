#include <cstdlib>

#include <iostream>

#include <vector>

#include <functional>

#include <fenv.h>

#include <memory>

#include <map>

#include <complex>


typedef std::vector<double> DoubleVector;
typedef std::vector<std::complex<double> > ComplexVector;
typedef int QuadType;

QuadType getInt(const int i) {return i;}


template<typename T, QuadType quadType>
class QuadValuesManager {

  public:

};

template<typename T, QuadType quadType> using QuadParameters = std::vector< QuadValuesManager<T, quadType> >;


template<typename T, QuadType quadType>
QuadParameters<T, quadType>
makeQuadParameters(const int i, T &t, QuadType quadType)
{

  QuadParameters<T, quadType> quadParameters(i, QuadValuesManager<T,quadType>() );
  return quadParameters;

}


template<class Derived>
class Adapter {

  public:
  

};


class AdapterDerived : public Adapter<AdapterDerived> {

  public:


  template<QuadType quadType = 2, typename T = DoubleVector>
  QuadParameters<T, quadType>
  getQuadParameters();

};

template<QuadType quadType, typename T>
QuadParameters<T, quadType>
AdapterDerived::getQuadParameters()
{

  return QuadParameters<T, quadType>(1, QuadValuesManager<T, quadType>());

}

class FEMFunctional {

  public:

  template<class T>
  double
  getValue(Adapter<T> & adapter); 

  template<typename U, QuadType quadType>
  void
  getParametersAtQuadPoint(QuadParameters<U,quadType> & quadParameters); 

};

template<typename U, QuadType quadType>
void
FEMFunctional::getParametersAtQuadPoint(QuadParameters<U,quadType> & quadParameters)
{
  std::cout << "Voila" << std::endl;
}

template<class T>
double
FEMFunctional::getValue(Adapter<T> & adapter)
{

  T tObj; 

  auto parameters = tObj.getQuadParameters();
  getParametersAtQuadPoint(parameters);
  return 1.1;
}


int main()

{

  Adapter<AdapterDerived> adapter;
  
  FEMFunctional().getValue(adapter);


}



