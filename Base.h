#include <vector>

namespace test {

class BaseClass {

  public:

  template <typename T>
  std::vector<T> makeT(T value)
  {

    return std::vector<T>(2,value);

  }


   
  std::vector<int> getValue()
  { return makeT(3); }

 // std::vector<double> getValue()
  //{ return makeT(3.33);}
    
};

}
