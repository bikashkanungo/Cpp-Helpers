#include <cstdlib>

#include <iostream>

#include <vector>

#include <functional>

#include <fenv.h>

#include <memory>

#include <map>

typedef std::vector<int> vecB;
typedef vecB::value_type value_type;

vecB vector(0);

class dummy {

 public:

   static const std::map<int, int> create_map()
   {
      std::map<int,int> m;
      m[1] = 2;
      m[3] = 4;
      m[5] = 6;
      return m;
    }

  static const std::map<int,int> myMap;
  static int i;

  virtual int func() {return 5;}



};

int dummy::i = 2;

const std::map<int, int> dummy::myMap =  dummy::create_map();


template <typename T, int I>
class Test{


};

 
template<typename T, int I> using vecA = Test<T, I>;

int main()
{

// BaseClass *basePointer = new BaseClass();
// int i = 3;
// auto x = i;

//  auto varA = basePointer->getValue();
//

// void * a;

// int b = 5;

// a = &(b);

// std::cout << "\n\n" << *((int*)a) << std::endl;

// std::shared_ptr<void> ptr (new int(5));

 //ptr = &(b);

// std::cout << *((int*)ptr) << std::endl;


//auto a = vector[0];


dummy dum;

auto a = dum.func();

auto b = dummy::myMap.find(1)->second;

std::cout << (int) b << std::endl;

int c = dummy::i;

vecA<int,c> vecX;

//vecA<decltype(a), c> vecConcrete = Test<int, 2>(); 


}




