#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


class A {

  public:

  int a;

  A() { a = 5; }

  ~A(){}

  void
  get(int A::* ptr)
  {
    
    ptr = &A::a;

  }

};

int main()
{

  A objA;

  int A::* ptr;

  objA.get(ptr);

  std::cout << "Value: " << objA.*ptr << std::endl;

}
