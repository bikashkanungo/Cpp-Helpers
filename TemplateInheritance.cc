#include <cstdlib>
#include <iostream>
#include "TemplateBase.h"

//template <typename T>
//class Base {
//
//  public:
//
//  Base()
//  {
//
//    std::cout << "Base Constructor" << std::endl;
//
//  }
//
//  ~Base() {}
//
//  virtual
//  void
//  setValue(T val)
//  {
//
//    d_val = val;
//
//  }
//
//  virtual
//  T
//  getValue()
//  {
//
//    return d_val;
//
//  }
//
//  virtual void foo() = 0;
//
//
//  private:
//
//  T d_val;
//
//};

template <typename T>
class Derived: public Base<T> 
{

  public:
  Derived(){}
  ~Derived(){}
  void foo() { std::cout << "Derived foo" << std::endl; }
  
};

template <typename T>
class A {

  A(Base<T> & base):
  d_baseObj(base)
  {

    std::cout << "A constructor" << std::endl;
	
  }
 
  ~A(){}
  private:
  Base<T> d_base;

};


int main()
{

  Base<int> * basePtr;
  basePtr = new Derived<int>();

  basePtr->foo();

}
