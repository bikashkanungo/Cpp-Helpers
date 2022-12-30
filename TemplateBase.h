#include <cstdlib>
#include <iostream>

template <typename T>
class Base {

  public:

  Base();
  ~Base();

  virtual void setValue(T val);
  virtual T getValue();
  virtual void foo() = 0;

  private:

  T d_val;

};
#include "TemplateBase.cc"
