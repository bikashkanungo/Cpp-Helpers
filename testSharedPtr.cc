#include <iostream>
#include<memory>
#include<vector>

using namespace std;

class A 
{
 
 public:
 A() {}
 ~A() {}
  int foo()
  {
      d = 10;
      return d;
  }
  
 private:
  int d;
    
};

int main()
{
    
    cout<<"Hello World" << endl;
    vector<shared_ptr<A> > x(10, make_shared<A>());
	
    vector<A> y(10);
    
    vector<shared_ptr<A> >::iterator it = x.begin();
    vector<A>::iterator it1 = y.begin();
    bool a = true;
    a = (x.begin() == x.end());
    cout << a;
    for(; it != x.end(); ++it)
    {
        int p = (*it)->foo();
        cout << p << endl;
    }
    
    return 0;
}
