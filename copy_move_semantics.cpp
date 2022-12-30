#include <iostream>
#include <string>
#include <memory>
#include <vector>

/**
 * @note Some compilers use an optimization called copy elision which can bypass copy and move constructors/assigments.
 * In order to explicitly invoke copy and move constructors/assignments use the -fno-elide-constructors compiler flag
 */
 #include <iostream>
#include <string>
#include <memory>
#include <vector>

class A {
 
 public:
 
 A(int N, int init=0):
 d_N(N), 
 d_m(new int[N])
 {   
     std::cout << "A ctor" << std::endl;
     std::fill(&d_m[0], &d_m[d_N], init);
 }
 
 A(const A & a):d_N(a.d_N),
 d_m(new int[a.d_N])
 {
    std::cout << "A copy ctor" << std::endl;
    std::copy(&(a.d_m[0]), &(a.d_m[d_N]), &d_m[0]);
 }
 
 A(A && a) noexcept:
 d_N(std::move(a.d_N)),
 d_m(std::move(a.d_m))
 {     
   std::cout << "A move ctor" << std::endl;
 }
 
 
 A&
 operator=(const A& a)
 {
     std::cout << "A copy assignment" << std::endl;
     if(&a != this)
     {
       if(d_N != a.d_N)
       {
         delete d_m;
         d_N = a.d_N;
         d_m = new int[d_N];
       }
       std::copy(&(a.d_m[0]), &(a.d_m[d_N]), &d_m[0]);
     }
     return *this;
 }
 
 A&
 operator=(A&& a)
 {
     std::cout << "A move assignment" << std::endl;
     d_N=std::move(a.d_N);
     d_m=std::move(a.d_m); 
     return *this;
 }
 
  ~A(){delete d_m;}
  
 void
 print() 
 {
    for(int i = 0; i < d_N; ++i)
      std::cout << d_m[i] << " ";
      
    std::cout << std::endl;
 }
 
 int size() {return d_N;}
 
 private:
 int d_N;
 int * d_m;
    
};

class B{
    
 public:
   B(int N, int init=0): 
   d_A(std::make_shared<A>(N, init))
   {
       std::cout << "B ctor" << std::endl;
   }
   
   B(const B & b)
   :d_A(std::make_shared<A>((b.d_A)->size()))
   {
     *d_A=*(b.d_A); 
     std::cout << "B copy ctor" << std::endl;
   }
    
   B(B && b) noexcept:
   d_A(std::move(b.d_A))
   {
       std::cout << "B move ctor" << std::endl;
   }
   
   B&
   operator=(const B &b)
   {
       std::cout << "B copy assignment" << std::endl;
       d_A = std::make_shared<A>((b.d_A)->size());
       *d_A = *(b.d_A);
       return *this;
   }

   B&
   operator=(B &&b)
   {
       std::cout << "B move assignment" << std::endl;
       d_A = std::move(b.d_A);
       return *this;
   }
   
   void
   foo()
   {
     if(d_A != nullptr)
       d_A->print();
     else
       std::cout << "Attention! Accessing a null ptr." << std::endl;
   }

  private:
    std::shared_ptr<A> d_A;
};

int main()
{

    B b1(2,5);
    b1.foo();
    
    B b2(b1);
    b2.foo();
    
    B b3(B(3,4));
    b3.foo();
    
    B b4(0);
    b4 = b1;
    b4.foo();
    
    B b5(0);
    b5 = B(6,3);
    b5.foo();
    
    B b6 = std::move(b1);
    b6.foo();
    
    b1.foo();
    
    B b7(0);
    b7 = std::move(b2);
    b7.foo();
    
    b2.foo();
}

