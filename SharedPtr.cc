#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <cerrno>
#include <cassert>
#include <set>
#include <memory>


class A {

	public:
	
	A(std::shared_ptr<const std::vector<std::vector<double> > > ptr):
	d_ptr(ptr)
	{

	}

	~A()
	{
	
	}

void
	someFunction()
	{
			const std::vector<std::vector<double> > & d = *d_ptr;
			int M = d.size();//(*d_ptr).size();
			int N = d[0].size(); //(*d_ptr)[0].size();
			for(unsigned int i = 0; i < M; ++i)
			{
				for(unsigned int j = 0; j < N; ++j)
					std::cout << d[i][j] << " ";
				std::cout << std::endl;
			}
	}

	private:
	std::shared_ptr<const std::vector<std::vector<double> > > d_ptr;

};


void callLocal(std::shared_ptr<const std::vector<std::vector<double> > > ptr)
{

	A a(ptr);
	std::cout << "Count : " << ptr.use_count() << std::endl;
	a.someFunction();
	//delete &a;
	std::cout << "Count : " << ptr.use_count() << std::endl;


}

int main()
{

	std::vector<std::vector<double> > v(10,std::vector<double>(3,1.1));
	std::shared_ptr<const std::vector<std::vector<double> > > ptr = 
		std::make_shared<const std::vector<std::vector<double> > >(std::move(v));
	
	callLocal(ptr);
	
	std::cout << "Count : " << ptr.use_count() << std::endl;
	
}
