#include<cstdlib>
#include<iostream>

int b[10];

void movePointer(int * x)
{

	for(unsigned int i = 0; i < 5; ++i)
		x++;
	std::cout << "x: " << x[0] << std::endl;
}

int main()
{

	int a[10];
	for(unsigned int i = 0; i < 10; ++i)
	{	a[i] = i;
		b[i] = i;
	}
	movePointer(&a[0]);
	movePointer(&b[0]);
	std::cout << "a: " << a[0] << std::endl;
	std::cout << "b: " << b[0] << std::endl;

}
