#include<iostream>
#include<fstream>
#include <iomanip>
int main()
{

  std::ofstream file;
  file.open("Test");
  file << std::setprecision(16) << std::fixed;

  int a =1;
  double b= 3.0;
  file << a << b << std::endl;
  file.close();
}
