#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>


void foo(std::ofstream & file)
{

  file << "World" << std::endl;

}

int main()
{

   
  std::ofstream file;
  std::string name("testFile");
  file.open(name.c_str());

  file << "Hello" << std::endl;
  foo(file);

}
