#include <cstdlib>
#include <iostream>
#include <string>
#include <limits>
#include <cerrno>

bool isInteger(int &i, std::string s, const int base = 10)
{
    char *end;
    long  l;
    errno = 0;
    l = strtol(s.c_str(), &end, base);
    if ((errno == ERANGE && l == std::numeric_limits<long>::max()) || l > std::numeric_limits<int>::max()) 
    {
      return false;
    }
  
    if ((errno == ERANGE && l == std::numeric_limits<long>::min()) || l < std::numeric_limits<int>::min()) 
    {
      return false;
    }
    
    if (s.size() == 0 || *end != 0) 
    {
      return false;
    }
  
    i = l;
    return true;
}

int main()
{

    std::string s;
    std::cout << "Enter string:";
    std::cin >> s;
    std::cout << std::endl;
    int i;
    bool flag = isInteger(i,s);
    if(flag)
        std::cout << "Int: " << i << std::endl;
    else
        std::cout << "Not an Int" << std::endl;

    return 0;
}
