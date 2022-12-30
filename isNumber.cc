#include <cstdlib>
#include <iostream>
#include <string>
#include <limits>
#include <cerrno>

bool isNumber(double &i, std::string s)
{
    char *end;
    double  d;
    errno = 0;
    d = strtod(s.c_str(), &end);
    if ((errno == ERANGE && d == std::numeric_limits<double>::max())) 
    {
      return false;
    }
  
    if ((errno == ERANGE && d == std::numeric_limits<double>::min())) 
    {
      return false;
    }
    
    if (s.size() == 0 || *end != 0) 
    {
      return false;
    }
  
    i = d;
    return true;
}

int main()
{

    std::string s;
    std::cout << "Enter string:";
    std::cin >> s;
    std::cout << std::endl;
    double i;
    bool flag = isNumber(i,s);
    if(flag)
        std::cout << "double: " << i << std::endl;
    else
        std::cout << "Not an double" << std::endl;

    return 0;
}
