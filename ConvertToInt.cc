#include <cstdlib>

#include <iostream>

#include <string>


//
// A utility function to check whether x is numeric
//
bool isNumericChar(char x)
{

  return (x >= '0' && x <= '9')? true: false;

}
 
// A simple atoi() function. If the given string contains
// any invalid character, then this function throws error
int convertToInt(const char *str)
{

  if (str == NULL)
  {

    std::string message("String passed to convertToInt is empty or NULL");
    std::cout << message << std::endl;  
  
  }

 
  // 
  // Initialize result
  //
  int res= 0;       

  //
  // Initialize the sign
  //
  int sign = 1;

  //
  // Initialize the index of first character
  //
  int i = 0;

  //
  // remove white spaces
  //
  while(str[i]==' ')
   i++;

  // 
  // If number is negative, then update sign
  // and index of first digit
  //
  if (str[i] == '-')
  {
    sign = -1;

    i++; 

  }

  else if(str[i]=='+')
    i++;
 
  // Iterate through all digits of input string and update result
  for (; str[i] != '\0'; ++i)
  {

    if (isNumericChar(str[i]) == false)
    {

       std::string message("One of the character in the string passed to convertToInt is non-numeric");
       std::cout << message << std::endl;	
    
    }
    
    res = res*10 + str[i] - '0';

  }
 
  // Return result with sign
  return sign*res;

}


int main()
{

  std::string numberString;
  
  std::cout<< "Enter the number:";

  std::cin >> numberString;

  int number = convertToInt(numberString.c_str());

  std::cout << "Number: " << number << std::endl;
  std::cout << "Number: " << atoi(numberString.c_str()) << std::endl;

}
