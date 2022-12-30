#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

int main()
{
  int n;
  std::cout << "Enter integer: ";
  std::cin >> n;
  std::cout << std::endl;

  std::cout << "Factorial of " << n << " is " << factorial(n) << std::endl;

  return 0;

} 


