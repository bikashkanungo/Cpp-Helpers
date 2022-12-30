#include <iostream>
template<typename T, typename P>
bool pair_comparer(T a, P b) {
  // In real-world code, we wouldn't compare floating point values like
  // this. It would make sense to specialize this function for floating
  // point types to use approximate comparison.
  return a == b;
}

template<typename T, typename P, typename... Args>
bool pair_comparer(T a, P b, Args... args) {
  return pair_comparer(a,b) && pair_comparer(args...);
}

int main()
{
 
 std::cout << pair_comparer(1,1.0, 2.2, 2.2, 3.3, 3.3) << std::endl;
 
}
