#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#define SIZE 4

void printList(int *list, int m)
{
  for(unsigned int i = 0; i < m; ++i)
    std::cout << list[i];
  std::cout <<  std::endl;
}
void permute(int * list, int n, int m)
{
  if(n==2)
  {
    printList(&list[0] - (m-n), m);
    int x = list[n-1];
    list[n-1] = list[n-2];
    list[n-2] = x;
    printList(&list[0] - (m-n), m);
  }
  else
  {

    int listCopy[n];
    for(unsigned int i = 0; i < n; ++i)
      listCopy[i] = list[i];
    for(unsigned int i = 0; i < n; ++i)
    {
      list[0] = listCopy[i];
      int count  = 1;
      for(unsigned int j = 0; j < n; ++j)
      {
	if(j!=i)
	{
	  list[count] = listCopy[j];
	  count ++;
	}
      }
      permute(&list[1], n-1, m);
    }
  }
}

int main()
{
  int list[SIZE];
  for(unsigned int i = 0; i < SIZE; ++i)
    list[i] = i+1;
  permute(list, SIZE, SIZE);
}

    

    
