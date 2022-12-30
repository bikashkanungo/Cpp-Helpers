#include <cstdlib>

#include <iostream>

#include <string>

#include <vector>

void initX(double x[3], const int i)
{

  x[0] = 3*i + 1;
  x[1] = 3*i + 2;
  x[2] = 3*i + 3;

}

int main()
{

  const int polyOrder = 2;
  const int higherOrderCornerIds[8]  = {0, 
                                       polyOrder*(polyOrder+1)*(polyOrder+1), 
                                       polyOrder*(polyOrder+1)*(polyOrder+2), 
                                       polyOrder*(polyOrder+1), 
                                       polyOrder, 
                                       polyOrder*(polyOrder*polyOrder + 2*polyOrder + 2), 
                                       (polyOrder+1)*(polyOrder+1)*(polyOrder+1) - 1,
                                       polyOrder*(polyOrder+2)};

  for(unsigned int  i = 0; i < 8; ++i)
    std::cout <<  higherOrderCornerIds[i] << " " << std::endl;

  std::vector<std::vector<double> > x(10,std::vector<double>(3));

  for(unsigned int i = 0; i < x.size(); ++i)
  {   

    initX(&x[i][0], i);
    for(unsigned int j = 0; j < 3; ++j)
      std::cout << x[i][j] << " ";
  
    std::cout << std::endl;

  }

 

} 
