/**
 * This function evaluates the n-th order derivative of the function 
 * f(r)= r^l*e^(-alpha*r*r), where l is a non-negative integer and r>=0.
 */



#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>


// Returns value of Binomial Coefficient C(n, k)
int binomialCoeff(int n, int k)
{
  // Base Cases
  if (k==0 || k==n)
    return 1;
 
  // Recur
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k);
}

double derivativeRGaussian(const double r, const double alpha, const int l, const int derOrder)
{
    if(derOrder == 0 && l >= 0)
        return pow(r,l)*exp(-alpha*r*r);
    else if(derOrder == 0 && l < 0)
        return 0.0;
    else
        return l*derivativeRGaussian(r,alpha,l-1,derOrder-1) - 2*alpha*derivativeRGaussian(r,alpha,l+1,derOrder-1);
}

double functionVal(const double r, const double alpha, const int l)
{
    return pow(r,l)*exp(-alpha*r*r);
}

double forwardDifferenceDerivative(const double r, const double alpha, const int l, const double h, const int derOrder)
{
    double returnValue = 0.0;
    for(unsigned int i = 0; i <= derOrder; ++i)
    {

        double x = r + (derOrder-i)*h;
        returnValue += pow(-1,i)*binomialCoeff(derOrder,i)*functionVal(x,alpha,l);
    }       

    return returnValue/pow(h,derOrder);
}
int main()
{
    std::cout << std::setprecision(15);
    int l, n;
    double r, alpha;
    double h;
    std::cout << "Enter power (l) for r:";
    std::cin >> l;
    std::cout << std::endl;
    std::cout << "Enter order of derivative (n):";
    std::cin >> n;
    std::cout << std::endl;
    std::cout << "Enter value of r:";
    std::cin >> r;
    std::cout << std::endl;
    std::cout << "Enter value of alpha:";
    std::cin >> alpha;
    std::cout << std::endl;
    std::cout << "Enter value of h for finite difference:";
    std::cin >> h;
    std::cout << std::endl;

    double valExact = derivativeRGaussian(r,alpha,l,n);
    double valFDForward = forwardDifferenceDerivative(r,alpha,l,h,n);
    std::cout << "ValExact: " << valExact << "\t ValFDForward: " << valFDForward << std::endl; 
}

