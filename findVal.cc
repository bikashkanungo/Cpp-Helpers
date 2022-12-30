#include <cstdlib>
#include <iostream>
#include <cmath>

#define EPS_MACHINE 1e-16
#define EXPONENT_TOL log(EPS_MACHINE)
int main()
{
    double a = 1.68714400E-01;
    double b = 1.61277800E-01;
    double p = a+b;
    double r = 1.1089;
    double e = exp(-a*b*r*r/p);
    double f = pow(2*a*b/(M_PI*p),3.0/4.0);
    double kab = f*e;
    std::cout << "Exponent tolerance: " << EXPONENT_TOL << std::endl;
    std::cout << "Kab: " << kab << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout <<"prefactor: " << pow(2.0*p/M_PI,0.75) << std::endl; 
}
