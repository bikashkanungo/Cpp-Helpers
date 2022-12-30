#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

int main()
{
    int size = 9;
    double alphas[]={87141,
                     13052,
                     2970.6,
                     841.46,
                     274.59,
                     99.215,
                     38.694,
                     15.911,
                     6.7668};

    double cs[] = {2.82E-05,
                   0.00021819,
                   0.0011511,
                   0.0047859,
                   0.017134,
                   0.05112,
                   0.12805,
                   0.22773,
                   0.24758};

    const double alphaG = 1.4;
    const double normG = pow(2.0/M_PI,3.0/4.0)*pow(2.0,4.0)*pow(alphaG,11.0/4.0)/3.0;
    double integral = 0.0;
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
           const double alpha = alphas[i]+alphas[j];
           const double gint = pow(M_PI/alpha,3.0/2.0);
           const double c = cs[i]*cs[j];
           const double factor1 = 8*pow(alphas[i],3.0)/pow(M_PI,3.0);
           const double factor2 = 8*pow(alphas[j],3.0)/pow(M_PI,3.0);
           const double norm = pow(factor1*factor2,1.0/4.0);
           integral += c*norm*gint;
        }
    }
    std::cout << "Integral: " << integral << "\t SQRT: " << sqrt(integral) << std::endl;

    double integralSG = 0.0;
    for(unsigned int i = 0; i < size; ++i)
    {
        const double alphaCombined = alphas[i] + alphaG;
        const double factorS = 8*pow(alphas[i],3.0)/pow(M_PI,3.0);
        const double normS = pow(factorS,1.0/4.0);
        const double overlap = pow(M_PI,1.5)/(4.0*pow(alphaCombined,3.5));
        integralSG += normS*normG*cs[i]*overlap;
        std::cout << "[" << i << "] " << cs[i] << " " << normS << " " << normG << " " << overlap << std::endl; 
    }

    integralSG *= 1.0/sqrt(integral);

    std::cout <<"IntegralSG: " << integralSG << std::endl;
}
