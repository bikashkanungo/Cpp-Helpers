#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

struct GaussianFunction
{
    double alpha;
    unsigned int numPolyTerms;
    std::vector<std::vector<int> > exponents;
    std::vector<double> polyCoeffs;
    double normConst;
};

double I(const int exponent, const double alpha)
{
    if(exponent == 0)
        return sqrt(M_PI/2.0)*pow(alpha,-0.5);
    else if(exponent == 2)
        return sqrt(M_PI/2.0)*pow(alpha,-1.5)/4.0;
    else if(exponent == 4)
        return sqrt(M_PI/2.0)*pow(alpha,-2.5)*3.0/16.0;
    else if(exponent == 6)
        return sqrt(M_PI/2.0)*pow(alpha,-3.5)*15.0/64.0;
    else if(exponent == 8)
        return sqrt(M_PI/2.0)*pow(alpha,-4.5)*105.0/256.0;
    else if(exponent == 10)
        return sqrt(M_PI/2.0)*pow(alpha,-5.5)*945.0/1024.0;
    else if((exponent %2) != 0)
        return 0.0;
    else
        std::cout << "Invalid argument" << std::endl;
}

double overlapGaussianFunctions(const GaussianFunction & gf1,
                                const GaussianFunction & gf2)
{

    unsigned int numPolyTerms1 = gf1.numPolyTerms;
    unsigned int numPolyTerms2 = gf2.numPolyTerms;
    const double alpha1 = gf1.alpha;
    const double alpha2 = gf2.alpha;
    const double alpha = 0.5*(alpha1+alpha2);

    double returnValue = 0.0;
    for(unsigned int i = 0; i < numPolyTerms1; ++i)
    {
        for(unsigned int j = 0; j < numPolyTerms2; ++j)
        {
            double componentIntegral = 1.0;
            for(unsigned int k = 0; k < 3; ++k)
            {
                const int exponent = gf1.exponents[i][k] + gf2.exponents[j][k];
                componentIntegral *= I(exponent,alpha);
            }
            returnValue += gf1.polyCoeffs[i]*gf2.polyCoeffs[j]*componentIntegral;
        }
    }

    return returnValue;
}


int main()
{

    unsigned int numFunctions;
    std::cout << "Enter number of gaussian functions:" ;
    std::cin >> numFunctions;
    std::vector<GaussianFunction*> gfVec(numFunctions);
    for(unsigned int i = 0; i < numFunctions; ++i)
    {

        std::cout << "\nEntering Values for Gaussian " << i << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        gfVec[i] = new GaussianFunction;
        GaussianFunction * gf = gfVec[i];
        std::cout << "Enter alpha: " ;
        std::cin >> gf->alpha;
        std::cout << std::endl;
        std::cout << "Enter number terms in the polynomial: " ;
        std::cin >> gf->numPolyTerms;
        std::cout << std::endl;
        gf->exponents.resize(gf->numPolyTerms,std::vector<int>(3,0));
        gf->polyCoeffs.resize(gf->numPolyTerms,0.0);
        gf->normConst = 1.0;
        for(unsigned int j = 0; j < gf->numPolyTerms; ++j) 
        {   
            std::cout << "Entering Values for term" << j << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "Enter x exponent: ";
            std::cin >> gf->exponents[j][0];
            std::cout << std::endl;
            std::cout << "Enter y exponent: ";
            std::cin >> gf->exponents[j][1];
            std::cout << std::endl;
            std::cout << "Enter z exponent: ";
            std::cin >> gf->exponents[j][2];
            std::cout << std::endl;
            std::cout << "Enter coeff: ";
            std::cin >> gf->polyCoeffs[j];
            std::cout << std::endl;
        }
            
    }

    for(unsigned int i = 0; i < numFunctions; ++i)
    {
        double overlap = overlapGaussianFunctions(*gfVec[i],*gfVec[i]);
        gfVec[i]->normConst = 1.0/sqrt(overlap);
        std::cout << "Overlap for Gaussian " << i << " is  " << overlap << std::endl; 
    }

    for(unsigned int i = 0; i < numFunctions; ++i)
    {
        for(unsigned int j = i+1; j < numFunctions; ++j)
        {

            double overlap = overlapGaussianFunctions(*gfVec[i],*gfVec[j]);
            overlap *= (gfVec[i]->normConst)*(gfVec[j]->normConst);
            std::cout << "Overlap for Gaussian pair (" << i << "," << j << "): " << overlap << std::endl; 

        }
    }
    
    for(unsigned int i = 0; i < numFunctions; ++i)
        delete gfVec[i];
    
//    int numTerms1, numTerms2;
//    double alpha1,alpha2;
//    std::cout << "Enter alpha for 1st polynomial: " ;
//    std::cin >> alpha1;
//    std::cout << std::endl;
//    std::cout << "Enter number terms in 1st polynomial: " ;
//    std::cin >> numTerms1;
//    std::cout << std::endl;
//    std::vector<std::vector<int> >exponents1(numTerms1,std::vector<int>(3,0));
//    std::vector<double> coeffs1(numTerms1,0.0);
//    std::cout << "\nEntering Values for 1st polynomial\n\n"  << std::endl;
//    for(unsigned int i = 0; i < numTerms1; ++i)
//    {   
//        std::cout << "Entering Values for term" << i << ":\n" << std::endl;
//        std::cout << "Term[" << i << "] Enter x exponent: ";
//        std::cin >> exponents1[i][0];
//        std::cout << std::endl;
//        std::cout << "Term[" << i << "] Enter y exponent: ";
//        std::cin >> exponents1[i][1];
//        std::cout << std::endl;
//        std::cout << "Term[" << i << "] Enter z exponent: ";
//        std::cin >> exponents1[i][2];
//        std::cout << std::endl;
//        std::cout << "Term[" << i << "] Enter coeff: ";
//        std::cin >> coeffs1[i];
//        std::cout << std::endl;
//    }
//    
//    std::cout << "Enter alpha for 2nd polynomial: " ;
//    std::cin >> alpha2;
//    std::cout << std::endl;
//    std::cout << "Enter number terms in 2nd polynomial: " ;
//    std::cin >> numTerms2;
//    std::cout << std::endl;
//    std::vector<std::vector<int> >exponents2(numTerms2,std::vector<int>(3,0));
//    std::vector<double> coeffs2(numTerms2,0.0);
//    std::cout << "\nEntering Values for 2nd polynomial\n\n"  << std::endl;
//    for(unsigned int i = 0; i < numTerms2; ++i)
//    {   
//        std::cout << "Entering Values for term" << i << ":\n" << std::endl;
//        std::cout << "Term[" << i << "] Enter x exponent: ";
//        std::cin >> exponents2[i][0];
//        std::cout << std::endl;
//        std::cout << "Term[" << i << "] Enter y exponent: ";
//        std::cin >> exponents2[i][1];
//        std::cout << std::endl;
//        std::cout << "Term[" << i << "] Enter z exponent: ";
//        std::cin >> exponents2[i][2];
//        std::cout << std::endl;
//        std::cout << "Term[" << i << "] Enter coeff: ";
//        std::cin >> coeffs2[i];
//        std::cout << std::endl;
//    }
//
//    double normConst1 = 0.0;
//    double normConst2 = 0.0;
//    for(unsigned int i = 0; i < numTerms1; ++i)
//    {
//        for(unsigned int j = 0; j < numTerms1; ++j)
//        {
//            double integral = 1.0;
//            for(unsigned int k = 0; k < 3; ++k)
//            {
//                int exponent = exponents1[i][k]+exponents1[j][k];
//                integral *= I(exponent,alpha1);
//            }
//            normConst1 += coeffs1[i]*coeffs1[j]*integral;
//        }
//    }
//
//    std::cout << "Overlap 1: " << normConst1 << std::endl;
//    normConst1 = sqrt(1.0/normConst1);
//
//    for(unsigned int i = 0; i < numTerms2; ++i)
//    {
//        for(unsigned int j = 0; j < numTerms2; ++j)
//        {
//            double integral = 1.0;
//            for(unsigned int k = 0; k < 3; ++k)
//            {
//                int exponent = exponents2[i][k]+exponents2[j][k];
//                integral *= I(exponent,alpha2);
//            }
//            normConst2 += coeffs2[i]*coeffs2[j]*integral;
//        }
//    }
//
//    std::cout << "Overlap 2: " << normConst2 << std::endl;
//    normConst2 = sqrt(1.0/normConst2);
//
//    double overlap = 0.0;
//    for(unsigned int i = 0; i < numTerms1; ++i)
//    {
//        for(unsigned int j = 0; j < numTerms2; ++j)
//        {
//            double integral = 1.0;
//            for(unsigned int k = 0; k < 3; ++k)
//            {
//                int exponent = exponents1[i][k]+exponents2[j][k];
//                double alpha = 0.5*(alpha1+alpha2);
//                integral *= I(exponent,alpha);
//            }
//            overlap += coeffs1[i]*coeffs2[j]*integral;
//        }
//    }
//
//    overlap *= normConst1*normConst2;
//    std::cout << "Overlap: " << overlap << std::endl;
}

