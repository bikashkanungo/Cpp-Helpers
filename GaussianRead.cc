#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#define LENGTH 5.0

typedef double (*gaussianFunction)(const double , const double *, const double *);
typedef double (*func)(const double *, void *);

struct basis {
    int L; // Number of Gaussians used for contraction
    const double * alpha; // exponents of the L Gaussians used for contraction
    const double * c; // Coefficients of for each of the 
    const double * origin; // Origin or position of the atom about which the Gaussian is defined
    gaussianFunction gf;
    double normConst;
};

const double alpha_C6s [] = {3047.5249,457.36951,103.94869,29.210155,9.286663,3.163927};
const double c_C6s [] = {0.0018347,0.0140373,0.0688426,0.2321844,0.4679413,0.362312};

const double alpha_C3s [] = {7.8682724,1.8812885,0.5442493};
const double c_C3s [] = {-0.1193324,-0.1608542,1.1434564};

const double alpha_C3p [] = {7.8682724,1.8812885,0.5442493};
const double c_C3p [] = {0.0689991,0.316424,0.7443083};

const double alpha_C1s [] = {1.68714400E-01};
const double c_C1s [] = {1.0};

const double alpha_C1p [] = {1.68714400E-01};
const double c_C1p [] = {1.0};


const double alpha_H3s [] = {18.731137,2.8253937,0.6401217};
const double c_H3s [] = {0.0334946,0.23472695,0.81375733};

const double alpha [] = {1.61277800E-01};
const double c [] = {1.0};

//struct C_6s {
//    int size = 6;
//    double alpha [] = {3047.5249,457.36951,103.94869,29.210155,9.286663,3.163927};
//    double c []= {0.0018347,0.0140373,0.0688426,0.2321844,0.4679413,0.362312};
//}
//
//struct C_3s {
//    int size = 3;
//    double alpha [] = {7.8682724,1.8812885,0.5442493};
//    double c [] = {-0.1193324,-0.1608542,1.1434564};
//};
//
//struct C_3p {
//    int size = 3;
//    double alpha [] = {7.8682724,1.8812885,0.5442493};
//    double c [] = {0.0689991,0.316424,0.7443083};
//};
//
//struct C_1s {
//    int size = 1;
//    double alpha [] = {1.68714400E-01};
//    double c [] = {1.0};
//};
//
//struct C_1p {
//    int size = 1;
//    double alpha [] = {1.68714400E-01};
//    double c [] = {1.0};
//};
//
//struct H_3s {
//    int size = 3;
//    double alpha [] = {18.731137,2.8253937,0.6401217};
//    double c [] = {0.0334946,0.23472695,0.81375733};
//};
//
//struct H_1s {
//    int size = 1
//    double alpha [] = {1.61277800E-01};
//    double c [] = {1.0};
//};

double getDistance(const double * x, const double * y)
{
    double r = 0.0;
    for(unsigned int i = 0; i < 3; ++i)
        r += pow(x[i]-y[i],2.0);
    return sqrt(r);
}

double g1s(const double alpha, const double * x, const double * R)
{
    const double factor = 8*pow(alpha,3.0)/pow(M_PI,3.0);
    const double constant = pow(factor,1/4.0);
    const double r = getDistance(x,R);
    const double returnValue = constant*exp(-alpha*r*r);
    return returnValue;
}

double g2px(const double alpha, const double * x, const double * R)
{
    const double factor = 128*pow(alpha,5.0)/pow(M_PI,3.0);
    const double constant = pow(factor,1/4.0);
    const double r = getDistance(x,R);
    const double returnValue = constant*exp(-alpha*r*r);
    return (x[0]-R[0])*returnValue;
}

double g2py(const double alpha, const double * x, const double * R)
{
    const double factor = 128*pow(alpha,5.0)/pow(M_PI,3.0);
    const double constant = pow(factor,1/4.0);
    const double r = getDistance(x,R);
    const double returnValue = constant*exp(-alpha*r*r);
    return (x[1]-R[1])*returnValue;
}

double g2pz(const double alpha, const double * x, const double * R)
{
    const double factor = 128*pow(alpha,5.0)/pow(M_PI,3.0);
    const double constant = pow(factor,1/4.0);
    const double r = getDistance(x,R);
    const double returnValue = constant*exp(-alpha*r*r);
    return (x[2]-R[2])*returnValue;
}

double evaluateBasisValue(const basis * b, const double * x)
{
    const int L = b->L;
    const double normConst = b->normConst;
    double returnValue = 0.0;
    for(unsigned int i = 0; i < L; ++i)
    {
        const double alpha = b->alpha[i];
        const double c = b->c[i];
        const double * origin = b->origin;
        const gaussianFunction gf = b->gf;
        returnValue += c*((*gf)(alpha, x, origin));
    }
    return normConst*returnValue;
}

double trapezoidal3d(const double L,
                     const double h,
                     std::vector<const basis*> funcs)
{

    int N = L/h;
    double integral = 0.0;
    int numFuncs = funcs.size();
    for(int i = -N; i <= N; ++i)
    {
        std::vector<double> point(3);
        point[0] = i*h;
        double f2 = 0.0;
        for(int j = -N; j <= N; ++j)
        {
            point[1] = j*h;
            double f3 = 0.0;
            for(int k = -N; k <= N; ++k)
            {
                point[2] = k*h;
                double val = 1.0;
                for(unsigned int l = 0; l < numFuncs; ++l)
                    val *= evaluateBasisValue(funcs[l], &point[0]);
                if(k == -N || k == N)
                    f3 += 0.5*h*val;
                else 
                    f3 += h*val;
            }

            if(j == -N || j == N)
                f2 += 0.5*h*f3;
            else
                f2 += h*f3;
        }

        if(i == -N || i == N)
            integral += 0.5*h*f2;
        else
            integral += h*f2;

    }

    return integral;
                        
}

int main()
{
    double C [] = {0.0000000000,0.0000000000,0.0000000000};
    double H1 [] = {0.0000000000,0.0000000000,1.1089000000};
    double H2 [] = {1.0851085057,0.0000000000,-0.2284704375};

    std::vector<basis*> basisFunctions(0);

    double h;
    double len;
    std::cout <<"Enter length:";
    std::cin >> len;
    std::cout << std::endl;
    std::cout <<"Enter h:";
    std::cin >> h;
    std::cout << std::endl;
    //double N = len/h;
    //double integral = 0.0;
    //double alphaTest = 1.0;
    //for(int i = -N; i <= N; ++i)
    //{
    //    std::vector<double> point(3);
    //    point[0] = i*h;
    //    double f2 = 0.0;
    //    for(int j = -N; j <= N; ++j)
    //    {
    //        point[1] = j*h;
    //        double f3 = 0.0;
    //        for(int k = -N; k <= N; ++k)
    //        {
    //            point[2] = k*h;
    //            if(k == -N || k == N)
    //                f3 += 0.5*h*pow(point[0],2.0);//pow(g2px(alphaTest, &point[0], &C[0]),2.0);
    //            else
    //                f3 += h*pow(point[0],2.0);//pow(g2px(alphaTest, &point[0], &C[0]),2.0);
    //        }

    //        if(j == -N || j == N)
    //            f2 += 0.5*h*f3;
    //        else
    //            f2 += h*f3;
    //    }

    //    if(i == -N || i == N)
    //        integral += 0.5*h*f2;
    //    else
    //        integral += h*f2;
    //}

    //std::cout << "Integral Gaussian: " << integral << std::endl;

    // C_6s
    basis b1;
    b1.L = 6;
    b1.alpha = &alpha_C6s[0];
    b1.c = &c_C6s[0];
    b1.origin = &C[0];
    b1.gf = &g1s;
    b1.normConst = 1.0;

    // C_3p
    basis b3;
    b3.L = 3;
    b3.alpha = &alpha_C3p[0];
    b3.c = &c_C3p[0];
    b3.origin = &C[0];
    b3.gf = &g2px;
    b3.normConst = 1.0;


    // C_1s
    basis b6;
    b6.L = 1;
    b6.alpha = &alpha_C1s[0];
    b6.c = &c_C1s[0];
    b6.origin = &C[0];
    b6.gf = &g1s;
    b6.normConst = 1.0;

    // C_1p
    basis b7;
    b7.L = 1;
    b7.alpha = &alpha_C1p[0];
    b7.c = &c_C1p[0];
    b7.origin = &C[0];
    b7.gf = &g2px;
    b7.normConst = 1.0;
    
    // H_3s
    basis b10;
    b10.L = 3;
    b10.alpha = &alpha_H3s[0];
    b10.c = &c_H3s[0];
    b10.origin = &H1[0];
    b10.gf = &g1s;
    b10.normConst = 1.0;
    std::vector<const basis *> funcs(2);
    funcs[0] = &b1;
    funcs[1] = &b1;
    double integral;
    //integral = trapezoidal3d(len,h, funcs);
    //std::cout << "Integral Gaussian b1 and b1: " << integral << std::endl;
    //b1.normConst = 1.0/sqrt(integral);
    //
    //funcs[0] = &b10;
    //funcs[1] = &b10;
    //integral = trapezoidal3d(len,h, funcs);
    //std::cout << "Integral Gaussian b10 and b10: " << integral << std::endl;
    //b10.normConst = 1.0/sqrt(integral);

    funcs[0] = &b1;
    funcs[1] = &b10;
    integral = trapezoidal3d(len,h, funcs);
    std::cout << "Integral Gaussian b1 and b10: " << integral << std::endl;
    
    //funcs[1] = &b7;
    //integral = trapezoidal3d(len,h,funcs);
    //std::cout << "Integral 3 & 7: " << integral << std::endl;


}

