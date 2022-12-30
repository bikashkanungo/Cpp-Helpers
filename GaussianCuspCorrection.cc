#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "mkl.h"

double
getGaussian(const double r, const double alpha)
{
    return exp(-alpha*r*r);
}

double
getGaussianFirstDerivative(const double r, const double alpha)
{
    return -2.0*alpha*r*exp(-alpha*r*r);
}

    double
getGaussianSecondDerivative(const double r, const double alpha)
{
    return exp(-alpha*r*r)*(4.0*alpha*alpha*r*r - 2.0*alpha);
}

    double
getP(const double r,
    const std::vector<double> & alphas)
{

    double returnValue = 0.0;
    for(unsigned int i = 0; i < alphas.size(); ++i)
        returnValue += alphas[i]*pow(r,i);
    return returnValue;
}

double
getCorrectedFunction(const double r,
    const std::vector<double> & alphas)
{
    const double exponent = getP(r,alphas);
    return exp(exponent);
}

std::vector<double>
getMatVec(const std::vector<double> & A,
          const std::vector<double> & v,
          const int M,
          const int N)
{

    std::vector<double> returnValue(M,0.0);
    for(unsigned int i = 0; i < M; ++i)
    {
        for(unsigned int j = 0; j < N; ++j)
            returnValue[i] += A[i*N + j]*v[j];
    }

    return returnValue;
}

double 
integralGaussian(const double rc,
                 const double alpha)
{

    int N = 100000;
    double h = rc/N;
    double returnValue = 0.0;
    for(unsigned int i = 1; i < N; ++i)
    {
        const double r = i*h;
        returnValue += r*r*getGaussian(r,alpha);
    }

    returnValue *= h;
    
    const double r0 = 0.0; 
    const double f0 = r0*r0*getGaussian(r0,alpha);
    const double r1 = rc; 
    const double f1 = r1*r1*getGaussian(r1,alpha);
    returnValue += 0.5*h*(f0+f1);
    
    return returnValue;

}
double 
integralCorrectedFunction(const double rc,
                          const std::vector<double> & alphas)
{

    int N = 100000;
    double h = rc/N;
    double returnValue = 0.0;
    for(unsigned int i = 1; i < N; ++i)
    {
        const double r = i*h;
        returnValue += r*r*getCorrectedFunction(r,alphas);
    }

    returnValue *= h;
    const double r0 = 0.0; 
    const double f0 = r0*r0*getCorrectedFunction(r0,alphas);
    const double r1 = rc; 
    const double f1 = r1*r1*getCorrectedFunction(r1,alphas);
    returnValue += 0.5*h*(f0+f1);
    return returnValue;

}

std::vector<double>
getMatInverse(const std::vector<double> & A,
              const int M)

{

    std::vector<double> ACopy(A);
    //
    // Lapack uses column-major format. So invert A^T 
    //
    int numRows = M;
    int lda = numRows;
    std::vector<int> ipiv(numRows);
    int info;
    dgetrf_(&numRows, 
            &numRows, 
            &ACopy[0], 
            &lda, 
            &ipiv[0], 
            &info);

    int lwork = numRows*numRows;
    std::vector<double> work(lwork);
    dgetri_(&numRows,
            &ACopy[0],
            &lda,
            &ipiv[0],
            &work[0],
            &lwork,
            &info);

    return ACopy;

}


int main()
{

    std::vector<double> alphas(5,0.0);
    alphas[0] = 0.0;
    const double Z = 1.0;
    double gAlpha;
    std::cout << "Enter gAlpha: ";
    std::cin >> gAlpha;
    std::cout << std::endl;
    
    double rc;
    std::cout << "Enter rc: ";
    std::cin >> rc;
    std::cout << std::endl;
    
    const double a = -2.0*Z*getGaussian(0.0,gAlpha);
    const double b = getGaussian(rc, gAlpha);
    const double c = getGaussianFirstDerivative(rc,gAlpha);
    const double d = getGaussianSecondDerivative(rc,gAlpha);
    const double e = integralGaussian(rc,gAlpha);
    const double g = integralCorrectedFunction(rc,alphas);
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "d: " << d << std::endl;
    std::cout << "e: " << e << std::endl;
    std::cout << "g: " << g << std::endl;

    double A[3][3];
    A[0][0] = pow(rc,2.0);  A[0][1] = pow(rc,3.0);          A[0][2] = pow(rc,4.0);
    A[1][0] = 2.0*rc;       A[1][1] = 3.0*pow(rc,2.0);      A[1][2] = 4.0*pow(rc,3.0);
    A[2][0] = 2.0;          A[2][1] = 6.0*rc;               A[2][2] = 12.0*pow(rc,2.0);
    
    std::vector<double> R(3*3);
    for(unsigned int i = 0; i < 3; ++i)
    {
        for(unsigned int j = 0; j < 3; ++j)
        {
            R[i*3+j] = A[i][j];
        }
    }

    std::vector<double> RInverse = getMatInverse(R, 3);
    for(unsigned int i = 0; i < 3; ++i)
    {
        for(unsigned int j = 0; j < 3; ++j)
        {
            std::cout << RInverse[i*3+j] << " ";
        }
        std::cout << std::endl;
    }
    
    std::vector<double> rhs(3);
    for(unsigned int i = 0; i < 1000; ++i)
    {
        alphas[1] = a/exp(alphas[0]);
        
        rhs[0] = log(b)-alphas[0]-alphas[1]*rc;
        rhs[1] = c/b-alphas[1];
        rhs[2] = d/b-pow(c/b,2.0);
        
        std::vector<double> y = getMatVec(RInverse,rhs, 3, 3);
        alphas[2] = y[0];
        alphas[3] = y[1];
        alphas[4] = y[2];

        const double g = integralCorrectedFunction(rc,alphas);
        
        std::cout << "[" << i << "]" << " g: " << g << " |e-g|: " << std::fabs(g-e) << " alphas: " << alphas[0] << " " << alphas[1] 
            << " " << alphas[2] << " " << alphas[3] << " " << alphas[4] << std::endl;
        
        //if(std::fabs(g-e)/e < 1e-10)
        //    break;
        //else
        //    alphas[0] -= (1.0 - e/g);
       
			 std::vector<double> alphasCopy(alphas);
			 alphasCopy[0] = 0.0;
       const double h = integralCorrectedFunction(rc,alphasCopy);
			 const double alpha0 = log(e/h);
			 if(std::fabs(alpha0-alphas[0]) < 1e-10)
				 break;
			else
				alphas[0] = alpha0;
    }

    std::string fileName = "GaussianAndCusp";
    std::ostringstream oss;
    oss << gAlpha;
    fileName += "_" + oss.str();
    oss.str("");
    oss << rc;
    fileName += "_" + oss.str();
    std::ofstream outfile;
    outfile.open(fileName.c_str());
    const int N = 10000;
    double h = 2.0*rc/N;
    for(unsigned int i = 0; i < N; ++i)
    {
        const double r = i*h;
        const double f1 = getGaussian(r,gAlpha);
        const double f2 = getCorrectedFunction(r,alphas);
        double f3 = 0.0;
        for(unsigned int j = 1; j < 5; ++j)
            f3 += j*alphas[j]*pow(r,j-1);
        f3 *= f2;
        outfile << r << "\t" << f1 << "\t" << f2 << "\t" << f3 << std::endl;
    }

    outfile.close();
}
