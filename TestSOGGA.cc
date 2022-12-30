#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <cerrno>
#include <cassert>
#include <set>
#include <xc.h>

double ax_global[] ={0.50000, -2.95535, 15.7974, -91.1804, 96.2030, 0.18683};
double bx_global[] = {0.50000, 3.50743, -12.9523, 49.7870, -33.2545, -11.1396};
double Ax_global = -3.0/4.0*std::pow(3.0/(M_PI),1.0/3.0);
//double Ax_global = -3.0/2.0*std::pow(3.0/(4.0*M_PI),1.0/3.0);
double As_global = 1.0/std::pow(24.0*M_PI*M_PI,1.0/3.0);
//double As_global = 1.0/std::pow(48.0*M_PI*M_PI,1.0/3.0);
double mu_global = 0.1234567901234567901234567901234567901235;
double kappa_global = 0.552;


double dotProduct(const std::vector<double> & x,
        const std::vector<double> & y)
{

    double returnValue = 0.0;
    for(unsigned int i = 0; i < x.size(); ++i)
        returnValue += x[i]*y[i];

    return returnValue;
}
    
std::vector<double>
getFDCoeffsForwardDiff(const int numFDPoints)
{

    std::vector<double> returnValue(numFDPoints);
    if(numFDPoints == 2)
    {

        returnValue[0] = -1.0;
        returnValue[1] = 1.0;

    }

    else if(numFDPoints == 3)
    {

        returnValue[0] = -3.0/2.0;
        returnValue[1] = 4.0/2.0;
        returnValue[2] = -1.0/2.0; 	

    }

    else if(numFDPoints == 4)
    {

        returnValue[0] = -11.0/6.0;
        returnValue[1] = 18.0/6.0;
        returnValue[2] = -9.0/6.0;
        returnValue[3] = 2.0/6.0;

    }

    else if(numFDPoints == 5)
    {

        returnValue[0] = -25.0/12.0;
        returnValue[1] = 48.0/12.0;
        returnValue[2] = -36.0/12.0;
        returnValue[3] = 16.0/12.0;
        returnValue[4] = -3.0/12.0; 	

    }

    else if(numFDPoints == 6)
    {

        returnValue[0] = -137.0/60.0;
        returnValue[1] = 300.0/60.0;
        returnValue[2] = -300.0/60.0;
        returnValue[3] = 200.0/60.0;
        returnValue[4] = -75.0/60.0; 	
        returnValue[5] = 12.0/60.0; 	

    }
    
    else if(numFDPoints == 7)
    {

        returnValue[0] = -147.0/60.0;
        returnValue[1] = 360.0/60.0;
        returnValue[2] = -450.0/60.0;
        returnValue[3] = 400.0/60.0;
        returnValue[4] = -225.0/60.0; 	
        returnValue[5] = 72.0/60.0; 	
        returnValue[6] = -10.0/60.0; 	

    }
    
    else if(numFDPoints == 8)
    {

        returnValue[0] = -1089.0/420.0;
        returnValue[1] = 2940.0/420.0;
        returnValue[2] = -4410.0/420.0;
        returnValue[3] = 4900.0/420.0;
        returnValue[4] = -3675.0/420.0; 	
        returnValue[5] = 1764.0/420.0; 	
        returnValue[6] = -490.0/420.0; 	
        returnValue[7] = 60.0/420.0; 	

    }
    
    else if(numFDPoints == 9)
    {

        returnValue[0] = -2283.0/840.0;
        returnValue[1] = 6720.0/840.0;
        returnValue[2] = -11760.0/840.0;
        returnValue[3] = 15680.0/840.0;
        returnValue[4] = -14700.0/840.0; 	
        returnValue[5] = 9408.0/840.0; 	
        returnValue[6] = -3920.0/840.0; 	
        returnValue[7] = 960.0/840.0; 	
        returnValue[8] = -105.0/840.0; 	

    }
    
    else if(numFDPoints == 10)
    {
        returnValue[0] = -7129.0/2520.0;
        returnValue[1] = 22680.0/2520.0;
        returnValue[2] = -45360.0/2520.0;
        returnValue[3] = 70560.0/2520.0;
        returnValue[4] = -79380.0/2520.0; 	
        returnValue[5] = 63504.0/2520.0; 	
        returnValue[6] = -35280.0/2520.0; 	
        returnValue[7] = 12960.0/2520.0; 	
        returnValue[8] = -2835.0/2520.0; 	
        returnValue[9] = 280.0/2520.0; 	

    }
    
    else if(numFDPoints == 11)
    {
        returnValue[0] = -7381.0/2520.0;
        returnValue[1] = 25200.0/2520.0;
        returnValue[2] = -56700.0/2520.0;
        returnValue[3] = 100800.0/2520.0;
        returnValue[4] = -132300.0/2520.0; 	
        returnValue[5] = 127008.0/2520.0; 	
        returnValue[6] = -88200.0/2520.0; 	
        returnValue[7] = 43200.0/2520.0; 	
        returnValue[8] = -14175.0/2520.0; 	
        returnValue[9] = 2800.0/2520.0; 	
        returnValue[10] = -252/2520.0; 	
    }

    else if(numFDPoints == 12)
    {
        returnValue[0] = -83711.0/27720.0;
        returnValue[1] = 304920.0/27720.0;
        returnValue[2] = -762300.0/27720.0;
        returnValue[3] = 1524600.0/27720.0;
        returnValue[4] = -2286900.0/27720.0; 	
        returnValue[5] = 2561328.0/27720.0; 	
        returnValue[6] = -2134440.0/27720.0; 	
        returnValue[7] = 1306800.0/27720.0; 	
        returnValue[8] = -571725.0/27720.0; 	
        returnValue[9] = 169400.0/27720.0; 	
        returnValue[10] = -30492.0/27720.0; 	
        returnValue[11] = 2520.0/27720.0; 	
    }
    
    else if(numFDPoints == 13)
    {
        returnValue[0] = -86021.0/27720.0;
        returnValue[1] = 332640.0/27720.0;
        returnValue[2] = -914760.0/27720.0;
        returnValue[3] = 2032800.0/27720.0;
        returnValue[4] = -3430350.0/27720.0; 	
        returnValue[5] = 4390848.0/27720.0; 	
        returnValue[6] = -4268880.0/27720.0; 	
        returnValue[7] = 3136320.0/27720.0; 	
        returnValue[8] = -1715175.0/27720.0; 	
        returnValue[9] = 677600.0/27720.0; 	
        returnValue[10] = -182952.0/27720.0; 	
        returnValue[11] = 30240.0/27720.0; 	
        returnValue[12] = -2310.0/27720.0; 	
    }
    
    else if(numFDPoints == 14)
    {
        returnValue[0] = -1144013/360360.0;
        returnValue[1] = 4656960.0/360360.0;
        returnValue[2] = -13873860.0/360360.0;
        returnValue[3] = 33633600.0/360360.0;
        returnValue[4] = -62432370.0/360360.0; 	
        returnValue[5] = 88792704.0/360360.0; 	
        returnValue[6] = -97117020.0/360360.0; 	
        returnValue[7] = 81544320.0/360360.0; 	
        returnValue[8] = -52026975.0/360360.0; 	
        returnValue[9] = 24664640.0/360360.0; 	
        returnValue[10] = -8324316.0/360360.0; 	
        returnValue[11] = 1834560.0/360360.0; 	
        returnValue[12] = -210210.0/360360.0; 	
        returnValue[13] = 1980.0/360360.0; 	
    }
    else
    {

        std::cout << "Invalid number of FD points. Please enter number of FD points between 2-14." << std::endl;
    }

    return returnValue;
}

    std::vector<double>
getFDCoeffs(const int numFDPoints)
{

    std::vector<double> returnValue(numFDPoints);
    if(numFDPoints == 3)
    {

        returnValue[0] = -1.0/2.0;
        returnValue[1] = 0.0;
        returnValue[2] = 1.0/2.0;

    }

    else if(numFDPoints == 5)
    {

        returnValue[0] = 1.0/12.0;
        returnValue[1] = -2.0/3.0;
        returnValue[2] = 0.0;
        returnValue[3] = 2.0/3.0;
        returnValue[4] = -1.0/12.0; 	

    }

    else if(numFDPoints == 7)
    {

        returnValue[0] = -1.0/60.0;
        returnValue[1] = 3.0/20.0;
        returnValue[2] = -3.0/4.0;
        returnValue[3] = 0.0;
        returnValue[4] = 3.0/4.0; 	
        returnValue[5] = -3.0/20.0;
        returnValue[6] = 1.0/60.0; 	

    }

    else if(numFDPoints == 9)
    {

        returnValue[0] = 1.0/280.0;
        returnValue[1] = -4.0/105.0;
        returnValue[2] = 1.0/5.0;
        returnValue[3] = -4.0/5.0;
        returnValue[4] = 0.0; 	
        returnValue[5] = 4.0/5.0;
        returnValue[6] = -1.0/5.0; 	
        returnValue[7] = 4.0/105.0; 	
        returnValue[8] = -1.0/280.0; 	

    }

    else if(numFDPoints == 11)
    {

        returnValue[0] = -2.0/2520.0;
        returnValue[1] = 25.0/2520.0;
        returnValue[2] = -150.0/2520.0;
        returnValue[3] = 600.0/2520.0;
        returnValue[4] = -2100.0/2520.0; 	
        returnValue[5] = 0.0/2520.0;
        returnValue[6] = 2100.0/2520.0; 	
        returnValue[7] = -600.0/2520.0; 	
        returnValue[8] = 150.0/2520.0; 	
        returnValue[9] = -25.0/2520.0; 	
        returnValue[10] = 2.0/2520.0; 	

    }
    
    else if(numFDPoints == 13)
    {

        returnValue[0] = 5.0/27720.0;
        returnValue[1] = -72.0/27720.0;
        returnValue[2] = 495.0/27720.0;
        returnValue[3] = -2200.0/27720.0;
        returnValue[4] = 7425.0/27720.0; 	
        returnValue[5] = -23760/27720.0;
        returnValue[6] = 0.0/27720.0; 	
        returnValue[7] = 23760.0/27720.0; 	
        returnValue[8] = -7425.0/27720.0; 	
        returnValue[9] = 2200.0/27720.0; 	
        returnValue[10] = -495.0/27720.0; 	
        returnValue[11] = 72.0/27720.0; 	
        returnValue[12] = -5.0/27720.0; 	

    }

    else
    {

        std::cout << "Invalid number of FD points. Please enter number of FD points as 3, 5, 7, 9, 11 or 13." << std::endl;
    }

    return returnValue;
}

    double
getHigherOrderDerivativeFD(const std::vector<double> & dependentVariables,
                           double (*functionPtr)(const std::vector<double> &),                    
                           std::vector<int> & indices,
                           const int numFDPoints,
                           const double h)
{

    int numIndices = indices.size();
    double returnValue = 0.0;
    if(numIndices == 0)
    {
        returnValue = functionPtr(dependentVariables);
        return returnValue;
    }

    else
    {
        int currentIndex = indices[numIndices-1];
        std::vector<int> indicesNext(numIndices-1);
        for(unsigned int i = 0; i < numIndices - 1; ++i)
            indicesNext[i] = indices[i];

        std::vector<double> coeffs = getFDCoeffs(numFDPoints);
        for(unsigned int iPoint = 0; iPoint < numFDPoints; ++iPoint)
        {
            std::vector<double> FDPoint(dependentVariables);
            int shiftIndex = (int)(iPoint-numFDPoints/2);
            FDPoint[currentIndex] += shiftIndex*h;
            const double factor = coeffs[iPoint]/h;
            returnValue += factor*getHigherOrderDerivativeFD(FDPoint,
                                                             functionPtr,
                                                             indicesNext,
                                                             numFDPoints,
                                                             h);
        }

        return returnValue;
    }
}
    
  double
getHigherOrderDerivativeFDForwardDiff(const std::vector<double> & dependentVariables,
                           double (*functionPtr)(const std::vector<double> &),                    
                           std::vector<int> & indices,
                           const int numFDPoints,
                           const double h)
{

    int numIndices = indices.size();
    double returnValue = 0.0;
    if(numIndices == 0)
    {
        returnValue = functionPtr(dependentVariables);
        return returnValue;
    }

    else
    {
        int currentIndex = indices[numIndices-1];
        std::vector<int> indicesNext(numIndices-1);
        for(unsigned int i = 0; i < numIndices - 1; ++i)
            indicesNext[i] = indices[i];

        std::vector<double> coeffs = getFDCoeffsForwardDiff(numFDPoints);
        for(unsigned int iPoint = 0; iPoint < numFDPoints; ++iPoint)
        {
            std::vector<double> FDPoint(dependentVariables);
            int shiftIndex = (int)(iPoint);
            FDPoint[currentIndex] += shiftIndex*h;
            const double factor = coeffs[iPoint]/h;
            returnValue += factor*getHigherOrderDerivativeFDForwardDiff(FDPoint,
                                                             functionPtr,
                                                             indicesNext,
                                                             numFDPoints,
                                                             h);
        }

        return returnValue;
    }
}


  double
getHigherOrderDerivativeFDForwardDiffSelective(const std::vector<double> & dependentVariables,
                           const int index,
                           double (*functionPtr)(const std::vector<double> &, const int),                    
                           std::vector<int> & indices,
                           const int numFDPoints,
                           const double h)
{

    int numIndices = indices.size();
    double returnValue = 0.0;
    if(numIndices == 0)
    {
        returnValue = functionPtr(dependentVariables, index);
        return returnValue;
    }

    else
    {
        int currentIndex = indices[numIndices-1];
        std::vector<int> indicesNext(numIndices-1);
        for(unsigned int i = 0; i < numIndices - 1; ++i)
            indicesNext[i] = indices[i];

        std::vector<double> coeffs = getFDCoeffsForwardDiff(numFDPoints);
        for(unsigned int iPoint = 0; iPoint < numFDPoints; ++iPoint)
        {
            std::vector<double> FDPoint(dependentVariables);
            int shiftIndex = (int)(iPoint);
            FDPoint[currentIndex] += shiftIndex*h;
            const double factor = coeffs[iPoint]/h;
            returnValue += factor*getHigherOrderDerivativeFDForwardDiffSelective(FDPoint,
                                                             index,
                                                             functionPtr,
                                                             indicesNext,
                                                             numFDPoints,
                                                             h);
        }

        return returnValue;
    }
}
double
getRho(const std::vector<double> & point)
{
    const double x = point[0];
    const double y = point[1];
    const double z = point[2];
    const double x2 = x*x;
    const double y2 = y*y;
    const double z2 = z*z;
    const double r = std::sqrt(x2 + y2 + z2);
    return exp(-r) + x2*exp(-r) + x2*y*exp(-r) + x2*y2*exp(-r);
}

std::vector<double>
getGradRho(const std::vector<double> & point)
{

    std::vector<double> returnValue(3,0.0);
    const double h = 1.0e-6;
    const int numFDPoints = 13;
    for(unsigned int i = 0; i < 3; ++i)
    {
        std::vector<int> indices(1);
        indices[0] = i;
        returnValue[i] = getHigherOrderDerivativeFD(point,
                                                    &getRho,
                                                    indices,
                                                    numFDPoints,
                                                    h);
    }

    return returnValue;
}

std::vector<double>
getRhoDoubleDerivatives(const std::vector<double> & point)
{

    std::vector<double> returnValue(9,0.0);
    const double h = 1.0e-6;
    const int numFDPoints = 13;
    for(unsigned int i = 0; i < 3; ++i)
    {
        for(unsigned int j = 0; j < 3; ++j)
        {
            const int index = i*3 + j;
            std::vector<int> indices(2);
            indices[0] = i;
            indices[1] = j;
            returnValue[index] = getHigherOrderDerivativeFD(point,
                                                            &getRho,
                                                            indices,
                                                            numFDPoints,
                                                            h);
        }
    }

    return returnValue;
}


double
getLaplacianRho(const std::vector<double> & point)
{

    double returnValue = 0.0;
    const double h = 1.0e-6;
    const int numFDPoints = 13;
    for(unsigned int i = 0; i < 3; ++i)
    {
        std::vector<int> indices(2);
        indices[0] = i;
        indices[1] = i;
        returnValue += getHigherOrderDerivativeFD(point,
                                                  &getRho,
                                                  indices,
                                                  numFDPoints,
                                                  h);
    }

    return returnValue;
}

std::vector<double>
getGradModGradRho(const std::vector<double> & point)
{

    std::vector<double> returnValue(3,0.0);
    std::vector<double> gradRho = getGradRho(point);
    const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
    std::vector<double> rhoDoubleDerivatives = getRhoDoubleDerivatives(point);
    for(unsigned int i = 0; i < 3; ++i)
    {

        for(unsigned int j = 0; j < 3; ++j)
        {
            const int index = i*3 + j;
            returnValue[i] += rhoDoubleDerivatives[index]*gradRho[j];
        }

        returnValue[i] /= modGradRho;
    }

    return returnValue;

}

double
getexUniform(const double rho)
{

    return Ax_global*std::pow(rho,4.0/3.0);
}

double
getdexUniformDRho(const double rho)
{

    return 4.0/3.0*Ax_global*std::pow(rho,1.0/3.0);
}

double 
getS(const double rho,
     const double modGradRho)
{
    return As_global*modGradRho/std::pow(rho,4.0/3.0);
}

double
getF(const std::vector<double> & s)
{
    double returnValue = 0.0;
    const double muByKappa = mu_global/kappa_global;
    const double sVal = s[0];
    const double x = muByKappa*sVal*sVal;
    const double y = 1.0 - 1.0/(1.0+x);
    const double z = 1.0 - exp(-x);
    returnValue += ax_global[0] + bx_global[0];
    for(unsigned int i = 1; i <= 5; ++i)
    {
        returnValue += ax_global[i]*std::pow(y,i);
        returnValue += bx_global[i]*std::pow(z,i);
    }

    return returnValue;
}

double
getdFds(const double s)
{

    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(1);
    dependentVariables[0] = s;
    std::vector<int> indices(1);
    indices[0] = 0;
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
                                             &getF,
                                             indices,
                                             numFDPoints,
                                             h);
    return returnValue;
}

    double
getd2Fds2(const double s)
{

    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(1);
    dependentVariables[0] = s;
    std::vector<int> indices(2);
    indices[0] = 0;
    indices[1] = 0;
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
                                             &getF,
                                             indices,
                                             numFDPoints,
                                             h);
    return returnValue;
}

double
getFSelective(const std::vector<double> & s,
             const int index)
{
    double returnValue = 0.0;
    const double muByKappa = mu_global/kappa_global;
    const double sVal = s[0];
    const double x = muByKappa*sVal*sVal;
    const double y = 1.0 - 1.0/(1.0+x);
    const double z = 1.0 - exp(-x);
    if(index < 6)
    {
        if(index == 0)
            returnValue = ax_global[0];
        else
            returnValue = ax_global[index]*std::pow(y,index);
    }
    if(index >=6 && index < 12)
    {
        int indexAdjusted = index - 6;

        if(index == 6)
            returnValue = bx_global[indexAdjusted];
        else
            returnValue = bx_global[indexAdjusted]*std::pow(z,indexAdjusted);

    }

    return returnValue;
}

double
getdFdsSelective(const double s, const int index)
{

    double returnValue = 0.0;
    const double h = 1.0e-4;
    const int numFDPoints = 10;
    std::vector<double> dependentVariables(1);
    dependentVariables[0] = s;
    std::vector<int> indices(1);
    indices[0] = 0;
    returnValue = getHigherOrderDerivativeFDForwardDiffSelective(dependentVariables,
                                             index,
                                             &getFSelective,
                                             indices,
                                             numFDPoints,
                                             h);
    return returnValue;
}

    double
getd2Fds2Selective(const double s, const int index)
{

    double returnValue = 0.0;
    const double h = 1.0e-4;
    const int numFDPoints = 10;
    std::vector<double> dependentVariables(1);
    dependentVariables[0] = s;
    std::vector<int> indices(2);
    indices[0] = 0;
    indices[1] = 0;
    returnValue = getHigherOrderDerivativeFDForwardDiffSelective(dependentVariables,
                                             index,
                                             &getFSelective,
                                             indices,
                                             numFDPoints,
                                             h);
    return returnValue;
}
//double
//getdFds(const double s)
//{
//    double returnValue = 0.0;
//    const double muByKappa = mu_global/kappa_global;
//    const double x = muByKappa*s*s;
//    const double y = 1.0 - 1.0/(1.0+x);
//    const double z = 1.0 - exp(-x);
//    const double dyds = 2.0*muByKappa*s/std::pow(1.0 + muByKappa*s*s, 2.0);
//    const double dzds = 2.0*muByKappa*s*exp(-muByKappa*s*s);
//    for(unsigned int i = 1; i <= 5; ++i)
//    {
//        returnValue += ax_global[i]*i*std::pow(y,i-1)*dyds;
//        returnValue += bx_global[i]*i*std::pow(z,i-1)*dzds;
//    }
//
//    return returnValue;
//}

//double
//getd2Fds2(const double s)
//{
//    double returnValue = 0.0;
//    const double muByKappa = mu_global/kappa_global;
//    const double x = muByKappa*s*s;
//    const double y = 1.0 - 1.0/(1.0+x);
//    const double z = 1.0 - exp(-x);
//    const double dyds = 2.0*muByKappa*s/std::pow(1.0 + muByKappa*s*s, 2.0);
//    const double d2yds2 = 2.0*muByKappa/std::pow(1.0 + muByKappa*s*s, 2.0) - 2.0*std::pow(2.0*muByKappa*s, 2.0)/std::pow(1.0 + muByKappa*s*s, 3.0);
//    const double dzds = 2.0*muByKappa*s*exp(-muByKappa*s*s);
//    const double d2zds2 = exp(-muByKappa*s*s)*(2.0*muByKappa - std::pow(2.0*muByKappa*s, 2.0));
//    for(unsigned int i = 1; i <= 5; ++i)
//    {
//        returnValue += ax_global[i]*i*((i-1)*std::pow(y,i-2)*dyds*dyds + std::pow(y,i-1)*d2yds2);
//        returnValue += bx_global[i]*i*((i-1)*std::pow(z,i-2)*dzds*dzds + std::pow(z,i-1)*d2zds2);
//    }
//
//    return returnValue;
//}

double getDexDSigma(const std::vector<double> & point)
{
    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
        gradRho[i] *= 2.0;
    double sigma = dotProduct(gradRho,gradRho);

     int exceptParamX;
     xc_func_type funcX;
     exceptParamX = xc_func_init(&funcX,XC_GGA_X_SOGGA11,XC_UNPOLARIZED);
     xc_func_set_dens_threshold(&funcX, 1e-10);
     double dexdRho, dexdSigma;
     xc_gga_vxc(&funcX, 1, &rho, &sigma, &dexdRho, &dexdSigma);
     
     return dexdSigma;

}

    std::vector<double>
getGradientDexDSigmaFD(const std::vector<double> & point)
{

    std::vector<double> returnValue(3,0.0);
    const int numFDPoints = 13;
    const double h = 1.0e-4;
    for(int i = 0; i < 3; ++i)
    {
        std::vector<int> indices(1);
        indices[0] = i;
        returnValue[i] = getHigherOrderDerivativeFD(point,
                                                    &getDexDSigma,
                                                    indices,
                                                    numFDPoints,
                                                    h);
    }

    return returnValue;

}


double getFAX(const std::vector<double> & dependentVariables)
{
    const double sVal = dependentVariables[0];
    std::vector<double> params(dependentVariables.begin() + 1, dependentVariables.end());
    double returnValue = 0.0;
    const double muByKappa = mu_global/kappa_global;
    const double x = muByKappa*sVal*sVal;
    const double y = 1.0 - 1.0/(1.0+x);
    const int numParams = params.size();
    returnValue += params[0];
    for(unsigned int i = 1; i < numParams; ++i)
    {
        returnValue += params[i]*std::pow(y,i);
    }
    return returnValue;

}

double getFBX(const std::vector<double> & dependentVariables)
{
    const double sVal = dependentVariables[0];
    std::vector<double> params(dependentVariables.begin() + 1, dependentVariables.end());
    double returnValue = 0.0;
    const double muByKappa = mu_global/kappa_global;
    const double x = muByKappa*sVal*sVal;
    const double z = 1.0 - exp(-x);
    const int numParams = params.size();
    returnValue += params[0];
    for(unsigned int i = 1; i < numParams; ++i)
    {
        returnValue += params[i]*std::pow(z,i);
    }
    
    return returnValue;

}


    double 
getFAXHigherOrderDerivative(const std::vector<double> & dependentVariables,
        const std::vector<int> & indices)
{
    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<int> indicesCopy(indices);
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
            &getFAX,
            indicesCopy,
            numFDPoints,
            h);
    
    return returnValue;


}

    double 
getFBXHigherOrderDerivative(const std::vector<double> & dependentVariables,
        const std::vector<int> & indices)
{
    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<int> indicesCopy(indices);
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
            &getFBX,
            indicesCopy,
            numFDPoints,
            h);
    
    return returnValue;

}
    double
getFXHigherOrderDerivative(const double sVal,
        const std::vector<double> & axParams,
        const std::vector<double> & bxParams,
        const std::vector<int> & indices,
        const int ABId)
{

    double returnValue = 0.0;
    if(ABId == 0) // Include both A and B terms
    {
        std::vector<double> FAXVariables(0);
        FAXVariables.push_back(sVal);
        FAXVariables.insert(FAXVariables.end(), axParams.begin(), axParams.end());
        returnValue += getFAXHigherOrderDerivative(FAXVariables,
                indices);

        std::vector<double> FBXVariables(0);
        FBXVariables.push_back(sVal);
        FBXVariables.insert(FBXVariables.end(), bxParams.begin(), bxParams.end());
        returnValue += getFBXHigherOrderDerivative(FBXVariables,
                indices);
    }

    else if(ABId == 1) // Include only A terms
    {

        std::vector<double> FAXVariables(0);
        FAXVariables.push_back(sVal);
        FAXVariables.insert(FAXVariables.end(), axParams.begin(), axParams.end());

        returnValue += getFAXHigherOrderDerivative(FAXVariables,
                indices);
    }

    else if(ABId == 2) // Include only B terms
    {
        std::vector<double> FBXVariables(0);
        FBXVariables.push_back(sVal);
        FBXVariables.insert(FBXVariables.end(), bxParams.begin(), bxParams.end());
        returnValue += getFBXHigherOrderDerivative(FBXVariables,
                indices);
    }
    else
    {

        const std::string message("Invalid ABId provided.");
        std::cout << message << std::endl;
    }

    return returnValue;
}

    double
getVXCFromLibXC(const std::vector<double> & point)
{

    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    double laplacianRho = getLaplacianRho(point);

    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
        gradRho[i] *= 2.0;

    double sigma = dotProduct(gradRho,gradRho);
    laplacianRho *= 2.0;

    std::vector<double> gradDexDSigma = 
        getGradientDexDSigmaFD(point);

    int exceptParamX;
    xc_func_type funcX;
    exceptParamX = xc_func_init(&funcX,XC_GGA_X_SOGGA11,XC_UNPOLARIZED);
    xc_func_set_dens_threshold(&funcX, 1e-10);
    double dexdRho, dexdSigma;
    xc_gga_vxc(&funcX, 1, &rho, &sigma, &dexdRho, &dexdSigma);
    const double term1 = dexdRho; 
    const double term2 = -2.0*dotProduct(gradDexDSigma, gradRho);
    const double term3 = -2.0*laplacianRho*dexdSigma;
    const double returnValue = term1 + term2 + term3; 

    return returnValue;
}


    double
getVXTest(const std::vector<double> & point)
{

    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    double laplacianRho = getLaplacianRho(point);
    std::vector<double> gradModGradRho = getGradModGradRho(point);
    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
    {
        gradRho[i] *= 2.0;
        gradModGradRho[i] *= 2.0;
    }
    laplacianRho *= 2.0;

    std::vector<double> d_paramsAX(ax_global, ax_global+6);
    std::vector<double> d_paramsBX(bx_global, bx_global+6);
    const double rho4by3 = std::pow(rho,4.0/3.0);
    const double rho8by3 = std::pow(rho,8.0/3.0);
    const double rho11by3 = std::pow(rho,11.0/3.0);
    const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
    const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
    const double s = getS(rho, modGradRho);
    const double exUniform = getexUniform(rho);
    const double dexUniformdRho = getdexUniformDRho(rho);

    int ABId = 0;
    std::vector<int> noIndices(0);
    const double F = getFXHigherOrderDerivative(s, d_paramsAX, d_paramsBX, noIndices, ABId);
    std::vector<int> indicesDs(1);
    indicesDs[0] = 0;
    const double dFds = getFXHigherOrderDerivative(s, d_paramsAX, d_paramsBX, indicesDs, ABId);
    std::vector<int> indicesDs2(2);
    indicesDs2[0] = 0;
    indicesDs2[1] = 0;
    const double d2Fds2 =  getFXHigherOrderDerivative(s, d_paramsAX, d_paramsBX, indicesDs2, ABId);

    const double t1 = dexUniformdRho*F;
    const double t21 = -dexUniformdRho*modGradRho/rho4by3;
    const double t22 = exUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t23 = -exUniform*laplacianRho/(rho4by3*modGradRho);
    const double t2 = As_global*dFds*(t21 + t22 + t23);
    const double t31 = (4.0/3.0)*exUniform*modGradRho*modGradRho/rho11by3;
    const double t32 = -exUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t3 = As_global*As_global*d2Fds2*(t31 + t32);
    const double returnValue = t1 + t2 + t3;

    return returnValue;
}

    double
getVXC(const std::vector<double> & point)
{

    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    double laplacianRho = getLaplacianRho(point);
    std::vector<double> gradModGradRho = getGradModGradRho(point);
    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
    {
        gradRho[i] *= 2.0;
        gradModGradRho[i] *= 2.0;
    }
    laplacianRho *= 2.0;

    const double rho4by3 = std::pow(rho,4.0/3.0);
    const double rho8by3 = std::pow(rho,8.0/3.0);
    const double rho11by3 = std::pow(rho,11.0/3.0);
    const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
    const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
    const double s = getS(rho, modGradRho);
    const double exUniform = getexUniform(rho);
    const double dexUniformdRho = getdexUniformDRho(rho);
    std::vector<double> sDummyVec(1);
    sDummyVec[0] = s;
    const double F = getF(sDummyVec);
    const double dFds = getdFds(s);
    const double d2Fds2 = getd2Fds2(s);

    const double t1 = dexUniformdRho*F;
    const double t21 = -dexUniformdRho*modGradRho/rho4by3;
    const double t22 = exUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t23 = -exUniform*laplacianRho/(rho4by3*modGradRho);
    const double t2 = As_global*dFds*(t21 + t22 + t23);
    const double t31 = (4.0/3.0)*exUniform*modGradRho*modGradRho/rho11by3;
    const double t32 = -exUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t3 = As_global*As_global*d2Fds2*(t31 + t32);
    const double returnValue = t1 + t2 + t3;
    //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")" << " vxc: " << returnValue << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << std::endl; 
    return returnValue;

} 
    double
getdVXCdParam(const std::vector<double> & point, const int paramId)
{

    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    double laplacianRho = getLaplacianRho(point);
    std::vector<double> gradModGradRho = getGradModGradRho(point);
    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
    {
        gradRho[i] *= 2.0;
        gradModGradRho[i] *= 2.0;
    }
    laplacianRho *= 2.0;

    const double rho4by3 = std::pow(rho,4.0/3.0);
    const double rho8by3 = std::pow(rho,8.0/3.0);
    const double rho11by3 = std::pow(rho,11.0/3.0);
    const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
    const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
    const double s = getS(rho, modGradRho);
    const double exUniform = getexUniform(rho);
    const double dexUniformdRho = getdexUniformDRho(rho);
    std::vector<double> sDummyVec(1);
    sDummyVec[0] = s;
    const double F = getFSelective(sDummyVec, paramId);
    const double dFds = getdFdsSelective(s, paramId);
    const double d2Fds2 = getd2Fds2Selective(s, paramId);

    const double t1 = dexUniformdRho*F;
    const double t21 = -dexUniformdRho*modGradRho/rho4by3;
    const double t22 = exUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t23 = -exUniform*laplacianRho/(rho4by3*modGradRho);
    const double t2 = As_global*dFds*(t21 + t22 + t23);
    const double t31 = (4.0/3.0)*exUniform*modGradRho*modGradRho/rho11by3;
    const double t32 = -exUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t3 = As_global*As_global*d2Fds2*(t31 + t32);
    double factor;
    if(paramId < 6)
        factor = ax_global[paramId];
    else
        factor = bx_global[paramId-6];
    const double returnValue = (t1 + t2 + t3)/factor;
    //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")" << " vxc: " << returnValue << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << std::endl; 
    return returnValue;

}


    bool
isAParamX(const int xParamId)
{
    if(xParamId < 6)
        return true;
    else if(xParamId >= 6 && xParamId < 12)
        return false;
    else
    {
        const std::string message("Invalid xParamId provided.");
        std::cout << message << std::endl;
    }
}

    int 
getXABAdjustedParamIdFromXParamId(const int xParamId)
{
    if(xParamId < 6)
        return xParamId;
    else if(xParamId >= 6 && xParamId < 12)
        return (xParamId - 6);
    else
    {
        const std::string message("Invalid xParamId provided.");
        std::cout << message << std::endl;
    }
}



    double
getdVXdParamTest(const std::vector<double> & point,
        const int xParamId)
{
    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    double laplacianRho = getLaplacianRho(point);
    std::vector<double> gradModGradRho = getGradModGradRho(point);
    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
    {
        gradRho[i] *= 2.0;
        gradModGradRho[i] *= 2.0;
    }
    laplacianRho *= 2.0;

    int d_paramIndexOffsetX = 1;
    std::vector<double> d_paramsAX(ax_global, ax_global+6);
    std::vector<double> d_paramsBX(bx_global, bx_global+6);

    const double rho4by3 = std::pow(rho,4.0/3.0);
    const double rho8by3 = std::pow(rho,8.0/3.0);
    const double rho11by3 = std::pow(rho,11.0/3.0);
    const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
    const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
    const double s = getS(rho, modGradRho);
    const double exUniform = getexUniform(rho);
    const double dexUniformdRho = getdexUniformDRho(rho);

    int ABId;
    bool isAPresent = isAParamX(xParamId);
    if(isAPresent)
        ABId = 1;
    else
        ABId = 2;
    int ABAdjustedParamId = getXABAdjustedParamIdFromXParamId(xParamId);

    std::vector<int> indicesDParam(1);
    indicesDParam[0] = ABAdjustedParamId + d_paramIndexOffsetX;
    const double dFdParam = getFXHigherOrderDerivative(s, d_paramsAX, d_paramsBX, indicesDParam, ABId);
    std::vector<int> indicesDsDParam(2);
    indicesDsDParam[0] = 0;
    indicesDsDParam[1] = ABAdjustedParamId + d_paramIndexOffsetX;
    const double d2FdsdParam = getFXHigherOrderDerivative(s, d_paramsAX, d_paramsBX, indicesDsDParam, ABId);
    std::vector<int> indicesDs2DParam(3);
    indicesDs2DParam[0] = 0;
    indicesDs2DParam[1] = 0;
    indicesDs2DParam[2] = ABAdjustedParamId + d_paramIndexOffsetX;
    const double d3Fds2dParam =  getFXHigherOrderDerivative(s, d_paramsAX, d_paramsBX, indicesDs2DParam, ABId);

    const double t1 = dexUniformdRho*dFdParam;
    const double t21 = -dexUniformdRho*modGradRho/rho4by3;
    const double t22 = exUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t23 = -exUniform*laplacianRho/(rho4by3*modGradRho);
    const double t2 = As_global*d2FdsdParam*(t21 + t22 + t23);
    const double t31 = (4.0/3.0)*exUniform*modGradRho*modGradRho/rho11by3;
    const double t32 = -exUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t3 = As_global*As_global*d3Fds2dParam*(t31 + t32);
    const double returnValue = t1 + t2 + t3;

    return returnValue;
}


int main()
{
    int numPoints = 10;
    srand(0);
    for(int LIndex = -3; LIndex <= 2; ++LIndex)
    {
        const double L = std::pow(10.0,LIndex);
        std::cout << "\nPrinting for L: " << L << std::endl;
        std::cout << "----------------------------------------------------------------------------------------------------------" << std::endl;
        for(unsigned int i = 0; i < numPoints; ++i)
        {
            std::vector<double> point(3);
            for(unsigned int j = 0; j < 3; ++j)
                point[j] = L*(((rand()+0.0)/RAND_MAX)-0.5);

            const double rho = getRho(point);
            const std::vector<double> gradRho = getGradRho(point);
            const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
            const double s = getS(rho, modGradRho);
            const double laplacianRho = getLaplacianRho(point);
            const std::vector<double> gradModGradRho = getGradModGradRho(point);
            const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
            std::vector<double> sDummyVec(1);
            sDummyVec[0] = s;
            const double F = getF(sDummyVec);
            const double dFds = getdFds(s);
            const double d2Fds2 = getd2Fds2(s);
            double vxc = getVXC(point);
            double vxcTest = getVXTest(point);
            double vxcLIBXC = getVXCFromLibXC(point);
            int xParamId = rand()%12;
            double factor;
            if(xParamId < 6)
                factor = ax_global[xParamId];
            else
                factor = bx_global[xParamId-6];
            double dvxcdParam = getdVXCdParam(point, xParamId);
            double dvxcdParamTest = getdVXdParamTest(point, xParamId);
            //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")\trho: " << rho << "\tmodGradRho: " << modGradRho << "\t laplacianRho: " << laplacianRho << "\tgradModGradRhoDotGradRho: " << gradModGradRhoDotGradRho <<"\ts: " << s << "\t F: " << F << "\tF': " << dFds << "\tF'': " << d2Fds2 << "\tvxc: " << vxc << "\t" << vxcLIBXC << "\tRel. Err.: " << std::abs((vxc-vxcLIBXC)/vxcLIBXC) << std::endl;
            std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")\trho: " << rho  << "\tvxc: " << vxc << "\t" << "\tvxcTest: " << vxcTest << "\tvxc Error: " << std::abs(vxc-vxcLIBXC) << "\tvxcTest Error: " << std::abs((vxc-vxcTest)) << "\tParamId: " << xParamId << "\tdvxcdParam: " << dvxcdParam <<  "\tdvxcdParamTest: " << dvxcdParamTest << "\tdvxcdParamError: " << std::abs(dvxcdParam-dvxcdParamTest) << std::endl;
            //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")\tvxc: " << vxc << "\t" << "\tvxcTest: " << vxcTest << "\tvxcTestError: " << std::abs(vxc-vxcTest) << "\t vxcLibXC: " << vxcLIBXC << "\tRel. Err.: " << std::abs((vxc-vxcLIBXC)/vxcLIBXC) << std::endl;
        }
    }
}
