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

double ac_global[] = {0.50000, -4.62334, 8.00410, -130.226, 38.2685, 69.5599};
double bc_global[] = {0.50000, 3.62334, 9.36393, 34.5114, -18.5684, -0.16519};
double As_global = 1.0/(2.0*std::pow(3.0*M_PI*M_PI,1.0/3.0));
double beta_global =  0.066725; 
double alpha_global = std::pow(3,1.0/3.0)*std::pow(M_PI,5.0/3.0)/4.0;

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
getepsiloncUniform(const double rho)
{  
    int exceptParamC;
    xc_func_type funcC;
    exceptParamC = xc_func_init(&funcC,XC_LDA_C_PW,XC_UNPOLARIZED);
    xc_func_set_dens_threshold(&funcC, 1e-10);
    double epsilonc;
    double rhoCopy = rho;
    xc_lda_exc(&funcC, 1, &rhoCopy, &epsilonc);
    xc_func_end(&funcC);

    return epsilonc;
}

    double
getecUniform(const double rho)
{
    return rho*getepsiloncUniform(rho);
}

    double
getdecUniformDRho(const double rho)
{

    int exceptParamC;
    xc_func_type funcC;
    exceptParamC = xc_func_init(&funcC,XC_LDA_C_PW,XC_UNPOLARIZED);
    xc_func_set_dens_threshold(&funcC, 1e-10);
    double decdRho;
    double rhoCopy = rho;
    xc_lda_vxc(&funcC, 1, &rhoCopy, &decdRho);
    xc_func_end(&funcC);

    return decdRho;
}

    double 
getS(const double rho,
        const double modGradRho)
{
    return As_global*modGradRho/std::pow(rho,4.0/3.0);
}

    double
getF(const std::vector<double> & sAndRho)
{
    double returnValue = 0.0;
    const double sVal = sAndRho[0];
    const double rhoVal = sAndRho[1];
    const double rhoVal1By3 = std::pow(rhoVal,1.0/3.0);
    const double epsilonc = getepsiloncUniform(rhoVal);

    const double x = alpha_global*beta_global*sVal*sVal*rhoVal1By3/epsilonc;
    const double y = 1.0 - 1.0/(1.0-x);
    const double z = 1.0 - exp(x);
    returnValue += ac_global[0] + bc_global[0];
    for(unsigned int i = 1; i <= 5; ++i)
    {
        returnValue += ac_global[i]*std::pow(y,i);
        returnValue += bc_global[i]*std::pow(z,i);
    }

    return returnValue;
}

    double
getdFds(const std::vector<double> & sAndRho)
{

    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
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
getdFdRho(const std::vector<double> & sAndRho)
{

    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
    std::vector<int> indices(1);
    indices[0] = 1;
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
            &getF,
            indices,
            numFDPoints,
            h);
    return returnValue;
}

    double
getd2Fds2(const std::vector<double> & sAndRho)
{

    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
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
getd2FdsdRho(const std::vector<double> & sAndRho)
{

    double returnValue = 0.0;
    const double h = 1.0e-2;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
    std::vector<int> indices(2);
    indices[0] = 0;
    indices[1] = 1;
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
            &getF,
            indices,
            numFDPoints,
            h);
    return returnValue;
}



double
getFSelective(const std::vector<double> & sAndRho,
             const int index)
{
    double returnValue = 0.0;
    const double sVal = sAndRho[0];
    const double rhoVal = sAndRho[1];
    const double rhoVal1By3 = std::pow(rhoVal,1.0/3.0);
    const double epsilonc = getepsiloncUniform(rhoVal);
    const double x = alpha_global*beta_global*sVal*sVal*rhoVal1By3/epsilonc;
    const double y = 1.0 - 1.0/(1.0-x);
    const double z = 1.0 - exp(x);
    if(index < 6)
    {
        if(index == 0)
            returnValue = ac_global[0];
        else
            returnValue = ac_global[index]*std::pow(y,index);
    }
    if(index >=6 && index < 12)
    {
        int indexAdjusted = index - 6;

        if(index == 6)
            returnValue = bc_global[indexAdjusted];
        else
            returnValue = bc_global[indexAdjusted]*std::pow(z,indexAdjusted);

    }

    return returnValue;
}

double
getdFdsSelective(const std::vector<double> & sAndRho, const int index)
{

    double returnValue = 0.0;
    const double h = 1.0e-3;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
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
getdFdRhoSelective(const std::vector<double> & sAndRho, const int index)
{

    double returnValue = 0.0;
    const double h = 1.0e-3;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
    std::vector<int> indices(1);
    indices[0] = 1;
    returnValue = getHigherOrderDerivativeFDForwardDiffSelective(dependentVariables,
                                             index,
                                             &getFSelective,
                                             indices,
                                             numFDPoints,
                                             h);
    return returnValue;
}

    double
getd2Fds2Selective(const std::vector<double> & sAndRho, const int index)
{

    double returnValue = 0.0;
    const double h = 1.0e-3;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
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

    double
getd2FdsdRhoSelective(const std::vector<double> & sAndRho, const int index)
{

    double returnValue = 0.0;
    const double h = 1.0e-3;
    const int numFDPoints = 12;
    std::vector<double> dependentVariables(sAndRho);
    std::vector<int> indices(2);
    indices[0] = 0;
    indices[1] = 1;
    returnValue = getHigherOrderDerivativeFDForwardDiffSelective(dependentVariables,
                                             index,
                                             &getFSelective,
                                             indices,
                                             numFDPoints,
                                             h);
    return returnValue;
}


double getFAC(const std::vector<double> & dependentVariables)
{
    const double sVal = dependentVariables[0];
    const double rhoVal = dependentVariables[1];
    std::vector<double> params(dependentVariables.begin() + 2, dependentVariables.end());
    double returnValue = 0.0;
    const double rhoVal1By3 = std::pow(rhoVal,1.0/3.0);
    const double epsilonc = getepsiloncUniform(rhoVal);
    const double x = alpha_global*beta_global*sVal*sVal*rhoVal1By3/epsilonc;
    const double y = 1.0 - 1.0/(1.0-x);
    const int numParams = params.size();
    returnValue += params[0];
    for(unsigned int i = 1; i < numParams; ++i)
    {
        returnValue += params[i]*std::pow(y,i);
    }
    return returnValue;

}

double getFBC(const std::vector<double> & dependentVariables)
{
    const double sVal = dependentVariables[0];
    const double rhoVal = dependentVariables[1];
    std::vector<double> params(dependentVariables.begin() + 2, dependentVariables.end());
    double returnValue = 0.0;
    const double rhoVal1By3 = std::pow(rhoVal,1.0/3.0);
    const double epsilonc = getepsiloncUniform(rhoVal);
    const double x = alpha_global*beta_global*sVal*sVal*rhoVal1By3/epsilonc;
    const double z = 1.0 - exp(x);
    const int numParams = params.size();
    returnValue += params[0];
    for(unsigned int i = 1; i < numParams; ++i)
    {
        returnValue += params[i]*std::pow(z,i);
    }

    return returnValue;

}


    double 
getFACHigherOrderDerivative(const std::vector<double> & dependentVariables,
        const std::vector<int> & indices)
{
    double returnValue = 0.0;
    const double h = 1.0e-3;
    const int numFDPoints = 10;
    std::vector<int> indicesCopy(indices);
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
            &getFAC,
            indicesCopy,
            numFDPoints,
            h);

    return returnValue;


}

    double 
getFBCHigherOrderDerivative(const std::vector<double> & dependentVariables,
        const std::vector<int> & indices)
{
    double returnValue = 0.0;
    const double h = 1.0e-3;
    const int numFDPoints = 10;
    std::vector<int> indicesCopy(indices);
    returnValue = getHigherOrderDerivativeFDForwardDiff(dependentVariables,
            &getFBC,
            indicesCopy,
            numFDPoints,
            h);

    return returnValue;

}
    double
getFCHigherOrderDerivative(const double sVal,
        const double rhoVal,
        const std::vector<double> & acParams,
        const std::vector<double> & bcParams,
        const std::vector<int> & indices,
        const int ABId)
{

    double returnValue = 0.0;
    if(ABId == 0) // Include both A and B terms
    {
        std::vector<double> FACVariables(0);
        FACVariables.push_back(sVal);
        FACVariables.push_back(rhoVal);
        FACVariables.insert(FACVariables.end(), acParams.begin(), acParams.end());
        returnValue += getFACHigherOrderDerivative(FACVariables,
                indices);

        std::vector<double> FBCVariables(0);
        FBCVariables.push_back(sVal);
        FBCVariables.push_back(rhoVal);
        FBCVariables.insert(FBCVariables.end(), bcParams.begin(), bcParams.end());
        returnValue += getFBCHigherOrderDerivative(FBCVariables,
                indices);
    }

    else if(ABId == 1) // Include only A terms
    {

        std::vector<double> FACVariables(0);
        FACVariables.push_back(sVal);
        FACVariables.push_back(rhoVal);
        FACVariables.insert(FACVariables.end(), acParams.begin(), acParams.end());
        returnValue += getFACHigherOrderDerivative(FACVariables,
                indices);
    }

    else if(ABId == 2) // Include only B terms
    {
        std::vector<double> FBCVariables(0);
        FBCVariables.push_back(sVal);
        FBCVariables.push_back(rhoVal);
        FBCVariables.insert(FBCVariables.end(), bcParams.begin(), bcParams.end());
        returnValue += getFBCHigherOrderDerivative(FBCVariables,
                indices);
    }
    else
    {

        const std::string message("Invalid ABId provided.");
        std::cout << message << std::endl;
    }

    return returnValue;
}


double getDecDSigma(const std::vector<double> & point)
{
    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
        gradRho[i] *= 2.0;
    double sigma = dotProduct(gradRho,gradRho);

    int exceptParamC;
    xc_func_type funcC;
    exceptParamC = xc_func_init(&funcC,XC_GGA_C_SOGGA11,XC_UNPOLARIZED);
    xc_func_set_dens_threshold(&funcC, 1e-10);
    double decdRho, decdSigma;
    xc_gga_vxc(&funcC, 1, &rho, &sigma, &decdRho, &decdSigma);
    xc_func_end(&funcC);

    return decdSigma;

}

    std::vector<double>
getGradientDecDSigmaFD(const std::vector<double> & point)
{

    std::vector<double> returnValue(3,0.0);
    const int numFDPoints = 13;
    const double h = 1.0e-4;
    for(int i = 0; i < 3; ++i)
    {
        std::vector<int> indices(1);
        indices[0] = i;
        returnValue[i] = getHigherOrderDerivativeFD(point,
                &getDecDSigma,
                indices,
                numFDPoints,
                h);
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
    laplacianRho *= 2.0;

    double sigma = dotProduct(gradRho,gradRho);
    std::vector<double> gradDecDSigma = 
        getGradientDecDSigmaFD(point);

    int exceptParamC;
    xc_func_type funcC;
    exceptParamC = xc_func_init(&funcC, XC_GGA_C_SOGGA11, XC_UNPOLARIZED);
    xc_func_set_dens_threshold(&funcC, 1e-10);
    double decdRho, decdSigma;
    xc_gga_vxc(&funcC, 1, &rho, &sigma, &decdRho, &decdSigma);
    const double term1 = decdRho; 
    const double term2 = -2.0*dotProduct(gradDecDSigma, gradRho);
    const double term3 = -2.0*laplacianRho*decdSigma;
    const double returnValue = term1 + term2 + term3; 
    xc_func_end(&funcC);

    return returnValue;
}

    double
getVXCTest(const std::vector<double> & point)
{

    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    std::vector<double> gradModGradRho = getGradModGradRho(point);
    double laplacianRho = getLaplacianRho(point);

    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
    {
        gradRho[i] *= 2.0;
        gradModGradRho[i] *= 2.0;
    }
    laplacianRho *= 2.0;

    std::vector<double> d_paramsAC(ac_global, ac_global+6);
    std::vector<double> d_paramsBC(bc_global, bc_global+6);
    const double rho4by3 = std::pow(rho,4.0/3.0);
    const double rho8by3 = std::pow(rho,8.0/3.0);
    const double rho11by3 = std::pow(rho,11.0/3.0);
    const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
    const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
    const double s = getS(rho, modGradRho);

    const double ecUniform = getecUniform(rho);
    const double decUniformdRho = getdecUniformDRho(rho);
    std::vector<double> sAndRho(2);
    sAndRho[0] = s;
    sAndRho[1] = rho;


    int ABId = 0;
    std::vector<int> noIndices(0);
    const double F = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, noIndices, ABId);
    std::vector<int> indicesDs(1);
    indicesDs[0] = 0;
    const double dFds = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDs, ABId);
    std::vector<int> indicesDRho(1);
    indicesDRho[0] = 1;
    const double dFdRho = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDRho, ABId);
    std::vector<int> indicesDs2(2);
    indicesDs2[0] = 0;
    indicesDs2[1] = 0;
    const double d2Fds2 = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDs2, ABId);
    std::vector<int> indicesDsDRho(2);
    indicesDsDRho[0] = 0;
    indicesDsDRho[1] = 1;
    const double d2FdsdRho = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDsDRho, ABId);

    const double t1 = decUniformdRho*F;
    const double t2 = ecUniform*dFdRho;
    const double t31 = -decUniformdRho*modGradRho/rho4by3;
    const double t32 = ecUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t33 = -ecUniform*laplacianRho/(rho4by3*modGradRho);
    const double t3 = As_global*dFds*(t31 + t32 + t33);
    const double t41 = (4.0/3.0)*ecUniform*modGradRho*modGradRho/rho11by3;
    const double t42 = -ecUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t4 = As_global*As_global*d2Fds2*(t41 + t42);
    const double t5 = -As_global*d2FdsdRho*ecUniform*modGradRho/rho4by3;
    const double returnValue = t1 + t2 + t3 + t4 + t5;
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
    const double ecUniform = getecUniform(rho);
    const double decUniformdRho = getdecUniformDRho(rho);
    std::vector<double> sAndRho(2);
    sAndRho[0] = s;
    sAndRho[1] = rho;
    const double F = getFSelective(sAndRho, paramId);
    const double dFds = getdFdsSelective(sAndRho, paramId);
    const double dFdRho = getdFdRhoSelective(sAndRho, paramId);
    const double d2Fds2 = getd2Fds2Selective(sAndRho, paramId);
    const double d2FdsdRho = getd2FdsdRhoSelective(sAndRho, paramId);

    const double t1 = decUniformdRho*F;
    const double t2 = ecUniform*dFdRho;
    const double t31 = -decUniformdRho*modGradRho/rho4by3;
    const double t32 = ecUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t33 = -ecUniform*laplacianRho/(rho4by3*modGradRho);
    const double t3 = As_global*dFds*(t31 + t32 + t33);
    const double t41 = (4.0/3.0)*ecUniform*modGradRho*modGradRho/rho11by3;
    const double t42 = -ecUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t4 = As_global*As_global*d2Fds2*(t41 + t42);
    const double t5 = -As_global*d2FdsdRho*ecUniform*modGradRho/rho4by3;
    double factor;
    if(paramId < 6)
        factor = ac_global[paramId];
    else
        factor = bc_global[paramId-6];
    const double returnValue = (t1 + t2 + t3)/factor;

    return returnValue;

}


    bool
isAParamC(const int cParamId)
{
    if(cParamId < 6)
        return true;
    else if(cParamId >= 6 && cParamId < 12)
        return false;
    else
    {
        const std::string message("Invalid cParamId provided.");
        std::cout << message << std::endl;
    }
}

    int 
getCABAdjustedParamIdFromCParamId(const int cParamId)
{
    if(cParamId < 6)
        return cParamId;
    else if(cParamId >= 6 && cParamId < 12)
        return (cParamId - 6);
    else
    {
        const std::string message("Invalid cParamId provided.");
        std::cout << message << std::endl;
    }
}

    double
getdVXCdParamTest(const std::vector<double> & point, const int cParamId)
{

    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    std::vector<double> gradModGradRho = getGradModGradRho(point);
    double laplacianRho = getLaplacianRho(point);

    rho *= 2.0;
    for(unsigned int i = 0; i < 3; ++i)
    {
        gradRho[i] *= 2.0;
        gradModGradRho[i] *= 2.0;
    }
    laplacianRho *= 2.0;

    int d_paramIndexOffsetC = 2;
    std::vector<double> d_paramsAC(ac_global, ac_global+6);
    std::vector<double> d_paramsBC(bc_global, bc_global+6);
    const double rho4by3 = std::pow(rho,4.0/3.0);
    const double rho8by3 = std::pow(rho,8.0/3.0);
    const double rho11by3 = std::pow(rho,11.0/3.0);
    const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
    const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
    const double s = getS(rho, modGradRho);

    const double ecUniform = getecUniform(rho);
    const double decUniformdRho = getdecUniformDRho(rho);
    std::vector<double> sAndRho(2);
    sAndRho[0] = s;
    sAndRho[1] = rho;


    int ABId;
    bool isAPresent = isAParamC(cParamId);
    if(isAPresent)
        ABId = 1;
    else
        ABId = 2;
    int ABAdjustedParamId = getCABAdjustedParamIdFromCParamId(cParamId);
    
    
    std::vector<int> indicesDParam(1);
    indicesDParam[0] = ABAdjustedParamId + d_paramIndexOffsetC;
    const double dFdParam = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDParam, ABId);

    std::vector<int> indicesDsDParam(2);
    indicesDsDParam[0] = 0;
    indicesDsDParam[1] = ABAdjustedParamId + d_paramIndexOffsetC;
    const double d2FdsdParam = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDsDParam, ABId);
    
    std::vector<int> indicesDRhoDParam(2);
    indicesDRhoDParam[0] = 1;
    indicesDRhoDParam[1] =  ABAdjustedParamId + d_paramIndexOffsetC;
    const double d2FdRhodParam = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDRhoDParam, ABId);

    std::vector<int> indicesDs2DParam(3);
    indicesDs2DParam[0] = 0;
    indicesDs2DParam[1] = 0;
    indicesDs2DParam[2] = ABAdjustedParamId + d_paramIndexOffsetC;
    const double d3Fds2dParam = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDs2DParam, ABId);

    std::vector<int> indicesDsDRhoDParam(3);
    indicesDsDRhoDParam[0] = 0;
    indicesDsDRhoDParam[1] = 1;
    indicesDsDRhoDParam[2] = ABAdjustedParamId + d_paramIndexOffsetC;
    const double d3FdsdRhodParam = getFCHigherOrderDerivative(s,rho, d_paramsAC, d_paramsBC, indicesDsDRhoDParam, ABId);

    const double t1 = decUniformdRho*dFdParam;
    const double t2 = ecUniform*d2FdRhodParam;
    const double t31 = -decUniformdRho*modGradRho/rho4by3;
    const double t32 = ecUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t33 = -ecUniform*laplacianRho/(rho4by3*modGradRho);
    const double t3 = As_global*d2FdsdParam*(t31 + t32 + t33);
    const double t41 = (4.0/3.0)*ecUniform*modGradRho*modGradRho/rho11by3;
    const double t42 = -ecUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t4 = As_global*As_global*d3Fds2dParam*(t41 + t42);
    const double t5 = -As_global*d3FdsdRhodParam*ecUniform*modGradRho/rho4by3;
    const double returnValue = t1 + t2 + t3 + t4 + t5;
    //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")" << " vxc: " << returnValue << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << std::endl; 
    return returnValue;

}
    double
getVXC(const std::vector<double> & point)
{

    double rho = getRho(point);
    std::vector<double> gradRho = getGradRho(point);
    std::vector<double> gradModGradRho = getGradModGradRho(point);
    double laplacianRho = getLaplacianRho(point);

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

    const double ecUniform = getecUniform(rho);
    const double decUniformdRho = getdecUniformDRho(rho);
    std::vector<double> sAndRho(2);
    sAndRho[0] = s;
    sAndRho[1] = rho;
    const double F = getF(sAndRho);
    const double dFds = getdFds(sAndRho);
    const double dFdRho = getdFdRho(sAndRho);
    const double d2Fds2 = getd2Fds2(sAndRho);
    const double d2FdsdRho = getd2FdsdRho(sAndRho);

    const double t1 = decUniformdRho*F;
    const double t2 = ecUniform*dFdRho;
    const double t31 = -decUniformdRho*modGradRho/rho4by3;
    const double t32 = ecUniform*gradModGradRhoDotGradRho/(rho4by3*modGradRho*modGradRho);
    const double t33 = -ecUniform*laplacianRho/(rho4by3*modGradRho);
    const double t3 = As_global*dFds*(t31 + t32 + t33);
    const double t41 = (4.0/3.0)*ecUniform*modGradRho*modGradRho/rho11by3;
    const double t42 = -ecUniform*gradModGradRhoDotGradRho/(rho8by3*modGradRho);
    const double t4 = As_global*As_global*d2Fds2*(t41 + t42);
    const double t5 = -As_global*d2FdsdRho*ecUniform*modGradRho/rho4by3;
    const double returnValue = t1 + t2 + t3 + t4 + t5;
    //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")" << " vxc: " << returnValue << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << std::endl; 
    return returnValue;

}

int main()
{
    int numPoints = 10;
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

            double rho = getRho(point);
            std::vector<double> gradRho = getGradRho(point);
            std::vector<double> gradModGradRho = getGradModGradRho(point);
            double laplacianRho = getLaplacianRho(point);
            rho *= 2.0;
            for(unsigned int i = 0; i < 3; ++i)
            {
                gradRho[i] *= 2.0;
                gradModGradRho[i] *= 2.0;
            }
            laplacianRho *= 2.0;

            const double modGradRho = std::sqrt(dotProduct(gradRho,gradRho));
            const double s = getS(rho, modGradRho);
            const double gradModGradRhoDotGradRho = dotProduct(gradModGradRho, gradRho);
            std::vector<double> sAndRho(2);
            sAndRho[0] = s;
            sAndRho[1] = rho;
            const double F = getF(sAndRho);
            const double dFds = getdFds(sAndRho);
            const double d2Fds2 = getd2Fds2(sAndRho);
            double vxc = getVXC(point);
            double vxcTest = getVXCTest(point);
            double vxcLIBXC = getVXCFromLibXC(point);
            int cParamId = rand()%12;
            double factor;
            if(cParamId < 6)
                factor = ac_global[cParamId];
            else
                factor = bc_global[cParamId-6];
            double dvxcdParam = getdVXCdParam(point, cParamId);
            double dvxcdParamTest = getdVXCdParamTest(point, cParamId);

            std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")\trho: " << rho  << "\tvxc: " << vxc << "\t" << "\tvxcTest: " << vxcTest << "\tvxc Error: " << std::abs(vxc-vxcLIBXC) << "\tvxcTest Error: " << std::abs((vxcTest-vxcLIBXC)) << "\tParamId: " << cParamId << "\tdvxcdParam: " << dvxcdParam <<  "\tdvxcdParamTest: " << dvxcdParamTest << "\tdvxcdParamError: " << std::abs(dvxcdParam-dvxcdParamTest) << std::endl;
            //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")\trho: " << rho << "\tmodGradRho: " << modGradRho << "\t laplacianRho: " << laplacianRho << "\tgradModGradRhoDotGradRho: " << gradModGradRhoDotGradRho <<"\ts: " << s << "\t F: " << F << "\tF': " << dFds << "\tF'': " << d2Fds2 << "\tvxc: " << vxc << "\t" << vxcLIBXC << "\tRel. Err.: " << std::abs((vxc-vxcLIBXC)/vxcLIBXC) << std::endl;
            //std::cout << "Point: (" << point[0] << ",\t" << point[1] << ",\t" << point[2] << ")\tvxc: " << vxc << "\t" << "\tvxcTest: " << vxcTest << "\tvxcLIBXC: " << vxcLIBXC << "\tVxc Error: " << std::abs((vxc-vxcLIBXC)) << "\tvxcTest Error: " << std::abs(vxcTest-vxcLIBXC) << std::endl;
        }
    }
}
