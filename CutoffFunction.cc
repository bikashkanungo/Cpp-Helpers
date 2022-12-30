#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#define a 1.0
#define b 1.0

    double
    getFermiDiracCutoff(const double radius, const double mollifierRadius)
    {
    
      //
      // Using a Frmi-Dirac distribution with kb*T = mu/10.0
      //
      double kbT = mollifierRadius/5.0;//std::max(mollifierRadius/5.0, 0.25);
      //double returnValue = 1.0/(exp((radius - mollifierRadius)/kbT) + 1.0);
      double returnValue = exp((mollifierRadius - radius)/kbT)/(exp((mollifierRadius - radius)/kbT) + 1.0);
    
      return returnValue;
    
    
    }
    
    double
    getFermiDiracCutoffDerivative(const double radius, const double mollifierRadius)
    {
    
      //
      // Using a Frmi-Dirac distribution with kb*T = mu/10.0
      //
      double kbT = mollifierRadius/5.0;//std::max(mollifierRadius/5.0, 0.25);
      //double returnValue = -1.0 * exp((radius - mollifierRadius)/kbT)/(kbT * pow(exp((radius - mollifierRadius)/kbT) + 1.0, 2.0));
      double returnValue = -1.0 * exp((mollifierRadius - radius)/kbT)/(kbT * pow(exp((mollifierRadius - radius)/kbT) + 1.0, 2.0));
    
      return returnValue;
    
    }

    double f1(const double x)
    {
      if(x <= 0.0)
        return 0.0;
      
      else
        return exp(-1.0/(a*x));
    
    }
    
    double f1Der(const double x)
    {
    
      return f1(x)/(a*x*x);
    
    }
    
    double f2(const double x)
    {
    
      return (f1(x)/(f1(x)+f1(1-x)));
    
    }
    
    
    double f2Der(const double x)
    {
    
      return ((f1Der(x)*f1(1-x) + f1(x)*f1Der(1-x))/(pow(f1(x) + f1(1-x), 2.0)));
    
    }
   
    double Y(const double x, const double r, const double d)
    {

      return  (1 - d*(x-r)/r);

    }

    double YDer(const double x, const double r, const double d)
    {

	return (-d/r);

    }

 
    double getSmoothCutoff(const double x, const double r, const double d)
    {
	    
      const double y = Y(x,r,d);
      return pow(f2(y), 1.0);
	          
    }
    
    double getSmoothCutoffDerivative(const double x, const double r, const double d)
    {
	    
      const double y = Y(x,r,d);
      return f2Der(y)*YDer(x,r,d);
		        
    }
  
    double logisticFunction(const double x, const double r, const double k)
    {

      const double exponent = -k*(x-r);
      if(exponent > log(1e15))
	return 1.0;
      else
	return 1.0/(1.0+exp(exponent));
    }
  
int main()
{

  double r;
  std::cout << "\nEnter cut-off: ";
  std::cin >> r;
  
  double d;
  std::cout << "\nEnter steepness for cutoff: ";
  std::cin >> d;

  double k;
  std::cout << "\nEnter steepness for logistic: ";
  std::cin >> k;
  
  double tol;
  std::cout << "\nEnter Rho tol: ";
  std::cin >> tol;

  std::ofstream outFile;
  outFile.open("SmoothCutoff.dat");

  const double start = 0.000;
  const double end = 20.0;
  const double delta = 0.001;

  double x = start;

  while(x <= end)
  {

    const double y = Y(x,r,d);
    outFile << x << "\t" << getSmoothCutoff(x,r,d) << "\t" << getSmoothCutoffDerivative(x,r,d) << "\t" << pow(f2(y), -0.5) << "\t" << getFermiDiracCutoff(x,r) << "\t" << getFermiDiracCutoffDerivative(x,r) << std::endl;
    x += delta;

  }

  std::ofstream outFileWeight;
  outFileWeight.open("RhoWeight.dat");
  x = start;
  while(x <= end)
  {

    const double rho = 100.0*exp(-2.0*x);
    const double numerator = (rho+tol);
    const double denominator = rho*getSmoothCutoff(x,r,d) + logisticFunction(x,r,k);
    outFileWeight << x << "\t" << rho << "\t" << numerator << "\t" << denominator << "\t" << logisticFunction(x,r,k) <<"\t" << pow(denominator/numerator,2.0) << std::endl;
    x += delta;

  }

}

    

