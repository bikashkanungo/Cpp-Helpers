#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>
#include <cmath>

#define SHIFT_VAL -1.0
#define AMP_VAL 0.5
#define YUKAWA_M 1.0
#define YUKAWA_AMP -3.0


double f1(const double x)
{
    if(x <= 0.0)
        return 0.0;

    else
        return exp(-1.0/x);

}

double f1Der(const double x)
{

    return f1(x)/(x*x);

}

double f2(const double x)
{

    return (f1(x)/(f1(x)+f1(1-x)));

}


double f2Der(const double x)
{

    if(fabs(x-0.0) < 1e-12 || fabs(1-x) < 1e-12)
        return 0.0;

    else
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


double getMollifierValue(const double x, const double r, const double d)
{

    const double y = Y(x,r,d);
    return pow(f2(y), 1.0);

}

double getMollifierDerivative(const double x, const double r, const double d)
{

    const double y = Y(x,r,d);
    return f2Der(y)*YDer(x,r,d);

}

double  getYukawaPotential(const double x)
{

    const double y = std::fabs(x);
    return YUKAWA_AMP*exp(-YUKAWA_M*y)/y;
}

double getVPSP(const double x, const double rc)
{

    const double f = getMollifierValue(x, rc, 2.0);
    const double a = getYukawaPotential(rc);
    return f*(SHIFT_VAL + AMP_VAL*cos(2.0*M_PI*x/(1.2*rc))) + (1.0-f)*a;
}


double poly(const double x, const std::vector<double> & nodes)
{
  double returnValue = 1.0;
  const int numNodes = nodes.size();
  for(unsigned int i = 0; i < numNodes; ++i)
    returnValue *= (x-nodes[i]);
  return returnValue;
}


int main()
{

    double rc;
    double h, xhi;
    int psiNumNodes;
    std::string outfileName;
    std::ofstream outfile;

    std::cout << "Enter cut-off radius: ";
    std::cin >> rc;
    std::cout << std::endl;
    std::cout << "Enter number nodes (excluding origin) for the wavefunction within cut-off radius: ";
    std::cin >> psiNumNodes;
    std::cout << std::endl;
    std::cout << "Enter grid-size: ";
    std::cin >> h;
    std::cout << std::endl;
    std::cout << "Enter x-max: ";
    std::cin >> xhi;
    std::cout << std::endl;
    std::cout << "Enter output filename: ";
    std::cin >> outfileName;
    std::cout << std::endl;

   outfile.open(outfileName.c_str());

   //psiNumNodes += 1;
   std::vector<double> psiNodes(psiNumNodes);
   psiNodes[0] = 0.0;
   for(unsigned int i = 1; i < psiNumNodes; ++i)
        psiNodes[i] = rc*((rand()+0.0)/(RAND_MAX + 1.0));
   for(unsigned int i = 0; i < psiNumNodes; ++i)
     std::cout << "Node[" << i << "]: " << psiNodes[i] << std::endl; 

    int numPoints = xhi/h;
    const double smoothFactor = 2.0;
    outfile << "# x \t PsiPseudo \t PsiAE \t VPSP \t VAE" << std::endl;
    for(unsigned int i = 0; i < numPoints; ++i)
    {

        const double x = i*h + 0.0001;
        const double f = getMollifierValue(x,rc,smoothFactor);
        const double f1 = getMollifierValue(x,rc-rc/5.0,smoothFactor);
        const double f2 = getMollifierValue(x,rc,0.5*smoothFactor);
        const double x0 = rc*(1.0+1.0/smoothFactor);
        const double psiPseudo = f*exp(-2.0*(x-x0)*(x-x0)) + (1.0-f)*exp(-std::fabs(x-rc));
        const double psiAE = 2.0*f1*exp(-std::fabs(x-rc)*3.0)*sin(x*2*M_PI*(psiNumNodes+0.5)/rc) + (1.0-f1)*psiPseudo;
        double vAE = -1.0/x;//getYukawaPotential(x);
        double vPSP = f2*getVPSP(x,rc) + (1.0-f2)*vAE;
        vAE = f*vAE + (1.0-f)*vPSP;
        vPSP = f2*vPSP + (1.0-f2)*vAE;
        vAE = -1.0/x;//f*getYukawaPotential(x) + (1.0-f)*(-3.0/x);//f*vAE + (1.0-f)*vPSP;
        outfile << x << "\t" << psiPseudo << "\t" << psiAE << "\t" << vPSP << "\t" << vAE << std::endl;
    }

}
