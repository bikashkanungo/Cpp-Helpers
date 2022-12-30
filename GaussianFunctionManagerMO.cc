#include "GaussianFunctionManagerMO.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <cerrno>
#include <cmath>
#include <cassert>
#include <iomanip>
//
//
//
/**
 * For the definition of the associated Legendre polynomials i.e. Plm() and their derivatives 
 * (as used for evaluating the real form of spherical harmonics and their derivatives) refer:
 * @article{bosch2000computation,
 *  	   title={On the computation of derivatives of Legendre functions},
 *     	   author={Bosch, W},
 *         journal={Physics and Chemistry of the Earth, Part A: Solid Earth and Geodesy},
 *         volume={25},
 *         number={9-11},
 *         pages={655--659},
 *         year={2000},
 *         publisher={Elsevier}
 *        }
 */

namespace dft {

#define ANGSTROM_TO_BOHR 1.889725989
#define DIST_TOL 1e-8
#define RADIUS_TOL 1e-15
#define POLAR_ANGLE_TOL 1e-12
#define EPS_MACHINE 1e-30
#define EXPONENT_TOL log(EPS_MACHINE)
  namespace
  {
	double doubleFactorial(int n)
	{
	  if (n == 0 || n==-1)
		return 1.0;
	  return n*doubleFactorial(n-2);
	}

	double dotProduct(const std::vector<double> & x,
		const std::vector<double> & y)
	{

	  double returnValue = 0.0;
	  for(unsigned int i = 0; i < x.size(); ++i)
		returnValue += x[i]*y[i];

	  return returnValue;
	}


	void
	  convertCartesianToSpherical(const std::vector<double> & x, 
		  double & r, 
		  double & theta, 
		  double & phi)
	  {

		r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

		if(r == 0)
		{

		  theta = 0.0;
		  phi = 0.0;

		}

		else
		{

		  theta = acos(x[2]/r);
		  //
		  // check if theta = 0 or PI (i.e, whether the point is on the Z-axis)
		  // If yes, assign phi = 0.0.
		  // NOTE: In case theta = 0 or PI, phi is undetermined. The actual value 
		  // of phi doesn't matter in computing the enriched function value or 
		  // its gradient. We assign phi = 0.0 here just as a dummy value
		  //
		  if(fabs(theta - 0.0) >= POLAR_ANGLE_TOL && fabs(theta - M_PI) >= POLAR_ANGLE_TOL)
			phi = atan2(x[1],x[0]);
		  else
			phi = 0.0;

		}

	  }


	template<typename T, std::size_t N>
	  std::size_t length_of( T const (&)[N] ) 
	  { 
		return N; 
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

		  const std::string message("Invalid number of FD points. Please enter number of FD points as 3, 5, 7, 9, 11 or 13.");
		  throw dft::InvalidArgument(message);
		}

		return returnValue;
	  }
	void initializeCartesianExponentsMoldenFormat(std::vector<std::vector<int> > & sCartExp,
		std::vector<std::vector<int> > & pCartExp,
		std::vector<std::vector<int> > & dCartExp,
		std::vector<std::vector<int> > & fCartExp,
		std::vector<std::vector<int> > & gCartExp)
	{

	  const char * pCartExpChars[] = {"x", "y", "z"}; 
	  const char * dCartExpChars[] = {"xx", "yy", "zz", "xy", "xz", "yz"}; 
	  const char * fCartExpChars[] = {"xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"}; 
	  const char * gCartExpChars[] = {"xxxx", "yyyy", "zzzz", "xxxy", "xxxz", "yyyx", "yyyz", "zzzx", "zzzy",
		"xxyy", "xxzz", "yyzz", "xxyz", "yyxz", "zzxy"}; 

	  int pSize = length_of(pCartExpChars);
	  int dSize = length_of(dCartExpChars);
	  int fSize = length_of(fCartExpChars);
	  int gSize = length_of(gCartExpChars);

	  assert(pSize == 3);
	  assert(dSize == 6);
	  assert(fSize == 10);
	  assert(gSize == 15);

	  std::vector<std::string> pCartExpStrings(pCartExpChars, pCartExpChars + pSize);
	  std::vector<std::string> dCartExpStrings(dCartExpChars, dCartExpChars + dSize);
	  std::vector<std::string> fCartExpStrings(fCartExpChars, fCartExpChars + fSize);
	  std::vector<std::string> gCartExpStrings(gCartExpChars, gCartExpChars + gSize);

	  sCartExp.resize(1, std::vector<int>(3,0));
	  pCartExp.resize(pSize, std::vector<int>(3));
	  dCartExp.resize(dSize, std::vector<int>(3));
	  fCartExp.resize(fSize, std::vector<int>(3));
	  gCartExp.resize(gSize, std::vector<int>(3));

	  for(unsigned int i = 0; i < pSize; ++i)
	  {
		pCartExp[i][0] = std::count(pCartExpStrings[i].begin(), pCartExpStrings[i].end(), 'x');    
		pCartExp[i][1] = std::count(pCartExpStrings[i].begin(), pCartExpStrings[i].end(), 'y');    
		pCartExp[i][2] = std::count(pCartExpStrings[i].begin(), pCartExpStrings[i].end(), 'z');    
	  }

	  for(unsigned int i = 0; i < dSize; ++i)
	  {
		dCartExp[i][0] = std::count(dCartExpStrings[i].begin(), dCartExpStrings[i].end(), 'x');    
		dCartExp[i][1] = std::count(dCartExpStrings[i].begin(), dCartExpStrings[i].end(), 'y');    
		dCartExp[i][2] = std::count(dCartExpStrings[i].begin(), dCartExpStrings[i].end(), 'z');    
	  }

	  for(unsigned int i = 0; i < fSize; ++i)
	  {
		fCartExp[i][0] = std::count(fCartExpStrings[i].begin(), fCartExpStrings[i].end(), 'x');    
		fCartExp[i][1] = std::count(fCartExpStrings[i].begin(), fCartExpStrings[i].end(), 'y');    
		fCartExp[i][2] = std::count(fCartExpStrings[i].begin(), fCartExpStrings[i].end(), 'z');    
	  }

	  for(unsigned int i = 0; i < gSize; ++i)
	  {
		gCartExp[i][0] = std::count(gCartExpStrings[i].begin(), gCartExpStrings[i].end(), 'x');    
		gCartExp[i][1] = std::count(gCartExpStrings[i].begin(), gCartExpStrings[i].end(), 'y');    
		gCartExp[i][2] = std::count(gCartExpStrings[i].begin(), gCartExpStrings[i].end(), 'z');    
	  }

	}

	void initializeCartesianExponentsNWChemFormat(std::vector<std::vector<int> > & sCartExp,
		std::vector<std::vector<int> > & pCartExp,
		std::vector<std::vector<int> > & dCartExp,
		std::vector<std::vector<int> > & fCartExp,
		std::vector<std::vector<int> > & gCartExp)
	{

	  for(unsigned int polyOrder = 0; polyOrder <= 4; ++polyOrder)
	  {          
		std::vector<std::vector<int> > expVec(0);
		for(int i = polyOrder; i >= 0; --i)
		{
		  for(int j = polyOrder; j >= 0; --j)
		  {
			for(int k = polyOrder; k >= 0; --k)
			{
			  if(i+j+k == polyOrder)
			  {
				std::vector<int> exponents(3);
				exponents[0] = i;
				exponents[1] = j;
				exponents[2] = k;
				expVec.push_back(exponents);
			  }
			}
		  }
		}

		if(polyOrder == 0)
		  sCartExp = expVec;
		if(polyOrder == 1)
		  pCartExp = expVec;
		if(polyOrder == 2)
		  dCartExp = expVec;
		if(polyOrder == 3)
		  fCartExp = expVec;
		if(polyOrder == 4)
		  gCartExp = expVec;

	  }
	}

	std::vector<std::vector<int> > 
	  getMoldenToNWChemCartExpMap(const std::vector<std::vector<std::vector<int> > > & cartExpMoldenVecs,
		  const std::vector<std::vector<std::vector<int> > > & cartExpNWChemVecs)
	  {

		const int numPoly = cartExpMoldenVecs.size();
		assert(numPoly == cartExpNWChemVecs.size());
		std::vector<std::vector<int> > returnValue(numPoly);
		for(unsigned int i = 0; i < numPoly; ++i)
		{
		  const std::vector<std::vector<int> > & icartExpMoldenVec = cartExpMoldenVecs[i];
		  const std::vector<std::vector<int> > & icartExpNWChemVec = cartExpNWChemVecs[i];

		  const int iNumExp = icartExpMoldenVec.size();
		  assert(iNumExp == icartExpNWChemVec.size());
		  returnValue[i].resize(iNumExp);
		  for(unsigned j = 0; j < iNumExp; ++j)
		  {
			const std::vector<int> moldenExp = icartExpMoldenVec[j];
			for(unsigned int k = 0; k < iNumExp; ++k)
			{
			  const std::vector<int> & nwchemExp = icartExpNWChemVec[k];
			  bool flag = true;
			  for(unsigned l = 0; l < nwchemExp.size(); ++l)
			  {
				if(nwchemExp[l] != moldenExp[l])
				{
				  flag = false;
				  break;
				}
			  }

			  if(flag)
			  {
				returnValue[i][j] = k;
			  }
			}
		  }
		}

		return returnValue;

	  }

	int factorial(int n)
	{

	  if(n==0) 
		return 1;
	  else 
		return n*factorial(n-1);

	}

	std::vector<double> getNormConsts(const double * alpha, const int l,const int L)
	{

	  std::vector<double> returnValue(L);
	  for(unsigned int i = 0; i < L; ++i)
	  {

		const double term1 = doubleFactorial(2*l+1)*sqrt(M_PI);
		const double term2 = pow(2.0,2*l+3.5)*pow(alpha[i],l+1.5);
		const double overlapIntegral = term1/term2;
		returnValue[i] = 1.0/sqrt(overlapIntegral);
	  }

	  return returnValue;

	}

	std::vector<std::vector<int> > generateCartesianExponents(const int polyOrder)
	{

	  std::vector<std::vector<int> > returnValue(0);
	  for(int i = polyOrder; i >= 0; --i)
	  {
		for(int j = polyOrder; j >= 0; --j)
		{
		  for(int k = polyOrder; k >= 0; --k)
		  {
			if(i+j+k == polyOrder)
			{
			  std::vector<int> exponents(3);
			  exponents[0] = i;
			  exponents[1] = j;
			  exponents[2] = k;
			  returnValue.push_back(exponents);
			}
		  }
		}
	  }

	  return returnValue;
	}


	bool isNumber(double &i, std::string s)
	{
	  char *end;
	  double  d;
	  errno = 0;
	  d = strtod(s.c_str(), &end);
	  if ((errno == ERANGE && d == std::numeric_limits<double>::max())) 
	  {
		return false;
	  }

	  if ((errno == ERANGE && d == std::numeric_limits<double>::min())) 
	  {
		return false;
	  }

	  if (s.size() == 0 || *end != 0) 
	  {
		return false;
	  }

	  i = d;
	  return true;
	}

	bool isInteger(int &i, std::string s, const int base = 10)
	{
	  char *end;
	  long  l;
	  errno = 0;
	  l = strtol(s.c_str(), &end, base);
	  if ((errno == ERANGE && l == std::numeric_limits<long>::max()) || l > std::numeric_limits<int>::max()) 
	  {
		return false;
	  }

	  if ((errno == ERANGE && l == std::numeric_limits<long>::min()) || l < std::numeric_limits<int>::min()) 
	  {
		return false;
	  }

	  if (s.size() == 0 || *end != 0) 
	  {
		return false;
	  }

	  i = l;
	  return true;
	}

	std::string trim(const std::string& str,
		const std::string& whitespace = " \t")
	{
	  std::size_t strBegin = str.find_first_not_of(whitespace);
	  if (strBegin == std::string::npos)
		return ""; // no content
	  std::size_t strEnd = str.find_last_not_of(whitespace);
	  std::size_t strRange = strEnd - strBegin + 1;
	  return str.substr(strBegin, strRange);
	}

	double getDistance(const double * x, const double * y)
	{
	  double r = 0.0;
	  for(unsigned int i = 0; i < 3; ++i)
		r += pow(x[i]-y[i],2.0);
	  return sqrt(r);
	}

	inline double Dm(const int m)
	{

	  if(m == 0)
		return 1.0/sqrt(2*M_PI);
	  else  
		return 1.0/sqrt(M_PI);
	}

	inline double Clm(const int l, const int m)
	{

	  assert(m >= 0);
	  assert(std::abs(m) <= l);
	  //const int modM = std::abs(m);
	  return sqrt(((2.0*l + 1)*boost::math::factorial<double>(l-m))/(2.0*boost::math::factorial<double>(l+m)));
	}

	double Qm(const int m, const double phi)
	{

	  if(m > 0)
		return cos(m*phi);
	  if(m == 0)
		return 1.0;
	  if(m < 0)
		return sin(std::abs(m)*phi);
	}

	double dQmDPhi(const int m, const double phi)
	{

	  if(m > 0)
		return -m*sin(m*phi);
	  if(m == 0)
		return 0.0;
	  if(m < 0)
		return std::abs(m)*cos(std::abs(m)*phi);

	}

	double Plm(const int l, const int m, const double x)
	{
	  if(std::abs(m) > l)
		return 0.0;
	  else
		//
		// NOTE: Multiplies by {-1}^m to remove the 
		// implicit Condon-Shortley factor in the associated legendre 
		// polynomial implementation of boost
		// This is done to be consistent with the QChem's implementation
		return pow(-1.0,m)*boost::math::legendre_p(l, m, x);
	}


	double dPlmDTheta(const int l, const int m, const double theta)
	{

	  const double cosTheta = cos(theta);

	  if(std::abs(m) > l)
		return 0.0;

	  else if(l == 0)
		return 0.0;

	  else if(m < 0)
	  {
		const int modM = std::abs(m);
		const double factor = pow(-1,m)*boost::math::factorial<double>(l-modM)/boost::math::factorial<double>(l+modM);
		return factor*dPlmDTheta(l, modM, theta);
	  }

	  else if(m == 0)
	  {

		return -1.0*Plm(l,1,cosTheta);
	  }

	  else if(m == l)
		return l*Plm(l,l-1, cosTheta);

	  else
	  {
		const double term1  = (l+m)*(l-m+1)*Plm(l, m-1, cosTheta);
		const double term2 = Plm(l, m+1, cosTheta);
		return 0.5*(term1-term2);
	  }

	}


	double d2PlmDTheta2(const int l, const int m, const double theta)
	{

	  const double cosTheta = cos(theta);

	  if(std::abs(m) > l)
		return 0.0;

	  else if(l == 0)
		return 0.0;

	  else if(m < 0)
	  {
		const int modM = std::abs(m);
		const double factor = pow(-1,m)*boost::math::factorial<double>(l-modM)/boost::math::factorial<double>(l+modM);
		return factor*d2PlmDTheta2(l, modM, theta);
	  }

	  else if(m == 0)
		return -1.0*dPlmDTheta(l, 1, theta);

	  else if(m == l)
		return l*dPlmDTheta(l,l-1, theta);

	  else
	  {
		double term1 = (l+m)*(l-m+1)*dPlmDTheta(l, m-1, theta);
		double term2 = dPlmDTheta(l, m+1, theta);
		return 0.5*(term1-term2);
	  }
	}


	double gfRadialPart(const double r, const int l, const double alpha)
	{

	  return pow(r,l)*exp(-alpha*r*r);
	}


	double gfRadialPartDerivative(const double r, const double alpha, const int l, const int derOrder)
	{
	  if(derOrder == 0 && l >= 0)
		return pow(r,l)*exp(-alpha*r*r);
	  else if(derOrder == 0 && l < 0)
		return 0.0;
	  else
		return l*gfRadialPartDerivative(r,alpha,l-1,derOrder-1) - 2*alpha*gfRadialPartDerivative(r,alpha,l+1,derOrder-1);
	}

	double getLimitingValueLaplacian(const int l, const int m, const double theta)
	{
	  double returnValue = 0.0;

	  if(std::fabs(theta-0.0) < POLAR_ANGLE_TOL)
	  {

		if(m == 0)
		  returnValue = -0.5*l*(l+1);
		if(m == 2)
		  returnValue = 0.25*(l-1)*l*(l+1)*(l+2);
	  }

	  if(std::fabs(theta-M_PI) < POLAR_ANGLE_TOL)
	  {

		if(m == 0)
		  returnValue = -0.5*l*(l+1)*pow(-1.0,l);
		if(m == 2)
		  returnValue = 0.25*(l-1)*l*(l+1)*(l+2)*pow(-1.0,l);;
	  }

	  return returnValue;
	}


	double evaluateBasisValue(const dft::GaussianFunctionManagerMO::basis * b, const double * x)
	{
	  const int L = b->L;
	  const double * R = b->origin;
	  const int l = b->l;
	  const int m = b->m;

	  std::vector<double> dx(3,0.0);
	  for(unsigned int i = 0; i < 3; ++i)
		dx[i] = x[i] - R[i];

	  double r, theta, phi;
	  convertCartesianToSpherical(dx, r, theta, phi);

	  double returnValue = 0.0;
	  for(unsigned int i = 0; i < L; ++i)
	  {
		const double alphaVal = b->alpha[i];
		const double cVal = b->c[i];
		const double norm = b->normConsts[i];
		returnValue += cVal*norm*gfRadialPart(r, l, alphaVal);
	  }

	  const int modM = std::abs(m);
	  const double C = Clm(l,modM)*Dm(m);
	  const double cosTheta = cos(theta);
	  const double P = Plm(l,modM,cosTheta);
	  const double Q = Qm(m,phi);

	  returnValue *= C*P*Q;
	  returnValue *= b->basisNormConst;
	  return returnValue;
	}


	std::vector<double> getSphericalGaussianGradient(const std::vector<double> & x,
		const int l,
		const int m,
		const double alpha)
	{
	  double r, theta, phi;
	  convertCartesianToSpherical(x, r, theta, phi);
	  std::vector<double> returnValue(3);
	  const int modM = std::abs(m);
	  const double C = Clm(l,modM)*Dm(m);
	  if(r < RADIUS_TOL)
	  {

		if(l==1)
		{
		  if(m==-1)
		  {
			returnValue[0] = 0.0;
			returnValue[1] = C;
			returnValue[2] = 0.0;
		  }

		  if(m==0)
		  {
			returnValue[0] = 0.0;
			returnValue[1] = 0.0;
			returnValue[2] = C;
		  }

		  if(m==1)
		  {
			returnValue[0] = C;
			returnValue[1] = 0.0;
			returnValue[2] = 0.0;
		  }

		}

		else
		{
		  returnValue[0] = 0.0;
		  returnValue[1] = 0.0;
		  returnValue[2] = 0.0;
		}

	  }

	  else if(std::fabs(theta-0.0) < POLAR_ANGLE_TOL)
	  {
		const double R = gfRadialPart(r,l,alpha);
		const double dRDr = gfRadialPartDerivative(r, alpha, l, 1);
		const double cosTheta = cos(theta);
		const double P = Plm(l, modM, cosTheta);
		if(m == 0)
		{
		  returnValue[0] = 0.0;
		  returnValue[1] = 0.0;
		  returnValue[2] = C*dRDr*P*cosTheta;

		}

		else if(m == 1)
		{

		  returnValue[0] = C*(R/r)*l*(l+1)/2.0;
		  returnValue[1] = 0.0;
		  returnValue[2] = 0.0;


		}

		else if(m == -1)
		{
		  returnValue[0] = 0.0;
		  returnValue[1] = C*(R/r)*l*(l+1)/2.0;
		  returnValue[2] = 0.0;

		}

		else
		{
		  returnValue[0] = 0.0;
		  returnValue[1] = 0.0;
		  returnValue[2] = 0.0;


		}

	  }

	  else if(std::fabs(theta-M_PI) < POLAR_ANGLE_TOL)
	  {
		const double R = gfRadialPart(r,l,alpha);
		const double dRDr = gfRadialPartDerivative(r, alpha, l, 1);
		const double cosTheta = cos(theta);
		const double P = Plm(l, modM, cosTheta);
		if(m == 0)
		{
		  returnValue[0] = 0.0;
		  returnValue[1] = 0.0;
		  returnValue[2] = C*dRDr*P*cosTheta;

		}

		else if(m == 1)
		{

		  returnValue[0] = C*(R/r)*l*(l+1)/2.0*pow(-1,l+1);
		  returnValue[1] = 0.0;
		  returnValue[2] = 0.0;


		}

		else if(m == -1)
		{
		  returnValue[0] = 0.0;
		  returnValue[1] = C*(R/r)*l*(l+1)/2.0*pow(-1,l+1);
		  returnValue[2] = 0.0;

		}

		else
		{
		  returnValue[0] = 0.0;
		  returnValue[1] = 0.0;
		  returnValue[2] = 0.0;
		}

	  }

	  else
	  {

		const double R = gfRadialPart(r,l,alpha);
		const double dRDr = gfRadialPartDerivative(r, alpha, l, 1);
		const double cosTheta = cos(theta);
		const double P = Plm(l, modM, cosTheta);
		const double dPDTheta = dPlmDTheta(l, modM, theta);
		const double Q = Qm(m,phi);
		const double dQDPhi = dQmDPhi(m,phi);      
		double jacobianInverse[3][3];
		jacobianInverse[0][0] = sin(theta)*cos(phi); jacobianInverse[0][1] = cos(theta)*cos(phi)/r; jacobianInverse[0][2] = -1.0*sin(phi)/(r*sin(theta));
		jacobianInverse[1][0] = sin(theta)*sin(phi); jacobianInverse[1][1] = cos(theta)*sin(phi)/r; jacobianInverse[1][2] = cos(phi)/(r*sin(theta));
		jacobianInverse[2][0] = cos(theta); 	   jacobianInverse[2][1] = -1.0*sin(theta)/r;     jacobianInverse[2][2] = 0.0;

		double partialDerivatives[3];
		partialDerivatives[0] = dRDr*P*Q;
		partialDerivatives[1] = R*dPDTheta*Q;
		partialDerivatives[2] = R*P*dQDPhi;
		for(unsigned int i = 0; i < 3; ++i)
		{
		  returnValue[i] = C*(jacobianInverse[i][0]*partialDerivatives[0]+
			  jacobianInverse[i][1]*partialDerivatives[1]+
			  jacobianInverse[i][2]*partialDerivatives[2]);
		}

	  }

	  return returnValue;

	}


	double getSphericalGaussianLaplacian(const std::vector<double> & x,
		const int l,
		const int m,
		const double alpha)
	{

	  double r, theta, phi;
	  convertCartesianToSpherical(x, r, theta, phi);
	  double returnValue = 0.0;
	  if(r < RADIUS_TOL)
	  {
		const int modM = std::abs(m);
		const double C = Clm(l,modM)*Dm(m);
		if(l == 0)
		  returnValue= -C*6.0*alpha;
		else
		  returnValue = 0.0;
	  }

	  else
	  {

		const int modM = std::abs(m);
		const double C = Clm(l,modM)*Dm(m);
		const double cosTheta = cos(theta);
		const double sinTheta = sin(theta);
		const double R = gfRadialPart(r, l , alpha);
		const double dRdr = gfRadialPartDerivative(r, alpha, l, 1);
		const double d2Rdr2 = gfRadialPartDerivative(r, alpha, l, 2);
		const double P = Plm(l,modM,cosTheta);
		const double Q = Qm(m,phi);
		const double term1 = C*P*Q*(2.0*dRdr/r + d2Rdr2);

		if(std::fabs(theta-0.0) < POLAR_ANGLE_TOL || std::fabs(theta-M_PI) < POLAR_ANGLE_TOL)
		{
		  const double limitingVal = getLimitingValueLaplacian(l, modM, theta);
		  const double term2 = C*(R/(r*r))*Q*(limitingVal+limitingVal);
		  const double term3 = -C*m*m*(R/(r*r))*Q*(limitingVal/2.0);
		  returnValue = (term1 + term2 + term3);
		}

		else
		{

		  const double a = dPlmDTheta(l, modM, theta);
		  const double b = d2PlmDTheta2(l, modM, theta);
		  const double term2 = C * (R/(r*r)) * Q * 
			((cosTheta/sinTheta)*a + b);
		  const double term3 = -C*m*m*(R/(r*r))*Q*P/(sinTheta*sinTheta);
		  returnValue = term1 + term2 + term3;

		}

	  }

	  return returnValue;
	}

	double
	  getBasisHigherOrderDerivativeFD(const double * point,
		  const dft::GaussianFunctionManagerMO::basis * basis,                    
		  std::vector<int> & indices,
		  const int numFDPoints,
		  const double h)
	  {

		int numIndices = indices.size();
		double returnValue = 0.0;
		if(numIndices == 0)
		{
		  returnValue = evaluateBasisValue(basis, point);
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
			std::vector<double> FDPoint(3);
			FDPoint[0] = point[0];
			FDPoint[1] = point[1];
			FDPoint[2] = point[2];
			int shiftIndex = (int)(iPoint-numFDPoints/2);
			FDPoint[currentIndex] += shiftIndex*h;
			const double factor = coeffs[iPoint]/h;
			returnValue += factor*getBasisHigherOrderDerivativeFD(&FDPoint[0],
				basis,
				indicesNext,
				numFDPoints,
				h);
		  }

		  return returnValue;
		}
	  }

	double trapezoidal3d(const double L,
		const double h,
		std::vector<const dft::GaussianFunctionManagerMO::basis *> funcs)
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

	void
	  readAtomicCoordsAndBasisFileNames(std::vector<std::vector<double> > & atomicCoords,
		  std::vector<std::string> & basisFileNames,
		  std::string fileName)
	  {

		std::ifstream readFile;

		//
		// string to read line
		//
		std::string readLine;

		readFile.open(fileName.c_str());
		assert(readFile.is_open()); 
		while(std::getline(readFile,readLine))
		{
		  if(readLine.empty())
		  {
			std::string message("Empty or invalid line while reading atomic coordinates in GaussianFunctionManagerMO");
			throw InvalidArgument(message);
		  }
		  std::istringstream lineString(readLine);
		  std::string word;
		  unsigned int count = 0;
		  std::vector<double> coord(0);
		  while(lineString >> word)
		  {
			if(count >= 5)
			{
			  std::string message("Invalid column entry while reading atomic coordinates in GaussianFunctionManagerMO. "
				  "Number of columns must be 5.");
			  throw InvalidArgument(message);
			}

			std::string wordTrimmed = trim(word);

			if(wordTrimmed.empty())
			{
			  std::string message("Empty column entry while reading atomic coordinates in GaussianFunctionManagerMO");
			  throw InvalidArgument(message);
			}

			if(count >=1 && count <= 3)
			{
			  double val;
			  if(isNumber(val, wordTrimmed))
				coord.push_back(val);
			  else
			  {
				std::string message("Coordinate entry in the atomic coordinates in GaussianFunctionManagerMO is not a number");
				throw InvalidArgument(message);

			  }
			}

			else if(count == 4)
			  basisFileNames.push_back(wordTrimmed);

			count++;

		  }

		  atomicCoords.push_back(coord);
		}

		assert(atomicCoords.size() == basisFileNames.size());
		readFile.close();
	  }

	void
	  readGaussianFiles(std::vector<dft::GaussianFunctionManagerMO::contractedGaussian *> & atomicContractedGaussians, const std::string fileName)
	  {

		std::ifstream readFile;

		//
		// string to read line
		//
		std::string readLine;

		readFile.open(fileName.c_str());
		assert(readFile.is_open()); 
		//
		// ignore the first line
		//
		std::getline(readFile, readLine);
		while(std::getline(readFile,readLine))
		{
		  std::istringstream lineString(readLine);
		  std::string word;
		  unsigned int count = 0;
		  while(lineString >> word)
		  {
			if(count >= 2)
			{   
			  count++;
			  continue;
			}
			if(!word.empty())
			  count++;
			//
			// check if it's a valid string
			// i.e., it contains one of the following string:
			// "S", "SP", "SPD", SPDF" ...
			// 
			std::size_t pos = word.find_first_not_of("SPDFGHspdfgh");
			if(pos == std::string::npos)
			{
			  std::string lChars = trim(word);
			  const int numLChars = lChars.size();
			  std::string strNContracted;
			  lineString >> strNContracted;
			  if(!strNContracted.empty())
				count++;
			  int nContracted;
			  if(isInteger(nContracted,strNContracted))
			  {
				double alpha[nContracted];
				double c[nContracted][numLChars];
				for(unsigned int i = 0; i < nContracted; ++i)
				{

				  if(std::getline(readFile,readLine))
				  {
					if(readLine.empty())
					  std::cout << "Undefined behavior in gaussian file: Empty line found" << std::endl;

					std::istringstream lineContracted(readLine);
					lineContracted >> word;
					if(!isNumber(alpha[i],word))
					  std::cout << "Undefined value read for the exponent in file: " << fileName << " " << word << std::endl;
					for(unsigned int j = 0; j < numLChars; ++j)
					{
					  lineContracted >> word;
					  if(!isNumber(c[i][j],word))
						std::cout << "Undefined value read for the coefficient in file: " << fileName << std::endl;
					}

				  }
				  else
					std::cout << "Undefined row for the contracted basis in file: " << fileName << std::endl;
				}

				for(unsigned int j = 0; j < numLChars; ++j)
				{
				  dft::GaussianFunctionManagerMO::contractedGaussian * cg = new dft::GaussianFunctionManagerMO::contractedGaussian;
				  cg->L = nContracted;
				  cg->alpha = new double[nContracted];
				  cg->c = new double[nContracted];
				  cg->lquantum = lChars.at(j);
				  for(unsigned int i = 0; i < nContracted; ++i)
				  {
					cg->alpha[i] = alpha[i];
					cg->c[i] = c[i][j];
				  }
				  atomicContractedGaussians.push_back(cg);
				}
			  }

			  else
				std::cout << "Undefined behavior: Number of contracted Gaussians is not an integer in " << fileName << std::endl;

			}

			else
			  std::cout << "Undefined L character(s) for the contracted basis read in file: " << fileName << std::endl;

		  }
		}

		readFile.close();
	  }

	void readMatrix(std::vector<std::vector<double> > & densityMat, const std::string fileName)
	{    
	  std::ifstream readFile;

	  //
	  // string to read line
	  //
	  std::string readLine;
	  readFile.open(fileName.c_str());
	  assert(readFile.is_open());
	  while(getline(readFile, readLine))
	  {

		if(!readLine.empty())
		{
		  std::vector<double> rowVals(0);
		  std::istringstream lineString(readLine);
		  std::string word;
		  while(lineString >> word)
		  {
			double val;
			if(isNumber(val,word))
			  rowVals.push_back(val);
			else
			  std::cout << "Undefined behavior: Value read in density matrix is not a number" << std::endl;
		  }
		  densityMat.push_back(rowVals);
		}
		else
		  std::cout << "Empty line found in density matrix file: " << fileName << std::endl;
	  }

	  readFile.close();

	}

  }

  //
  // default Constructor
  //
  GaussianFunctionManagerMO::GaussianFunctionManagerMO():
	d_atomicCoords(0),
	d_basisFileNames(0),
	d_basisFunctions(0),
	d_contractedGaussians(0)
  {
	std::string atomicFileName("AtomicCoords");
	std::string MOFileName("MOs");
	readAtomicCoordsAndBasisFileNames(d_atomicCoords, d_basisFileNames, atomicFileName);

	//
	// Convert from angstrom to bohr
	//
	const int numAtoms = d_atomicCoords.size();
	for(unsigned int i = 0; i < numAtoms; ++i)
	{
	  for(unsigned int j = 0; j < d_atomicCoords[i].size(); ++j)
		d_atomicCoords[i][j] *= ANGSTROM_TO_BOHR;
	}

	for(unsigned int i = 0; i < d_basisFileNames.size(); ++i)
	  d_uniqueBasisFileNames.insert(d_basisFileNames[i]);


	const int numUniqueBasisFiles = d_uniqueBasisFileNames.size();
	d_contractedGaussians.resize(numUniqueBasisFiles);
	unsigned int i = 0;
	for(std::set<std::string>::const_iterator iter = d_uniqueBasisFileNames.begin();
		iter!= d_uniqueBasisFileNames.end(); iter++)
	{
	  std::vector<dft::GaussianFunctionManagerMO::contractedGaussian *> & atomicContractedGaussians = 
		d_contractedGaussians[i];
	  std::string fileName = *iter;
	  readGaussianFiles(atomicContractedGaussians, fileName);
	  i++;
	}

	for(unsigned int i = 0; i < numUniqueBasisFiles; ++i)
	{
	  std::vector<dft::GaussianFunctionManagerMO::contractedGaussian *> & atomicContractedGaussians = 
		d_contractedGaussians[i];
	  const int numContractedGaussians = atomicContractedGaussians.size();
	  for(unsigned int j = 0; j < numContractedGaussians; ++j)
	  {
		const dft::GaussianFunctionManagerMO::contractedGaussian * cg = atomicContractedGaussians[j];
		const int L = cg->L;
		const char lquantum = cg->lquantum;
	  }
	}

	d_MOMat.resize(0);
	readMatrix(d_MOMat, MOFileName);
	std::cout << "Read MOs" << std::endl;
	for(unsigned int iAtom = 0; iAtom < d_atomicCoords.size(); ++iAtom)
	{
	  const std::string basisFileName = d_basisFileNames[iAtom];
	  std::set<std::string>::const_iterator iter = d_uniqueBasisFileNames.find(basisFileName);
	  unsigned int index = std::distance(d_uniqueBasisFileNames.begin(), iter);
	  std::vector<dft::GaussianFunctionManagerMO::contractedGaussian *> & atomicContractedGaussians = 
		d_contractedGaussians[index];
	  unsigned int numContractedGaussians = atomicContractedGaussians.size();
	  for(unsigned int j = 0; j < numContractedGaussians; ++j)
	  {
		const dft::GaussianFunctionManagerMO::contractedGaussian  * cg = atomicContractedGaussians[j];
		int polyOrder;
		char lquantum = cg->lquantum;

		if(lquantum == 'S'|| lquantum == 's')
		  polyOrder = 0;

		else if(lquantum == 'P'|| lquantum == 'p')
		  polyOrder = 1;

		else if(lquantum == 'D' || lquantum == 'd')
		  polyOrder = 2;

		else if(lquantum == 'F' || lquantum == 'f')
		  polyOrder = 3;

		else if(lquantum == 'G' || lquantum == 'g')
		  polyOrder = 4;

		else if(lquantum == 'H' || lquantum == 'h')
		  polyOrder = 5;

		else
		{
		  std::string message("Invalid L character detected while reading in GaussianFunctionManagerMO");
		  throw InvalidArgument(message);
		}

		const int nContracted = cg->L;

		//
		// NOTE: QChem even in the spherical form uses cartesian form for the s and p orbitals.
		// The ordering for cartesian orbitals are lexicographic - i.e., for p orbitals it's ordered
		// as x,y,z. While in the spherical form if one uses -l to +l ordering for the m quantum numbers for l=1 (p orbital), 
		// then it translates to ordering the p orbitals as y,z,x. 
		// To be consistent with QChem's ordering for p orbital, we do it in an 
		// ad-hoc manner by ordering the m quantum numbers as 1,-1,0 for l=1 (p orbital).
		//
		if(polyOrder == 1)
		{

		  int qChem_p[] = {1,-1,0};
		  for(unsigned int iM = 0; iM < 3; ++iM)
		  {

			dft::GaussianFunctionManagerMO::basis * b = new dft::GaussianFunctionManagerMO::basis;
			b->L = nContracted;
			b->alpha = cg->alpha;
			b->c = cg->c;
			b->l = polyOrder;
			b->m = qChem_p[iM];

			b->normConsts = new double[nContracted];
			std::vector<double> normConstsTmp = getNormConsts(b->alpha, b->l, nContracted);
			for(unsigned int iContracted = 0; iContracted < nContracted; ++iContracted)
			  b->normConsts[iContracted] = normConstsTmp[iContracted];

			b->origin = &d_atomicCoords[iAtom][0];
			// Set the basis normalization factor to 1.0
			b->basisNormConst = 1.0;
			d_basisFunctions.push_back(b);
		  }

		}

		else
		{
		  for(int m = -polyOrder; m <= polyOrder; ++m)
		  {
			dft::GaussianFunctionManagerMO::basis * b = new dft::GaussianFunctionManagerMO::basis;
			b->L = nContracted;
			b->alpha = cg->alpha;
			b->c = cg->c;
			b->l = polyOrder;
			b->m = m;

			b->normConsts = new double[nContracted];
			std::vector<double> normConstsTmp = getNormConsts(b->alpha, b->l, nContracted);
			for(unsigned int iContracted = 0; iContracted < nContracted; ++iContracted)
			  b->normConsts[iContracted] = normConstsTmp[iContracted];

			b->origin = &d_atomicCoords[iAtom][0];
			// Set the basis normalization factor to 1.0
			b->basisNormConst = 1.0;
			d_basisFunctions.push_back(b);

		  }
		}

	  }

	}

	unsigned int numBasis = d_basisFunctions.size();
	std::cout << "Number basis: " << d_basisFunctions.size() << std::endl;

	for(unsigned int i = 0; i < numBasis; ++i)
	{
	  d_basisFunctions[i]->basisNormConst = 1.0; //1.0/sqrt(SMatEvaluated[i][i]);
	}

	outfile.open("GaussianTestData");
	outfile << std::setprecision(16) << std::endl;

	for(unsigned int i = 0; i < numBasis; ++i)
	{
	  outfile << "Info for basis: " << i << std::endl;
	  const int L = d_basisFunctions[i]->L; 
	  outfile << "(l,m): (" << d_basisFunctions[i]->l << "," << d_basisFunctions[i]->m << ")\t" 
		<< "NormConsts: ";
	  for(unsigned int j = 0; j < L; ++j)
	  {

		outfile << d_basisFunctions[i]->normConsts[j] << "\t";
	  }

	  outfile << "Coeffs: ";
	  for(unsigned int j = 0; j < L; ++j)
	  {

		outfile << d_basisFunctions[i]->c[j] << "\t";
	  }

	  outfile << std::endl;
	}


	outfile << "\nPrinting Gaussian basis info\n";
	for(unsigned int i = 0; i < numBasis; ++i)
	{   
	  outfile << "\nBasis Id: " << i << "\tOrigin\t" << d_basisFunctions[i]->origin[0] << "\t"
		<< d_basisFunctions[i]->origin[1] << "\t" << d_basisFunctions[i]->origin[2] << std::endl;
	  outfile << "Alpha\tC" << std::endl;
	  int L = d_basisFunctions[i]->L;
	  for(unsigned int j = 0; j < L; ++j)
		outfile << d_basisFunctions[i]->alpha[j]/**pow(ANGSTROM_TO_BOHR,2.0)*/ << "\t" << d_basisFunctions[i]->c[j] << std::endl;
	}

	outfile.close();

  }

  //
  //Destructor
  //
  GaussianFunctionManagerMO::~GaussianFunctionManagerMO()
  {
	unsigned int numAtoms = d_contractedGaussians.size();
	for(unsigned int i = 0; i < numAtoms; ++i)
	{
	  std::vector<dft::GaussianFunctionManagerMO::contractedGaussian *> & atomicContractedGaussians = 
		d_contractedGaussians[i];
	  const int numAtomicContractedGaussians = atomicContractedGaussians.size();
	  for(unsigned int j = 0; j < numAtomicContractedGaussians; ++j)
	  {
		dft::GaussianFunctionManagerMO::contractedGaussian  * cg = atomicContractedGaussians[j];
		delete cg->alpha;
		delete cg->c;
		delete cg;
	  }
	}

	unsigned int numBasis = d_basisFunctions.size();
	for(unsigned int i = 0; i < numBasis; ++i)
	{
	  dft::GaussianFunctionManagerMO::basis * b = d_basisFunctions[i];
	  delete b->normConsts;
	  delete b;
	}
  }


  int 
	GaussianFunctionManagerMO::getNumberMOs()
	{

	  return d_MOMat[0].size();
	}
  //
  // get density value
  //
  double
	GaussianFunctionManagerMO::getMOValue(const double * x, const int MOId)
	{
	  const int numBasis = d_basisFunctions.size();
	  double val = 0.0;
	  for(unsigned int i = 0; i < numBasis; ++i)
	  {
		val += d_MOMat[i][MOId]*evaluateBasisValue(d_basisFunctions[i], x);
	  }

	  return val;
	}


  //
  // get basis value
  //
  double
	GaussianFunctionManagerMO::getBasisFunctionValue(const int basisId, const double * x)
	{
	  return evaluateBasisValue(d_basisFunctions[basisId], x);
	}

  //
  // get gradient of the basis
  //
  std::vector<double>
	GaussianFunctionManagerMO::getBasisFunctionGradient(const int basisId, const double * x)
	{
	  const dft::GaussianFunctionManagerMO::basis * b = 
		d_basisFunctions[basisId];
	  const int L = b->L;
	  const double * R = b->origin;
	  const int l = b->l;
	  const int m = b->m;

	  std::vector<double> dx(3,0.0);
	  for(unsigned int iCart = 0; iCart < 3; ++iCart)
		dx[iCart] = x[iCart] - R[iCart];

	  std::vector<double> returnValue(3,0.0);
	  for(unsigned int i = 0; i < L; ++i)
	  {

		const double alphaVal = b->alpha[i];
		const double cVal = b->c[i];
		const double norm = b->normConsts[i];

		std::vector<double> gradientPrimitiveGaussian = 
		  getSphericalGaussianGradient(dx, l, m, alphaVal);

		for(unsigned int iCart = 0; iCart < 3; ++iCart)
		  returnValue[iCart] += cVal*norm*gradientPrimitiveGaussian[iCart];
	  }


	  for(unsigned int iCart = 0; iCart < 3; ++iCart)
		returnValue[iCart] *= b->basisNormConst;

	  return returnValue;
	}


  //
  // get gradient of the basis
  //
  std::vector<double> 
	GaussianFunctionManagerMO::getBasisFunctionDoubleDerivatives(const int basisId, 
		const double * x)
	{

	  const int finiteDifferenceOrder = 13;
	  const double finiteDifferenceSpacing = 0.001;
	  std::vector<double> returnValue(9,0.0);
	  for(unsigned int iDim = 0; iDim < 3; ++iDim)
	  {
		for(unsigned int jDim = 0; jDim < 3; ++jDim)
		{
		  std::vector<int> indices(2);
		  indices[0] = iDim;
		  indices[1] = jDim;
		  int logically2DIndex = iDim*3 + jDim;
		  returnValue[logically2DIndex] = 
			getBasisHigherOrderDerivativeFD(x,
				d_basisFunctions[basisId],
				indices,
				finiteDifferenceOrder,
				finiteDifferenceSpacing);

		}

	  }

	  return returnValue;
	}

  //
  // get gradient of the basis
  //
  double GaussianFunctionManagerMO::getBasisFunctionLaplacian(const int basisId, const double * x)
  {
	const dft::GaussianFunctionManagerMO::basis * b = 
	  d_basisFunctions[basisId];
	const int L = b->L;
	const double * R = b->origin;
	const int l = b->l;
	const int m = b->m;

	std::vector<double> dx(3,0.0);
	for(unsigned int iCart = 0; iCart < 3; ++iCart)
	  dx[iCart] = x[iCart] - R[iCart];

	double returnValue = 0.0;
	for(unsigned int i = 0; i < L; ++i)
	{

	  const double alphaVal = b->alpha[i];
	  const double cVal = b->c[i];
	  const double norm = b->normConsts[i];

	  double laplacianPrimitiveGaussian = 
		getSphericalGaussianLaplacian(dx, l, m, alphaVal);

	  returnValue += cVal*norm*laplacianPrimitiveGaussian;
	}


	returnValue *= b->basisNormConst;

	return returnValue;
  }

  //
  // get number of basis functions
  //
  int GaussianFunctionManagerMO::getNumberBasisFunctions()
  {
	return d_basisFunctions.size();
  }

  std::vector<std::vector<double> > 
	GaussianFunctionManagerMO::getEvaluatedSMat()
	{

	  int meshId = 0;

	  //
	  // get mesh manager
	  //
	  MeshManager & meshManager = MeshManagerSingleton::getInstance();

	  //
	  // get QuadratureRuleManager
	  //
	  QuadratureRuleManager & quadratureRuleManager = QuadratureRuleManagerSingleton::getInstance();        

	  //
	  // get handle to FieldQuadratureTypeManager
	  //
	  FieldQuadratureTypeManager & fieldQuadratureTypeManager = FieldQuadratureTypeManagerSingleton::getInstance();


	  //
	  // Get the quadratureType for the fieldId
	  //

	  QuadratureRuleManager::QuadratureNameId quadratureType = fieldQuadratureTypeManager.getFieldQuadratureType(dft::ArrayNameManager::RHO);  

	  //
	  // get handle to Adaptive quadrature rule container
	  //
	  const QuadratureRuleContainer & quadratureRuleContainer = quadratureRuleManager.getQuadratureRuleContainer(quadratureType);

	  //
	  // get the number of elements in the mesh
	  //
	  const int numberElements = meshManager.getNumberElements(meshId);

	  //
	  // instantiate return value by getting the QuadratureValuesContainer associated with quadratureValuesManager
	  // FIXME: quadrature id used is 0
	  QuadratureValuesContainer<DoubleVector> quadValues(meshId,
		  0,
		  quadratureType,
		  1, //numberComponents
		  0.0);

	  const int numBasis = d_basisFunctions.size();
	  std::vector<std::vector<double> > SMat(numBasis, std::vector<double>(numBasis));
	  for(unsigned int i = 0; i < numBasis; ++i)
	  {
		for(unsigned int j = 0; j < numBasis; ++j)
		{

		  //
		  // iterate over elements
		  //
		  for (vtkIdType iElem = 0; iElem < numberElements; ++iElem) {


			//
			// get handle to the quadrature rule for the element
			//
			const QuadratureRule & quadratureRule = 
			  quadratureRuleContainer.getElementQuadratureRule(iElem,
				  0,
				  meshId);

			//
			// get the number of quad points
			//
			const QuadratureRule::quad_point_size_type numberQuadraturePoints = 
			  quadratureRule.getNumberPoints();


			//
			// get the current quad values for the element
			//
			const DoubleVector & elementQuadValuesCurrent = quadValues.getElementValues(iElem);

			//
			// copy current element quad values to a temporary storage
			//
			DoubleVector elementQuadValues = elementQuadValuesCurrent;

			//
			// compute value and gradient at quad points
			//
			for (int iQuadPoint = 0; 
				iQuadPoint < numberQuadraturePoints;
				++iQuadPoint) {

			  //
			  // get quad point coordinates
			  //
			  const Point & quadraturePoint = 
				quadratureRule.getGlobalPoint(iQuadPoint);

			  elementQuadValues[iQuadPoint] = evaluateBasisValue(d_basisFunctions[i],&quadraturePoint[0])*
				evaluateBasisValue(d_basisFunctions[j], &quadraturePoint[0]);

			}

			quadValues.setElementValues(iElem, elementQuadValues);

		  }

		  SMat[i][j] = dft::FieldIntegrator().integrate(quadValues,
			  meshId);
		}
	  }

	  return SMat;
	}

  std::vector<std::vector<double> > 
	GaussianFunctionManagerMO::getDensityMat()
	{

	  return d_densityMat;
	}
}
