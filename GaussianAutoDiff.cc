#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <typeinfo>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace boost::math::constants;
using namespace boost::math::differentiation;

std::vector<double> 
getFDWeightsFirstDerivative(const int stencil) 
{

    std::vector<double> returnValue(2*stencil+1);
    if(stencil == 1)
    {
        returnValue[0] = -1.0/2.0;
        returnValue[1] = 0.0;
        returnValue[2] = 1.0/2.0;
        return returnValue;
    }

    if(stencil == 2)
    {

        returnValue[0] = 1.0/12.0;
        returnValue[1] = -2.0/3.0;
        returnValue[2] = 0.0;
        returnValue[3] = 2.0/3.0;
        returnValue[4] = -1.0/12.0;
        return returnValue;
    }

    if(stencil == 3)
    {

        returnValue[0] = -1.0/60.0;
        returnValue[1] = 3.0/20.0;
        returnValue[2] = -3.0/4.0;
        returnValue[3] = 0.0;
        returnValue[4] = 3.0/4.0;
        returnValue[5] = -3.0/20.0;
        returnValue[6] = 1.0/60.0;
        return returnValue;
    }

    if(stencil == 4)
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
        return returnValue;
    }


}
    std::vector<double> 
getFDWeightsSecondDerivative(const int stencil) 
{

    std::vector<double> returnValue(2*stencil+1);

    if(stencil == 1)
    {
        returnValue[0] = 1.0;
        returnValue[1] = -2.0;
        returnValue[2] = 1.0;
        return returnValue;
    }

    if(stencil == 2)
    {

        returnValue[0] = -1.0/12.0;
        returnValue[1] = 4.0/3.0;
        returnValue[2] = -5.0/2.0;
        returnValue[3] = 4.0/3.0;
        returnValue[4] = -1.0/12.0;
        return returnValue;
    }

    if(stencil == 3)
    {

        returnValue[0] = 1.0/90.0;
        returnValue[1] = -3.0/20.0;
        returnValue[2] = 3.0/2.0;
        returnValue[3] = -49.0/18.0;
        returnValue[4] = 3.0/2.0;
        returnValue[5] = -3.0/20.0;
        returnValue[6] = 1.0/90.0;
        return returnValue;
    }

    if(stencil == 4)
    {

        returnValue[0] = -1.0/560.0;
        returnValue[1] = 8.0/315.0;
        returnValue[2] = -1.0/5.0;
        returnValue[3] = 8.0/5.0;
        returnValue[4] = -205.0/72.0;
        returnValue[5] = 8.0/5.0;
        returnValue[6] = -1.0/5.0;
        returnValue[7] = 8.0/315.0;
        returnValue[8] = -1.0/560.0;
        return returnValue;
    }
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
    return sqrt(((2.0*l + 1)*boost::math::factorial<double>(l-m))/(2.0*boost::math::factorial<double>(l+m)));
}

template <typename T>
T
Qm(const int m, const T & phi)
{

    if(m > 0)
        return cos(m*phi);
    else if(m == 0)
        return 1.0;
    else 
        return sin(std::abs(m)*phi);
}

inline
double
dQmDPhi(const int m, const double phi)
{

    if(m > 0)
        return -m*sin(m*phi);
    else if(m == 0)
        return 0.0;
    else 
        return std::abs(m)*cos(std::abs(m)*phi);

}

inline
    double
getLimitingValue(const int l, const int m, const double theta)
{
    double returnValue = 0.0;

    if(std::fabs(theta-0.0) < 1e-12)
    {

        if(m == 0)
            returnValue = -0.5*l*(l+1);
        if(m == 2)
            returnValue = 0.25*(l-1)*l*(l+1)*(l+2);
    }

    if(std::fabs(theta-M_PI) < 1e-12)
    {

        if(m == 0)
            returnValue = -0.5*l*(l+1)*pow(-1.0,l);
        if(m == 2)
            returnValue = 0.25*(l-1)*l*(l+1)*(l+2)*pow(-1.0,l);;
    }

    return returnValue;
}

template <typename T> 
T
Plm(const int l, const int m, const T & x)
{
    if(std::abs(m) > l)
        return 0.0;
    else
        return pow(-1.0,m)*boost::math::legendre_p(l, m, x);
}

    inline
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


inline
    double
d2PlmDTheta2(const int l, const int m, const double theta)
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


template <typename T>
void
convertCartesianToSpherical(const T & x, const T & y, const T & z, T & r, T & theta, T & phi)
{

    r = sqrt(x*x + y*y + z*z);

    if(r == 0.0)
    {

        theta = 0.0;
        phi = 0.0;

    }

    else
    {

        theta = acos(z/r);
        //
        // check if theta = 0 or PI (i.e, whether the point is on the Z-axis)
        // If yes, assign phi = 0.0.
        // NOTE: In case theta = 0 or PI, phi is undetermined. The actual value 
        // of phi doesn't matter in computing the enriched function value or 
        // its gradient. We assign phi = 0.0 here just as a dummy value
        //
        if(std::fabs(theta - 0.0) >= 1e-12 && std::fabs(theta - M_PI) >= 1e-12)
            phi = atan2(y,x);

        else
            phi = 0.0;

    }

}


//inline 
//    double
//getPolarFunctionVal(const double theta, const int l, const int m)
//{
//
//    const int modM = std::abs(m);
//    double polarFunctionVal;
//
//    if(modM > l)
//        polarFunctionVal = 0.0;
//
//    else
//    {
//
//        const double cosTheta = cos(theta);
//        polarFunctionVal =  Plm(l,modM,cosTheta); 
//        polarFunctionVal *= Clm(l,modM);
//
//    }
//
//    return polarFunctionVal;
//
//}

//
// Compute the azimuthal part of the enriched function.
// Real form of spherical harmonic is used for the non-periodic problem.
// The azimuthal part is:
// cos(m*phi)/sqrt(PI) 	if  m > 0;
// sin(|m|*phi)/sqrt(PI)  	if m < 0
// 1/sqrt(2*PI) 		if m = 0;
//

//inline 
//    double
//getAzimuthalFunctionVal(const double phi, const int m)
//{
//    double azimuthalFunctionVal;
//    if(m > 0) 
//        azimuthalFunctionVal = cos(m*phi);
//
//    else if(m == 0)
//        azimuthalFunctionVal = 1.0;
//
//    else if(m < 0)
//        azimuthalFunctionVal = sin(std::abs(m)*phi);
//
//    azimuthalFunctionVal *= Dm(m);
//
//    return azimuthalFunctionVal;
//
//}

template <typename T>
T
getRadialFunctionVal(const T & r, const int l, const double alpha)
{

    return pow(r,l)*exp(-alpha*r*r);
}

double derivativeRGaussian(const double r, const double alpha, const int l, const int derOrder)
{
    if(derOrder == 0 && l >= 0)
        return pow(r,l)*exp(-alpha*r*r);
    else if(derOrder == 0 && l < 0)
        return 0.0;
    else
        return l*derivativeRGaussian(r,alpha,l-1,derOrder-1) - 2*alpha*derivativeRGaussian(r,alpha,l+1,derOrder-1);
}



template<typename X, typename Y, typename Z>
promote<X, Y, Z>
getSphericalGaussianFunctionVal(const X & x,
				const Y & y,
				const Z & z,
        const int l,
        const int m,
        const double alpha)
{

    promote<X,Y,Z> r, theta, phi;
    //convertCartesianToSpherical(x, y, z, r, theta, phi);

    r = sqrt(x*x + y*y + z*z);

    if(r == 0.0)
    {

        theta = 0.0;
        phi = 0.0;

    }

    else
    {

        theta = acos(z/r);
        //
        // check if theta = 0 or PI (i.e, whether the point is on the Z-axis)
        // If yes, assign phi = 0.0.
        // NOTE: In case theta = 0 or PI, phi is undetermined. The actual value 
        // of phi doesn't matter in computing the enriched function value or 
        // its gradient. We assign phi = 0.0 here just as a dummy value
        //
        if(theta != 0.0 && theta != M_PI)
            phi = atan2(y,x);

        else
            phi = 0.0;

    }
    const int modM = std::abs(m);
    const double C = Clm(l,modM)*Dm(m);
    const auto cosTheta = cos(theta);
    const auto R = getRadialFunctionVal(r, l , alpha);
    const auto P = Plm(l,modM,cosTheta);
    const auto Q = Qm(m,phi);
    return C*R*P*Q;

    //const double t1 = getRadialFunctionVal(r, l , alpha);
    //const double t2 = getPolarFunctionVal(theta, l, m);
    //const double t3 = getAzimuthalFunctionVal(phi,m);

    //return t1*t2*t3;
}

//template<typename T>
//T
//getSphericalGaussianFunctionVal(const T & x,
//				const T & y,
//				const T & z,
//        const int l,
//        const int m,
//        const double alpha)
//{
//
//    T r, theta, phi;
//    //convertCartesianToSpherical(x, y, z, r, theta, phi);
//
//    r = sqrt(x*x + y*y + z*z);
//
//    if(r == 0.0)
//    {
//
//        theta = 0.0;
//        phi = 0.0;
//
//    }
//
//    else
//    {
//
//        theta = acos(z/r);
//        //
//        // check if theta = 0 or PI (i.e, whether the point is on the Z-axis)
//        // If yes, assign phi = 0.0.
//        // NOTE: In case theta = 0 or PI, phi is undetermined. The actual value 
//        // of phi doesn't matter in computing the enriched function value or 
//        // its gradient. We assign phi = 0.0 here just as a dummy value
//        //
//        if(theta != 0.0 && theta != M_PI)
//            phi = atan2(y,x);
//
//        else
//            phi = 0.0;
//
//    }
//    const int modM = std::abs(m);
//    const double C = Clm(l,modM)*Dm(m);
//    const auto cosTheta = cos(theta);
//    const auto R = getRadialFunctionVal(r, l , alpha);
//    const auto P = Plm(l,modM,cosTheta);
//    const auto Q = Qm(m,phi);
//    return C*R*P*Q;
//
//    //const double t1 = getRadialFunctionVal(r, l , alpha);
//    //const double t2 = getPolarFunctionVal(theta, l, m);
//    //const double t3 = getAzimuthalFunctionVal(phi,m);
//
//    //return t1*t2*t3;
//}

//inline 
//    std::vector<double>
//getSphericalGaussianGradientFD(const std::vector<double> & x,
//        const int l,
//        const int m,
//        const double alpha,
//        const double delta,
//        const int stencil)
//{
//
//    std::vector<double> returnValue(3,0.0);
//    double f = getSphericalGaussianFunctionVal(x, l, m, alpha);
//    std::vector<double> FDWeights = getFDWeightsFirstDerivative(stencil);
//    for(int i = 0; i < 3; ++i)
//    {
//        for(int j = 0; j < 2*stencil+1; ++j)
//        {
//
//            std::vector<double> xCopy(x);
//            xCopy[i] += delta*(j-stencil);
//            double fVal = getSphericalGaussianFunctionVal(xCopy, l, m, alpha);
//            returnValue[i] += FDWeights[j]*fVal;
//        }
//
//        returnValue[i] /= delta;
//    }
//
//    return returnValue;
//}

inline
std::vector<double>
getSphericalGaussianGradient(const std::vector<double> & x,
        const int l,
        const int m,
        const double alpha)
{
    double r, theta, phi;
    convertCartesianToSpherical(x[0], x[1], x[2], r, theta, phi);
    std::vector<double> returnValue(3);
    const int modM = std::abs(m);
    const double C = Clm(l,modM)*Dm(m);
    if(r < 1e-15)
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

    else if(std::fabs(theta-0.0) < 1e-12)
    {
        const double R = getRadialFunctionVal(r,l,alpha);
        const double dRDr = derivativeRGaussian(r, alpha, l, 1);
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

    else if(std::fabs(theta-M_PI) < 1e-12)
    {
        const double R = getRadialFunctionVal(r,l,alpha);
        const double dRDr = derivativeRGaussian(r, alpha, l, 1);
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

        const double R = getRadialFunctionVal(r,l,alpha);
        const double dRDr = derivativeRGaussian(r, alpha, l, 1);
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


//inline 
//    double
//getSphericalGaussianLaplacianFD(const std::vector<double> & x,
//        const int l,
//        const int m,
//        const double alpha,
//        const double delta,
//        const int stencil)
//{
//
//    double returnValue = 0.0;
//    double f = getSphericalGaussianFunctionVal(x, l, m, alpha);
//    //double delta = 1e-4;
//    std::vector<double> FDWeights = getFDWeightsSecondDerivative(stencil);
//
//    //for(unsigned int j = 0; j < 2*stencil + 1; ++j)
//    //    std::cout << "FDWeights["<< j << "]: " << FDWeights[j] << std::endl;
//    for(int i = 0; i < 3; ++i)
//    {
//
//        double der = 0.0;
//        for(int j = 0; j < 2*stencil+1; ++j)
//        {
//
//            std::vector<double> xCopy(x);
//            xCopy[i] += delta*(j-stencil);
//            double fVal = getSphericalGaussianFunctionVal(xCopy, l, m, alpha);
//            der += FDWeights[j]*fVal;
//        }
//
//        returnValue += der/(delta*delta);
//
//        //std::vector<double> xPrev(x);
//        //std::vector<double> xNext(x);
//        //xPrev[i] -= delta;
//        //xNext[i] += delta;
//        //double fPrev = getSphericalGaussianFunctionVal(xPrev, l, m, alpha);
//        //double fNext = getSphericalGaussianFunctionVal(xNext, l, m, alpha);
//
//        //returnValue += (fPrev - 2*f + fNext)/(delta*delta);
//    }
//
//    return returnValue;
//}

inline 
    double
getSphericalGaussianLaplacian(const std::vector<double> & x,
        const int l,
        const int m,
        const double alpha)
{

    double r, theta, phi;
    convertCartesianToSpherical(x[0], x[1], x[2], r, theta, phi);
    double returnValue = 0.0;
    if(r < 1e-15)
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
        const double R = getRadialFunctionVal(r, l , alpha);
        const double dRdr = derivativeRGaussian(r, alpha, l, 1);
        const double d2Rdr2 = derivativeRGaussian(r, alpha, l, 2);
        const double P = Plm(l,modM,cosTheta);
        const double Q = Qm(m,phi);
        const double term1 = C*P*Q*(2.0*dRdr/r + d2Rdr2);

        if(std::fabs(theta-0.0) < 1e-12 || std::fabs(theta-M_PI) < 1e-12)
        {
            const double limitingVal = getLimitingValue(l, modM, theta);
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
            //std::cout << std::setprecision(18) << std::endl;
            //std::cout << "t1: " << term1 << " t2: " << term2 << " t3: " << term3 << std::endl;

        }

    }

    return returnValue;
}

//    void
//testDPlmDTheta(const int l, const int m, const double theta, const double delta)
//{
//
//    const double cosTheta = cos(theta);
//    const double f = Plm(l, m, cosTheta);
//
//    const double cosThetaPrev = cos(theta-delta);
//    const double fPrev = Plm(l,m,cosThetaPrev);
//
//    const double cosThetaNext = cos(theta+delta);
//    const double fNext = Plm(l,m,cosThetaNext);
//
//    const double FDVal = (fNext-fPrev)/(2.0*delta);
//    const double exactVal = dPlmDTheta(l,m,theta);
//
//    std::cout << std::setprecision(18) << std::endl;
//    std::cout << "DPlmDTheta FDVal: " << FDVal << std::endl;
//    std::cout << "DPlmDTheta ExactVal: " << exactVal << std::endl;
//
//}
//
//    void
//testD2PlmDTheta2(const int l, const int m, const double theta, const double delta)
//{
//
//    const double cosTheta = cos(theta);
//    const double f = Plm(l, m, cosTheta);
//
//    const double cosThetaPrev = cos(theta-delta);
//    const double fPrev = Plm(l,m,cosThetaPrev);
//
//    const double cosThetaNext = cos(theta+delta);
//    const double fNext = Plm(l,m,cosThetaNext);
//
//    const double FDVal = (fPrev -2.0*f + fNext)/(delta*delta);
//    const double exactVal = d2PlmDTheta2(l,m,theta);
//
//    std::cout << std::setprecision(18) << std::endl;
//    std::cout << "D2PlmDTheta2 FDVal: " << FDVal << std::endl;
//    std::cout << "D2PlmDTheta2 ExactVal: " << exactVal << std::endl;
//
//}

std::vector<double>
getSphericalGaussianGradientAutoDiff(const std::vector<double> & point,
																		 const int l,
																		 const int m,
																		 const double alpha)
{

		constexpr unsigned orderX = 1;
		constexpr unsigned orderY = 1;
		constexpr unsigned orderZ = 1;
		auto const variables = make_ftuple<double, orderX, orderY, orderZ>(point[0], point[1], point[2]);
		auto const& x = std::get<0>(variables);  
		auto const& y = std::get<1>(variables);  
		auto const& z = std::get<2>(variables);  
		auto const v = getSphericalGaussianFunctionVal(x, y, z, l, m, alpha);
		std::vector<double> returnValue(3,0.0);
		returnValue[0] =  v.derivative(1,0,0);
		returnValue[1]=  v.derivative(0,1,0);
		returnValue[2] =  v.derivative(0,0,1);
		return returnValue;
}

double
getSphericalGaussianLaplacianAutoDiff(const std::vector<double> & point,
																		 const int l,
																		 const int m,
																		 const double alpha)
{

		constexpr unsigned orderX = 2;
		constexpr unsigned orderY = 2;
		constexpr unsigned orderZ = 2;
		auto const variables = make_ftuple<double, orderX, orderY, orderZ>(point[0], point[1], point[2]);
		auto const& x = std::get<0>(variables);  
		auto const& y = std::get<1>(variables);  
		auto const& z = std::get<2>(variables);  
		auto const v = getSphericalGaussianFunctionVal(x, y, z, l, m, alpha);
		double returnValue = 0.0;
		returnValue +=  v.derivative(2,0,0);
		returnValue +=  v.derivative(0,2,0);
		returnValue +=  v.derivative(0,0,2);
		return returnValue;
}

void
testGradient(const int alphaPowMin, 
             const int alphaPowMax, 
             const int scalePowMin, 
             const int scalePowMax)
{

    double L1Error = 0.0;
    double L1RelError = 0.0;
    int numTrials = 100;
    int totalNumTrials = numTrials*(alphaPowMax-alphaPowMin+1)*(scalePowMax-scalePowMin+1);

    double diffMax = 0.0, alphaAtMaxDiff, FDValAtMaxDiff, exactValAtMaxDiff;
    std::vector<double> xAtMaxDiff;
    int lAtMaxDiff, mAtMaxDiff;
    std::cout << std::setprecision(18) << std::endl;
    std::vector<double> x(3);
    for(int alphaPow = alphaPowMin; alphaPow <= alphaPowMax; ++alphaPow)
    {
        double alpha = pow(10.0,alphaPow);
        for(int scalePow = scalePowMin; scalePow <= scalePowMax; ++scalePow)
        {
            double scale = pow(10.0,scalePow);
            for(unsigned int i = 0; i < numTrials; ++i)
            {

                x[0] = scale*(rand()/(RAND_MAX+0.0)-0.5);
                x[1] = scale*(rand()/(RAND_MAX+0.0)-0.5);
                x[2] = 0.0;//scale*(rand()/(RAND_MAX+0.0)-0.5);

                for(int l = 0; l <=5; ++l)
                {
                    for(int m = -l; m <= l; ++m)
                    {

                        const std::vector<double> l1 = getSphericalGaussianGradientAutoDiff(x, l, m, alpha);
                        const std::vector<double> l2 = getSphericalGaussianGradient(x, l, m, alpha);
                        double diff = 0.0;
                        double l1Abs = 0.0;
                        double l2Abs = 0.0;
                        for(unsigned int i = 0; i < 3; ++i)
                        {
                            diff += pow(l1[i]-l2[i],2.0);
                            l1Abs += l1[i]*l1[i];
                            l2Abs += l2[i]*l2[i];
                        }

                        diff = sqrt(diff);
                        l2Abs = sqrt(l2Abs);

                        double diffRel;
                        if(std::fabs(l2Abs) > 1e-12)
                            diffRel = diff/std::fabs(l2Abs);
                        else
                            diffRel = 0.0;

                        if(diffRel > diffMax)
                        {
                            diffMax = diffRel;
                            alphaAtMaxDiff = alpha;
                            FDValAtMaxDiff = l1Abs;
                            exactValAtMaxDiff = l2Abs;
                            xAtMaxDiff = x;
                            lAtMaxDiff = l;
                            mAtMaxDiff = m;
                        }


                        L1Error += diff;
                        L1RelError += diffRel;
                        if(diffRel > 1e-6)
                            std::cout << "Diff: " << diffRel << " l: " << l << " m: " << m << " x: " << x[0] << " y: " << x[1] << " z: " << x[2] << std::endl;
                        //std::cout << "Clm: " << Clm(l,std::abs(m)) << std::endl;
                        //std::cout <<"Exact: " << l2 <<  " FD: " << l1 << " Diff: " << std::fabs(l1-l2) << std::endl;
                        //std::cout << "Exact: " << l2 << std::endl;
                        //std::cout << "Difference: " << std::fabs(l1-l2) << std::endl;

                    }
                }
            }
        }
    }
    std::cout << "L1Error: " << L1Error/totalNumTrials << std::endl;
    std::cout << "L1RelError: " << L1RelError/totalNumTrials << std::endl;
    std::cout << "DiffMax: " << diffMax << std::endl;
    std::cout << "FDValAtMaxDiff: " << FDValAtMaxDiff << std::endl;
    std::cout << "exactValAtMaxDiff: " << exactValAtMaxDiff << std::endl;
    std::cout << "lAtMaxDiff:  "<< lAtMaxDiff << std::endl;
    std::cout << "mAtMaxDiff:  "<< mAtMaxDiff << std::endl;
    std::cout << "xAtMaxDiff:  " << xAtMaxDiff[0] << " " << xAtMaxDiff[1] << " " << xAtMaxDiff[2] << std::endl;
    std::cout << "alphaAtMaxDiff: " << alphaAtMaxDiff << std::endl;

}


void
testLaplacian(const int alphaPowMin, 
              const int alphaPowMax, 
              const int scalePowMin, 
              const int scalePowMax) 
{

    double L1Error = 0.0;
    double L1RelError = 0.0;
    int numTrials = 50;
    int totalNumTrials = numTrials*(alphaPowMax-alphaPowMin+1)*(scalePowMax-scalePowMin+1);

    double diffMax = 0.0, alphaAtMaxDiff, FDValAtMaxDiff, exactValAtMaxDiff;
    std::vector<double> xAtMaxDiff;
    int lAtMaxDiff, mAtMaxDiff;
    std::cout << std::setprecision(18) << std::endl;
    std::vector<double> x(3);
    for(int alphaPow = alphaPowMin; alphaPow <= alphaPowMax; ++alphaPow)
    {
        double alpha = pow(10.0,alphaPow);
        for(int scalePow = scalePowMin; scalePow <= scalePowMax; ++scalePow)
        {
            double scale = pow(10.0,scalePow);
            for(unsigned int i = 0; i < numTrials; ++i)
            {

                x[0] = scale*(rand()/(RAND_MAX+0.0)-0.5);
                x[1] = scale*(rand()/(RAND_MAX+0.0)-0.5);
                x[2] = scale*(rand()/(RAND_MAX+0.0)-0.5);

                for(int l = 0; l <=5; ++l)
                {
                    for(int m = -l; m <= l; ++m)
                    {

                        const double l1 = getSphericalGaussianLaplacianAutoDiff(x, l, m, alpha);
                        const double l2 = getSphericalGaussianLaplacian(x, l, m, alpha);
                        const double diff = std::fabs(l1-l2);
                        double diffRel;
                        if(std::fabs(l2) > 1e-12)
                            diffRel = diff/std::fabs(l2);
                        else
                            diffRel = 0.0;

                        if(diffRel > diffMax)
                        {
                            diffMax = diffRel;
                            alphaAtMaxDiff = alpha;
                            FDValAtMaxDiff = l1;
                            exactValAtMaxDiff = l2;
                            xAtMaxDiff = x;
                            lAtMaxDiff = l;
                            mAtMaxDiff = m;
                        }


                        L1Error += diff;
                        L1RelError += diffRel;
                        if(diffRel > 1e-6)
                            std::cout << "Diff: " << diffRel << " l: " << l << " m: " << m << " x: " << x[0] << " y: " << x[1] << " z: " << x[2] << std::endl;
                        //std::cout << "Clm: " << Clm(l,std::abs(m)) << std::endl;
                        //std::cout <<"Exact: " << l2 <<  " FD: " << l1 << " Diff: " << std::fabs(l1-l2) << std::endl;
                        //std::cout << "Exact: " << l2 << std::endl;
                        //std::cout << "Difference: " << std::fabs(l1-l2) << std::endl;

                    }
                }
            }
        }
    }
    std::cout << "L1Error: " << L1Error/totalNumTrials << std::endl;
    std::cout << "L1RelError: " << L1RelError/totalNumTrials << std::endl;
    std::cout << "DiffMax: " << diffMax << std::endl;
    std::cout << "FDValAtMaxDiff: " << FDValAtMaxDiff << std::endl;
    std::cout << "exactValAtMaxDiff: " << exactValAtMaxDiff << std::endl;
    std::cout << "lAtMaxDiff:  "<< lAtMaxDiff << std::endl;
    std::cout << "mAtMaxDiff:  "<< mAtMaxDiff << std::endl;
    std::cout << "xAtMaxDiff:  " << xAtMaxDiff[0] << " " << xAtMaxDiff[1] << " " << xAtMaxDiff[2] << std::endl;
    std::cout << "alphaAtMaxDiff: " << alphaAtMaxDiff << std::endl;

}

int main()
{

    //int l,m;
    //double alpha;
    //std::vector<double> point(3);
		//std::cout << std::endl;
    //std::cout << "Enter l: ";
    //std::cin >> l;
		//std::cout << std::endl;
    //std::cout << "Enter m: ";
    //std::cin >> m;
    //std::cout << std::endl;
    //std::cout << "Enter alpha: ";
    //std::cin >> alpha;
    //std::cout << std::endl;
    //std::cout << "Enter x: ";
    //std::cin >> point[0];
    //std::cout << std::endl;
    //std::cout << "Enter y: ";
    //std::cin >> point[1];
    //std::cout << std::endl;
    //std::cout << "Enter z: ";
    //std::cin >> point[2];
    //std::cout << std::endl;

		//std::cout << "l: " << l << " m: " << m << std::endl;
		//std::cout << "point: (" << point[0] << ", " << point[1] << ", " << point[2] << ")" << std::endl;
    int alphaPowMin, alphaPowMax;
    int scalePowMin, scalePowMax;
    std::cout << "Enter alphaPowMin: ";
    std::cin >> alphaPowMin;
    std::cout << std::endl;
    std::cout << "Enter alphaPowMax: ";
    std::cin >> alphaPowMax;
    std::cout << std::endl;
    std::cout << "Enter scalePowMin: ";
    std::cin >> scalePowMin;
    std::cout << std::endl;
    std::cout << "Enter scalePowMax: ";
    std::cin >> scalePowMax;
    std::cout << std::endl;

		std::cout << "\n\nTesting Gradient" << std::endl;
		std::cout << "-------------------------------------------------" << std::endl;
    testGradient(alphaPowMin, 
                 alphaPowMax, 
                 scalePowMin, 
                 scalePowMax);

		std::cout << "\n\nTesting Laplacian" << std::endl;
		std::cout << "-------------------------------------------------" << std::endl;
		testLaplacian(alphaPowMin, 
								  alphaPowMax, 
								  scalePowMin, 
								  scalePowMax);

			//constexpr unsigned orderX = 1;
			//constexpr unsigned orderY = 1;
			//constexpr unsigned orderZ = 1;
			////auto const x = make_fvar<double, orderX>(point[0]);				
			////auto const y = make_fvar<double, orderY>(point[1]);				
			////auto const z = make_fvar<double, orderZ>(point[2]);				
			//auto const variables = make_ftuple<double, orderX, orderY, orderZ>(point[0], point[1], point[2]);
			//auto const& x = std::get<0>(variables);  
			//auto const& y = std::get<1>(variables);  
			//auto const& z = std::get<2>(variables);  
			////std::cout << typeid(x).name() << std::endl;
			////std::cout << typeid(y).name() << std::endl;
			////std::cout << decltype(y) << std::endl;
			////std::cout << decltype(z) << std::endl;
			//auto const v = getSphericalGaussianFunctionVal(x, y, z, l, m, alpha);
			//const double fx =  v.derivative(1,0,0);
			//const double fy =  v.derivative(0,1,0);
			//const double fz =  v.derivative(0,0,1);
			//std::vector<double> gradient = getSphericalGaussianGradient(point, l , m, alpha);
			//std::cout << "fx: (" << fx << ", " << gradient[0] << ")" << std::endl;
			//std::cout << "fy: (" << fy << ", " << gradient[1] << ")" << std::endl;
			//std::cout << "fz: (" << fz << ", " << gradient[2] << ")" << std::endl;

}
