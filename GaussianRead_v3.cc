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

#define LENGTH 10.0
#define ANGSTROM_TO_BOHR 1.889725989

typedef double (*gaussianFunction)(const double , const double *, const double *);
typedef double (*func)(const double *, void *);

struct contractedGaussian {
    int L;
    double * alpha;
    double * c;
    char lquantum;
};

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

void
readGaussianFiles(std::vector<contractedGaussian*> & atomicContractedGaussians, const std::string fileName)
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
            std::size_t pos = word.find_first_not_of("SPDF");
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
                                std::cout << "Undefined value read for the exponent in file: " << fileName << std::endl;
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
                        contractedGaussian * cg = new contractedGaussian;
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

void readDensityMatrix(std::vector<std::vector<double> > & densityMat, const std::string fileName)
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
int main()
{
    double C [] = {0.0000000000,0.0000000000,0.0000000000};
    double H1 [] = {0.0000000000,0.0000000000,1.1089000000};
    double H2 [] = {1.0851085057,0.0000000000,-0.2284704375};

    std::vector<double *> atomicCoords;
    atomicCoords.push_back(&C[0]);
    atomicCoords.push_back(&H1[0]);
    atomicCoords.push_back(&H2[0]);

    std::vector<std::string> basisFileNames;
    basisFileNames.push_back("C_gaussian");
    basisFileNames.push_back("H_gaussian");
    basisFileNames.push_back("H_gaussian");
  
    std::string densityMatrixFile("DensityMatrix");
    std::string SMatrixFile("SMatrix");

    std::set<std::string> uniqueBasisFileNames;
    for(unsigned int i = 0; i < basisFileNames.size(); ++i)
        uniqueBasisFileNames.insert(basisFileNames[i]);

    const int numUniqueBasisFiles = uniqueBasisFileNames.size();
    std::vector<std::vector<contractedGaussian*> > contractedGaussians(numUniqueBasisFiles);
    unsigned int i = 0;
    for(std::set<std::string>::const_iterator iter = uniqueBasisFileNames.begin();
            iter!= uniqueBasisFileNames.end(); iter++)
    {
        std::vector<contractedGaussian*> & atomicContractedGaussians = 
            contractedGaussians[i];
        std::string fileName = *iter;
        readGaussianFiles(atomicContractedGaussians, fileName);
        
        i++;
    }

    for(unsigned int i = 0; i < numUniqueBasisFiles; ++i)
    {
        std::vector<contractedGaussian*> & atomicContractedGaussians = 
            contractedGaussians[i];
        const int numContractedGaussians = atomicContractedGaussians.size();
        for(unsigned int j = 0; j < numContractedGaussians; ++j)
        {
            const contractedGaussian * cg = atomicContractedGaussians[j];
            const int L = cg->L;
            const char lquantum = cg->lquantum;
            std::cout << "NumContracted: " <<  L << " l character: " << lquantum << std::endl;
            for(unsigned int k = 0; k < L; ++k)
            {
                std::cout << cg->alpha[k] << "\t" << cg->c[k] << std::endl;
                cg->alpha[k] *= 1.0/pow(ANGSTROM_TO_BOHR,2.0);
            }
        }
    }

    std::vector<std::vector<double> > densityMat(0);
    readDensityMatrix(densityMat, densityMatrixFile);
    for(unsigned int i = 0; i < densityMat.size(); ++i)
    {
        for(unsigned int j = 0; j < densityMat[i].size(); ++j)
            std::cout << densityMat[i][j] << "\t";
        std::cout << std::endl;
    }
    
    std::vector<std::vector<double> > SMat(0);
    readDensityMatrix(SMat, SMatrixFile);
    for(unsigned int i = 0; i < SMat.size(); ++i)
    {
        for(unsigned int j = 0; j < SMat[i].size(); ++j)
            std::cout << SMat[i][j] << "\t";
        std::cout << std::endl;
    }

    std::vector<basis*> basisFunctions(0);
    for(unsigned int iAtom = 0; iAtom < atomicCoords.size(); ++iAtom)
    {
        for(unsigned int j = 0; j < 3; ++j)
            atomicCoords[iAtom][j] *= ANGSTROM_TO_BOHR;
        const std::string basisFileName = basisFileNames[iAtom];
        std::set<std::string>::const_iterator iter = uniqueBasisFileNames.find(basisFileName);
        unsigned int index = std::distance(uniqueBasisFileNames.begin(), iter);
        std::vector<contractedGaussian*> & atomicContractedGaussians = 
            contractedGaussians[index];
        unsigned int numContractedGaussians = atomicContractedGaussians.size();
        for(unsigned int j = 0; j < numContractedGaussians; ++j)
        {
            const contractedGaussian * cg = atomicContractedGaussians[j];
            if(cg->lquantum == 'S')
            {
                basis * b = new basis;
                b->L = cg->L;
                b->alpha = cg->alpha;
                b->c = cg->c;
                b->origin = atomicCoords[iAtom];
                b->gf = &g1s;
                b->normConst = 1.0;
                basisFunctions.push_back(b);
            }
            else if(cg->lquantum == 'P')
            {
                for(unsigned int cart = 0; cart < 3; ++cart)
                {
                    basis * b = new basis;
                    b->L = cg->L;
                    b->alpha = cg->alpha;
                    b->c = cg->c;
                    b->origin = atomicCoords[iAtom];
                    if(cart == 0)
                        b->gf = &g2px;
                    else if(cart == 1)
                        b->gf = &g2py;
                    else
                        b->gf = &g2pz;
                    b->normConst = 1.0;
                    basisFunctions.push_back(b);
                }
            }

            else
                std::cout << "Invalid L character detected" << std::endl;
        }
          
    }

    unsigned int numBasis = basisFunctions.size();
    assert(numBasis == densityMat.size());
    std::cout << "Number basis: " << numBasis << std::endl;
    
    double h;
    double len;
    std::cout <<"Enter length:";
    std::cin >> len;
    std::cout << std::endl;
    std::cout <<"Enter h:";
    std::cin >> h;
    std::cout << std::endl;
    double N = len/h;
    for(unsigned int i = 0; i < numBasis; ++i)
    {

        std::vector<const basis *> funcs(2,basisFunctions[i]);
        double integral = trapezoidal3d(len,h, funcs);
        std::cout << "Norm-const for basis: " << i << " =  " << integral << std::endl;
        basisFunctions[i]->normConst = 1.0/sqrt(integral);
    }

    double rhoIntegral = 0.0;
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
                double val = 0.0;
                for(unsigned int l = 0; l < numBasis; ++l)
                {
                    for(unsigned int m = 0; m < numBasis; ++m)
                    {
                        val += densityMat[l][m]*evaluateBasisValue(basisFunctions[l], &point[0])*
                            evaluateBasisValue(basisFunctions[m],&point[0]);
                    }
                }
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
            rhoIntegral += 0.5*h*f2;
        else
            rhoIntegral += h*f2;

    }

    double rhoIntegralTest = 0.0;
    for(unsigned int i = 0; i < numBasis; ++i)
        for(unsigned int j = 0; j < numBasis; ++j)
            rhoIntegralTest += densityMat[i][j]*SMat[i][j];

    std::cout << "Rho integral exact: " << 2.0*rhoIntegralTest << std::endl;
    std::cout << "Rho integral: " << 2.0*rhoIntegral << std::endl;
    //double h;
    //double len;
    //std::cout <<"Enter length:";
    //std::cin >> len;
    //std::cout << std::endl;
    //std::cout <<"Enter h:";
    //std::cin >> h;
    //std::cout << std::endl;
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
    //basis b1;
    //b1.L = 6;
    //b1.alpha = &alpha_C6s[0];
    //b1.c = &c_C6s[0];
    //b1.origin = &C[0];
    //b1.gf = &g1s;
    //b1.normConst = 1.0;

    //// C_3p
    //basis b3;
    //b3.L = 3;
    //b3.alpha = &alpha_C3p[0];
    //b3.c = &c_C3p[0];
    //b3.origin = &C[0];
    //b3.gf = &g2px;
    //b3.normConst = 1.0;


    //// C_1s
    //basis b6;
    //b6.L = 1;
    //b6.alpha = &alpha_C1s[0];
    //b6.c = &c_C1s[0];
    //b6.origin = &C[0];
    //b6.gf = &g1s;
    //b6.normConst = 1.0;

    //// C_1p
    //basis b7;
    //b7.L = 1;
    //b7.alpha = &alpha_C1p[0];
    //b7.c = &c_C1p[0];
    //b7.origin = &C[0];
    //b7.gf = &g2px;
    //b7.normConst = 1.0;
    //
    //std::vector<const basis *> funcs(2,&b3);
    //double integral = trapezoidal3d(len,h, funcs);
    //b1.normConst = 1.0/sqrt(integral);
    //std::cout << "Integral Gaussian: " << integral << std::endl;

    //funcs[1] = &b7;
    //integral = trapezoidal3d(len,h,funcs);
    //std::cout << "Integral 3 & 7: " << integral << std::endl;


}

