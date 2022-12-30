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
#include <iomanip>

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
void
readNWChemOverlap(std::vector<std::vector<double> > & SMat, const int numBasis, const std::string fileName)
{

    std::ifstream readFile;
    std::string readLine;
    readFile.open(fileName.c_str());
    assert(readFile.is_open()); 
    
    SMat.resize(numBasis, std::vector<double>(numBasis,0.0));
    int lineCount = 0;

    while(std::getline(readFile,readLine))
    {
        std::istringstream lineString(readLine);
        std::string word;
        unsigned int wordCount = 0;
        int iIndex, jIndex;
        double overlapVal;
        while(lineString >> word)
        {
            if(wordCount == 1)
            {

                int val;
                if(!isInteger(val,word))
                    std::cout << "Undefined value read for the first index in the overlap file: " << fileName << std::endl;
                iIndex = val-1;
            }
            
            if(wordCount == 4)
            {

                int val;
                if(!isInteger(val,word))
                    std::cout << "Undefined value read for the second index in the overlap file: " << fileName << std::endl;
                jIndex = val-1;
            }

            if(wordCount == 7)
            {
                double val;
                if(!isNumber(val,word))
                    std::cout << "Undefined value read for the overlap value in file: " << fileName << std::endl;
                overlapVal = val;
            }
                wordCount++;
        }

        SMat[iIndex][jIndex] = overlapVal;

    }

    readFile.close();
}

int main()
{
    int numBasis, numIgnoreLines;
    std::string nwchemOverlapFileName;
    std::string outputFileName;
    std::cout << "Enter NWChem overlap matrix file name: ";
    std::cin >> nwchemOverlapFileName;
    std::cout << std::endl;
    std::cout << "Enter number of basis: ";
    std::cin >> numBasis;
    std::cout << std::endl;
    std::cout << "Enter ouput file name: ";
    std::cin >> outputFileName;
    std::cout << std::endl;
    
    std::vector<std::vector<double> > SMat(0);
    readNWChemOverlap(SMat, numBasis, nwchemOverlapFileName);

    std::ofstream outfile;
    outfile.open(outputFileName.c_str());
    outfile << std::setprecision(16);

    for(unsigned int i = 0; i < numBasis; ++i)
    {
        assert(SMat[i].size() == numBasis);
        for(unsigned int j = 0; j < numBasis; ++j)
        {   if(SMat[i][j] != SMat[j][i])
                std::cout << "Asymmetric values observed for pair " << "(" << i << "," << j << ")\t" << SMat[i][j] << "\t" << SMat[j][i] << std::endl;
            outfile << std::setprecision(16) << SMat[i][j] << "\t";
        }
        outfile <<  std::endl;
    }

    outfile.close();
}
