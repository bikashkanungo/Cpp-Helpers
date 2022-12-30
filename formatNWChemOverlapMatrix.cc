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

void
readNWChemOverlap(std::vector<std::vector<double> > & SMat, const int numBasis, const int numIgnoreLines, const std::string fileName)
{

    std::ifstream readFile;
    std::string readLine;
    readFile.open(fileName.c_str());
    assert(readFile.is_open()); 
    
    SMat.resize(numBasis, std::vector<double>(0));
    int lineCount = 0;

    while(std::getline(readFile,readLine))
    {
        if(lineCount < numIgnoreLines-1)
        {
            lineCount++;
            continue;
        }

        lineCount = 0;
        for(unsigned int i = 0; i < numBasis; ++i)
        {
            std::getline(readFile,readLine);
            std::istringstream lineString(readLine);
            std::string word;
            unsigned int wordCount = 0;

            while(lineString >> word)
            {

                double val;
                if(!isNumber(val,word))
                    std::cout << "Undefined value read for the overlap file: " << fileName << std::endl;
                if(wordCount > 0)
                    SMat[i].push_back(val);
                wordCount++;
            }
        }

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
    std::cout << "Enter number of lines to ignore: ";
    std::cin >> numIgnoreLines;
    std::cout << std::endl;
    std::cout << "Enter ouput file name: ";
    std::cin >> outputFileName;
    std::cout << std::endl;
    
    std::vector<std::vector<double> > SMat(0);
    readNWChemOverlap(SMat, numBasis, numIgnoreLines, nwchemOverlapFileName);

    std::ofstream outfile;
    outfile.open(outputFileName.c_str());
    outfile << std::setprecision(16);

    for(unsigned int i = 0; i < numBasis; ++i)
    {
        assert(SMat[i].size() == numBasis);
        for(unsigned int j = 0; j < numBasis; ++j)
            outfile << std::setprecision(16) << SMat[i][j] << "\t";
        outfile <<  std::endl;
    }

    outfile.close();
}
