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

void
readSOGGAParams(std::vector<std::vector<double> > & paramsVec,
        std::string fileName)
{

    std::ifstream readFile;

    //
    // string to read line
    //
    std::string readLine;

    readFile.open(fileName.c_str());
    assert(readFile.is_open());
    int paramSetId = 0; 
    while(std::getline(readFile,readLine))
    {
        std::istringstream lineString(readLine);
        std::string word;
        lineString >> word;
        std::string wordTrimmed = trim(word);
        int numParamsInSet; 
        if(!isInteger(numParamsInSet, wordTrimmed))
        {
            std::string message("Non-integer value found for number of parameters in file: ");
            message += " " + fileName;
            std::cout << message << " Word: " << wordTrimmed << std::endl;

        }
        else
        {

            std::vector<double> & paramCurrentSet = paramsVec[paramSetId];
            paramCurrentSet.resize(numParamsInSet);
        
            for(unsigned int iParam = 0; iParam < numParamsInSet; ++iParam)
            {
                std::getline(readFile,readLine);
                std::istringstream lineString1(readLine);
                lineString1 >> word;
                wordTrimmed = trim(word);
                double val;
                if(!isNumber(val, wordTrimmed))
                {
                    std::string message("Invalid value found for number of parameters in file: ");
                    message += " " + fileName;
                    std::cout << message << std::endl;
                }
                paramCurrentSet[iParam] = val;
            }

            paramSetId++;       
        }
        std::getline(readFile,readLine); // ignore line at the end of the param set 
    }
    readFile.close();
}

int main()
{

    std::vector<std::vector<double> > paramsVec(4, std::vector<double>(0));
    std::string fileName("SOGGAParamsInput");
    readSOGGAParams(paramsVec, fileName);   
    for(unsigned int i = 0; i < 4; ++i)
    {
        int numParamsInSet = paramsVec[i].size();
        std::cout << "\nPrinting paramset: " << i << std::endl;
        std::cout << "------------------------------" << std::endl;

        for(unsigned int j = 0; j < numParamsInSet; ++j)
        {
            std::cout << paramsVec[i][j] << std::endl;
        }
    } 

}
