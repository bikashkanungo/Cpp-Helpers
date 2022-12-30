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

struct atom {
    std::string symbol;
    double x[3];
};

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
readAtomicInfo(const std::string fileName,
               std::vector<atom> & atomicInfo)
{
    
    std::ifstream readFile;
    
    //
    // string to read line
    //
    std::string readLine;

    readFile.open(fileName.c_str());
    assert(readFile.is_open()); 

    atom atomCurrent;
    while(std::getline(readFile,readLine))
    {
        std::istringstream lineString(readLine);
        std::string word;
        unsigned int count = 0;
        while(lineString >> word)
        {
            if(count == 0)
                atomCurrent.symbol = word;
            else if(count > 0 && count < 4)
            {
                if(!isNumber(atomCurrent.x[count-1], word))
                        std::cout << "xyz coordinate not a number" << std::endl;
            }

            else
                std::cout << "Invalid column in coordinate file" << std::endl;
            
            count++;
        }
        atomicInfo.push_back(atomCurrent);
    }
    
    readFile.close();
}

int main()
{
    int numBondIncs;
    double bondInc;
    std::string coordsFile;
    std::cout << "Enter H2O coordinate file: ";
    std :: cin >> coordsFile;
    std::cout << std::endl;
    std::cout << "Enter number of bond increments: ";
    std :: cin >> numBondIncs;
    std::cout << std::endl;
    std::cout << "Enter bond increment fraction: ";
    std :: cin >> bondInc;
    std::cout << std::endl;
    
    std::vector<atom> atomicInfo(0);
    readAtomicInfo(coordsFile,atomicInfo);

    for(unsigned int i = 0; i < atomicInfo.size(); ++i)
        std::cout << atomicInfo[i].symbol << "\t" << atomicInfo[i].x[0] << "\t" << atomicInfo[i].x[1] << "\t" <<  atomicInfo[i].x[2] << std::endl;


    const double * center = &atomicInfo[0].x[0];

    for(unsigned int i = 0; i < numBondIncs; ++i)
    {
        double bondLengthScale = 1.0 + bondInc*(i);
        std::ostringstream oss;
        oss << bondLengthScale-1.0;
        std::string fileName = coordsFile + ".bondInc_" + oss.str();
        std::ofstream file;
        file.open(fileName.c_str());

        atom atomCurrent;
        const int numAtoms = atomicInfo.size();
        for(unsigned int j = 0; j < numAtoms; ++j)
        {
            atomCurrent = atomicInfo[j];
            double bondLength = 0.0;
            for(unsigned int k = 0; k < 3; ++k)
            {
                atomCurrent.x[k] = center[k] + bondLengthScale*(atomCurrent.x[k] - center[k]);
                bondLength += pow(atomCurrent.x[k]-center[k],2.0);
            }
            
            file << atomCurrent.symbol << "\t" << atomCurrent.x[0] << "\t" << atomCurrent.x[1] << "\t" <<  atomCurrent.x[2] << std::endl;
            
            bondLength = sqrt(bondLength);
            if(j == 2)
                std::cout << "[" << i << "]" << "O-H length: " << bondLength << std::endl;
        }

        file.close();
    }
    return 0;
}
