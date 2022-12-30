#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
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

inline
void
    readFile(std::vector<std::vector<double> > &data,
                           std::string fileName,
                           unsigned int numColumns)
    {
    
      std::vector<double> rowData(numColumns, 0.0);
      std::ifstream readFile;
      readFile.open(fileName.c_str());

      //
      // String to store line and word
      //
      std::string readLine;
      std::string word;

      //
      // column index
      //
      int columnCount;

      if(readFile.is_open())
      {
        while (std::getline(readFile, readLine))
        {
          std::istringstream iss(readLine);
        
          columnCount = 0; 
					double val;
          while(iss >> word && columnCount < numColumns)
					{

            if(isNumber(val,word))
							rowData[columnCount++] = val;
						else
              std::cout << "Undefined value read in file: " << fileName << std::endl;
					}
          data.push_back(rowData);
        }
      }

      readFile.close();
      return;
    }

int main()
{
    int numMOs;
    std::string MOFile, DensityMatFile;
    std::cout << "Enter MO coeff(s) file:";
    std::cin >> MOFile;
    std::cout << std::endl;
    std::cout <<"Enter number of MOs:";
    std::cin >> numMOs;
    std::cout << std::endl;
    std::cout << "Enter Density Matrix output filename:";
    std::cin >> DensityMatFile;
    std::cout << std::endl;
    std::ofstream outfile;
    outfile.open(DensityMatFile.c_str());
    outfile << std::setprecision(12);

    std::vector<std::vector<double> > MO(0);
    readFile(MO, MOFile, numMOs);

    int numBasis = MO.size();
    std::vector<std::vector<double> > densityMat(numBasis, std::vector<double>(numBasis,0.0));

    for(unsigned int i = 0; i < numBasis; ++i)
    {
        for(unsigned int j = 0; j < numBasis; ++j)
        {
            for(unsigned int k = 0; k < numMOs; ++k)
                densityMat[i][j] += MO[i][k]*MO[j][k];
            outfile << densityMat[i][j] << "\t";
        }
        outfile << std::endl;
    }
}

    
