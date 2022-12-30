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
void readMatrix(std::vector<std::vector<double> > & Mat, const std::string fileName)
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
            Mat.push_back(rowVals);
        }
        else
            std::cout << "Empty line found in density matrix file: " << fileName << std::endl;
    }
    
    readFile.close();

}

int main()
{
    std::string MatFileName1, MatFileName2;
    std::cout << "Enter filename for 1st Matrix: ";
    std::cin >> MatFileName1;
    std::cout << std::endl;
    std::cout << "Enter filename for 2nd Matrix: ";
    std::cin >> MatFileName2;
    std::cout << std::endl;
    std::vector<std::vector<double> > Mat1(0);
    std::vector<std::vector<double> > Mat2(0);

    readMatrix(Mat1,MatFileName1);
    readMatrix(Mat2,MatFileName2);
    
    const int rowSize = Mat1.size();
    const int colSize = Mat1[0].size();
    assert(rowSize == Mat2.size());
    assert(colSize == Mat2[0].size());
    std::cout << "RowSize: " << rowSize << " ColumnSize: " << colSize << std::endl;
    std::vector<std::vector<double> > Mat3(rowSize,std::vector<double>(colSize,0.0));
    for(unsigned int i = 0; i < rowSize; ++i)
    {
        for(unsigned int j = 0; j < colSize; ++j)
            Mat3[i][j] = Mat1[i][j]*Mat2[i][j];
    }

    double sum = 0.0;
    for(unsigned int i = 0; i < rowSize; ++i)
    {
        for(unsigned int j = 0; j < colSize; ++j)
           sum += Mat3[i][j];
    }

    std::cout << std::setprecision(18) <<"Sum: " << sum << std::endl;
}
