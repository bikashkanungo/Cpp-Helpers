#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cerrno>
#include <cassert>


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
readLData(std::vector<std::vector<std::vector<double> > > & data,
        std::string fileName,
        const unsigned int lmax,
        const unsigned int numLinesPerL, 
        const unsigned int numHeaderLines,
        const int numColumns)
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
        //
        // skip the header lines
        //
        for(unsigned int i = 0; i < numHeaderLines; ++i)
            std::getline(readFile, readLine);

        for(unsigned int i = 0; i <= lmax; ++i)
        {

            //
            // skip the first line for each l
            //
            std::getline(readFile, readLine);

            for(unsigned int j = 0; j < numLinesPerL; ++j)
            {
                if(std::getline(readFile, readLine))
                {
                    //std::cout << i <<  " " << j << " " << readLine << std::endl; 
                    std::istringstream iss(readLine);
                    columnCount = 0; 

                    while(iss >> word && columnCount < numColumns)
                    {
                        if(!isNumber(rowData[columnCount],word))
                            std::cout << "Undefined number: " << word << std::endl;
                        columnCount++;
                        //rowData[columnCount++] = atof(word.c_str());
                    }

                    data[i].push_back(rowData); 
                }
                else
                    std::cout << "Error reading line." << std::endl;
            }
        }

    }

    else
        std::cout << "Error opening file." << std::endl;

    readFile.close();
    return;
}

int main()
{
    int lmax;
    int numLines;
    int numHeaderLines;
    std::string FHIFile, elementSymbol;
    std::cout << "Enter FHI File name:";
    std::cin >> FHIFile;
    std::cout << std::endl;
    std::cout <<"Enter element symbol:";
    std::cin >> elementSymbol;
    std::cout << std::endl;
    std::cout << "Enter lmax:";
    std::cin >> lmax;
    std::cout << std::endl;
    std::cout << "Number radial values (number of lines):";
    std::cin >> numLines;
    std::cout << std::endl;
    std::cout << "Number header lines:";
    std::cin >> numHeaderLines;
    std::cout << std::endl;

    std::vector<std::vector<std::vector<double> > > LData(lmax+1,std::vector<std::vector<double> >(0));
    readLData(LData,
                FHIFile,
                lmax,
                numLines,
                numHeaderLines,
                4);
    std::cout << "Read file" << std::endl;
    for(unsigned int i = 0; i <= lmax; ++i)
    {
        std::ostringstream ss;
        ss << i;
        std::string lstring = ss.str();
        std::ofstream outPsiFile;
        std::ofstream outPotFile;
        std::string psiFileName = elementSymbol + "_psi_l" + lstring + ".dat";
        std::string potFileName = elementSymbol + "_pot_l" + lstring + ".dat";
        outPsiFile.open(psiFileName.c_str());
        outPsiFile << std::setprecision(15);
        outPotFile.open(potFileName.c_str());
        outPotFile << std::setprecision(15);
        
        for(unsigned int j = 0; j < numLines; ++j)
        {
            const double radius = LData[i][j][1];
            const double psi = LData[i][j][2];
            const double pot = LData[i][j][3];
            outPsiFile << radius << " " << psi/radius << std::endl;
            outPotFile << radius << " " << pot << std::endl;

        }

        outPsiFile.close();
        outPotFile.close();
    
    }

}

