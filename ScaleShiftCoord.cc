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
readData(std::vector<std::vector<double> > & data,
	 std::string filename)
{

    data.resize(0);
    std::ifstream readFile;
    readFile.open(filename.c_str());
    assert(readFile.is_open()); 
    std::string readLine;
    while(std::getline(readFile,readLine))
    {
	    std::istringstream lineString(readLine);
	    std::string word;
	    std::vector<double> rowData(0);
	    while(lineString >> word)
	    {

		    double val;
		    if(!isNumber(val,word))
			    std::cout << "Undefined value read in file: " << filename << std::endl;
		    else
			rowData.push_back(val);
	    }
	    data.push_back(rowData);
    }
    
    readFile.close();
}


int main()
{

	std::string inFilename;
	std::cout << "Enter input file name: ";
	std::cin >> inFilename;
	std::cout <<std::endl;

  //double L[3];
  //std::cout <<"Enter L along x: ";
	//std::cin >> L[0];
	//std::cout << std::endl;

  //std::cout <<"Enter L along y: ";
	//std::cin >> L[1];
	//std::cout << std::endl;
  //
	//std::cout <<"Enter L along z: ";
	//std::cin >> L[2];
	//std::cout << std::endl;
	
	std::string outFilename;
	std::cout << "Enter output file name: ";
	std::cin >> outFilename;
	std::cout <<std::endl;
	std::ofstream outFile;
	outFile.open(outFilename.c_str());

	double scale;
	std::cout << "Enter scale factor: ";
	std::cin >> scale;
	std::cout <<std::endl;

	char shiftChar;
	bool shiftFlag;
	std::cout << "Shift center of mass to origin? (Y/N): ";
	std::cin >> shiftChar;
	std::cout <<std::endl;

	if(shiftChar == 'Y')
		shiftFlag = true;
	else if(shiftChar == 'N')
		shiftFlag = false;
	else
		std::cout << "Invalid character for shift entered. Enter Y for yes or N for no." << std::endl;

	std::vector<std::vector<double> > data(0);
	readData(data, inFilename);
	const int numAtoms = data.size();

	std::vector<double> com(3,0.0);
	for(unsigned i = 0; i < numAtoms; ++i)
	{
		for(unsigned int j = 0; j < 3; ++j)
		{

			data[i][j] *= scale;
			com[j] += data[i][j];
		}
	}

	for(unsigned int j = 0; j < 3; ++j)
	{

		com[j] /= numAtoms;
	}

	
	for(unsigned i = 0; i < numAtoms; ++i)
	{
		for(unsigned int j = 0; j < 3; ++j)
		{

			if(shiftFlag)
				data[i][j] -= com[j];
			
			outFile << data[i][j] << "\t";

		}
		outFile << std::endl;
	}

}
