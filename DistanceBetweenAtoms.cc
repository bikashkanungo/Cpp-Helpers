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

double getDistance(const std::vector<double> & x,
		const std::vector<double> & y)
{

	double returnValue = 0.0;

	for(unsigned int i = 0; i < x.size(); ++i)
		returnValue += pow(x[i]-y[i], 2.0);

	returnValue = sqrt(returnValue);
	return returnValue;
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

	std::vector<std::vector<double> > data(0);
	readData(data, inFilename);
	const int numAtoms = data.size();
	for(unsigned int i = 0; i < numAtoms; ++i)
	{

		for(unsigned int j = i+1; j < numAtoms; ++j)
		{
		  std::vector<double> dx(3);
			for(unsigned int k = 0; k < 3; ++k)
			{
			  const double d1 = std::fabs(data[i][k]-data[j][k]);
			  dx[k] = d1;
				//const double d2 = std::fabs(d1-L[k]);
				//dx[k] = d1 < d2 ? d1 : d2; 
			}

			double dist = 0.0;;
			for(unsigned int k = 0; k < 3; ++k)
			{	
				dist += dx[k]*dx[k];
			}

			dist =std::sqrt(dist);

			//double dist = getDistance(data[i], data[j]);
			outFile << i+1 << "\t" << j+1 << "\t" << dist << std::endl;			
		}	
	}

	outFile.close();
}
