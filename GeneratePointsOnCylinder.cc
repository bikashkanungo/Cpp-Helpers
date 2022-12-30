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

	double r, zMin, zMax, NZ, NTheta;
	
	std::cout << "Enter radius: ";
	std::cin >> r;
	std::cout <<std::endl;
	
	std::cout << "Enter Z min: ";
	std::cin >> zMin;
	std::cout <<std::endl;
	
	std::cout << "Enter Z max: ";
	std::cin >> zMax;
	std::cout <<std::endl;
	
	std::cout << "Enter number points on axis: ";
	std::cin >> NZ;
	std::cout <<std::endl;
	
	std::cout << "Enter number of intervals on a circle: ";
	std::cin >> NTheta;
	std::cout <<std::endl;

	std::string outFilename;
	std::cout << "Enter output file name: ";
	std::cin >> outFilename;
	std::cout <<std::endl;
	std::ofstream outFile;
	outFile.open(outFilename.c_str());

	const double h = (zMax-zMin)/NZ;
	const double dtheta = 2*M_PI/NTheta;
	std::vector<std::vector<double> > points(0);
	for(unsigned i = 0; i <= NZ; ++i)
	{

		const double z = zMin + i*h;
		for(unsigned int j = 0; j <= NTheta; ++j)
		{

			std::vector<double> point(3);
			const double theta = j*dtheta;
			point[0] = r*std::cos(theta);
			point[1] = r*std::sin(theta);
			point[2] = z;
			points.push_back(point);
		}
	}

	const unsigned int nPoints = points.size();
	for(unsigned i = 0; i < nPoints; ++i)
	{
		for(unsigned int j = 0; j < 3; ++j)
		{

			outFile << points[i][j] << "\t";

		}
		outFile << std::endl;
	}

}
