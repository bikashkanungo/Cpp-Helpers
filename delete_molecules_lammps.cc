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
#include <algorithm>

class Compare {
public:
    Compare(int col) : col_(col) {}
    bool operator()(const std::vector<double>& lhs, const std::vector<double>& rhs) 
    {
        return lhs[col_] < rhs[col_];
    }
private:
    int col_; 
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
readData(std::vector<std::vector<double> > & data, const std::string fileName)
{

	std::ifstream readFile;
	std::string readLine;
	readFile.open(fileName.c_str());
	assert(readFile.is_open()); 
	data.resize(0);
	int lineCount = 0;
	while(std::getline(readFile,readLine))
	{
		std::istringstream lineString(readLine);
		std::string word;
		std::vector<double> rowData(0);
		while(lineString >> word)
		{
			double val;
			if(!isNumber(val,word))
				std::cout << "Undefined value read for the file: " << fileName << std::endl;
			else
				rowData.push_back(val);
		}
		data.push_back(rowData);
	}


	readFile.close();
}

bool isCoordInBox(const std::vector<double> & coord,
		  const std::vector<std::vector<double> > & boxBounds)
{
	bool returnValue = true;
	for(unsigned int i = 0; i < 3; ++i)
	{
		if(coord[i] < boxBounds[i][0] || coord[i] > boxBounds[i][1])
		{
			returnValue = false;
			break;
		}
	}

	return returnValue;

}

int main()
{
	std::vector<std::vector<double> > boxBounds(3, std::vector<double>(2,0.0));
	std::string inFilename;
	std::string outFilename;

	std::ofstream outFile;

	std::cout << "Enter input filename: ";
	std::cin >> inFilename;
	std::cout << "Enter Xlo: ";
	std::cin >> boxBounds[0][0];
	std::cout << std::endl;
	std::cout << "Enter Xhi: ";
	std::cin >> boxBounds[0][1];
	std::cout << std::endl;
	std::cout << "Enter Ylo: ";
	std::cin >> boxBounds[1][0];
	std::cout << std::endl;
	std::cout << "Enter Yhi: ";
	std::cin >> boxBounds[1][1];
	std::cout << std::endl;
	std::cout << "Enter Zlo: ";
	std::cin >> boxBounds[2][0];
	std::cout << std::endl;
	std::cout << "Enter Zhi: ";
	std::cin >> boxBounds[2][1];
	std::cout << std::endl;

	std::cout << "\nEnter output filename: ";
	std::cin >> outFilename;
	outFile.open(outFilename.c_str());

	int dataCoordIndexOffset = 4;
	int dataMolIdIndex = 1;
	
	std::vector<std::vector<double> > data(0);
	readData(data, inFilename);
	std::sort(data.begin(), data.end(), Compare(dataMolIdIndex));
		
	int numAtoms = data.size();
	//for(unsigned int i = 0; i < numAtoms; ++i)
	//{		
	//	for(unsigned int k = 0; k < data[0].size(); ++k)
	//		outFile << "\t" << data[i][k];
	//	outFile << std::endl;
	//}
	
	std::set<int> excludedMolIds;
	for(unsigned int i = 0; i < numAtoms; ++i)
	{
		std::vector<double> coord(3);
		for(unsigned int j = 0; j < 3; ++j)
			coord[j] = data[i][j+dataCoordIndexOffset];
		if(!isCoordInBox(coord, boxBounds))
		{
			int molId = (int)(data[i][dataMolIdIndex]);
			excludedMolIds.insert(molId);
		}
	}
	
	int numExcludedMols = excludedMolIds.size();
	int startIndex = 0;
        int atomCount = 0;
	std::set<int>::const_iterator iter;
	int excludedMolCount = 0;
	for(iter = excludedMolIds.begin(); iter != excludedMolIds.end(); ++iter)
	{
		int molId = *iter;//excludedMolIds[i];
		for(unsigned int j = startIndex; j < numAtoms; ++j)
		{
			if(molId == data[j][dataMolIdIndex])
			{
				int k = j;
				while(data[k][dataMolIdIndex] == molId) 
				{
					k++;
				}
				startIndex = k;
				break;
			}
			else
			{
				data[j][0] = atomCount + 1;
				data[j][dataMolIdIndex] -= excludedMolCount;
				for(unsigned int k = 0; k < data[0].size(); ++k)
					outFile << data[j][k] << "\t";
				outFile << std::endl;
				atomCount++;
			}	
		}
		excludedMolCount++;
	}	
	outFile.close();
}
