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
#include <map>
#include <iomanip>
#include <algorithm>

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

void
readInts(std::vector<int> & data, std::string fileName)
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
			int val;
			if(!isInteger(val,word))
				std::cout << "Undefined value read for the file: " << fileName << std::endl;
			else
				data.push_back(val);
		}
	}

	readFile.close();
}

int main()
{
	std::map<int, std::string> customMap;
	customMap.insert(std::pair<int,std::string>(1,"H"));
	customMap.insert(std::pair<int,std::string>(2,"O"));
	std::string inFilename;
	std::string outFilename;
	std::cout << "Input file: ";
	std::cin >> inFilename;
	std::cout << std::endl;
	std::cout << "Output file: ";
	std::cin >> outFilename;
	std::ofstream outFile;
	outFile.open(outFilename.c_str());
	
	std::vector<int> data(0);
	readInts(data, inFilename);
	for(unsigned int i = 0; i < data.size(); ++i)
	{
		int x = data[i];
		std::map<int,std::string>::const_iterator pos = customMap.find(x);
		if(pos != customMap.end())
			outFile << pos->second << std::endl;
		else
			std::cout << "ERROR: No string found for index: " << x << std::endl;	
	}
}
