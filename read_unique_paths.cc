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


void readLines(std::vector<std::string> & lines, const std::string fileName)
{    
    std::ifstream readFile;
		lines.resize(0);

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
					lines.push_back(readLine);
        }
    }
    
    readFile.close();
}

void extractPath(const std::vector<std::string> & lines, std::vector<std::string> & paths)
{

	const unsigned int N = lines.size();
	paths.resize(N);
	for(unsigned int i = 0; i < N; ++i)
	{
		std::cout <<" length: " << lines[i].size() << std::endl;
		std::size_t first = lines[i].find_first_of("\\/");
		std::size_t last = lines[i].find_last_of("\\/");
		std::cout << "first: " << first << std::endl;
		std::cout << "last: " << last << std::endl;
		if(first!=std::string::npos && last!=std::string::npos) 
			paths[i] = lines[i].substr(first,last-first+1);
	}

}


int main()
{
    std::string infile;
    std::cout << "Enter input file name: ";
    std::cin >> infile;
    std::cout <<std::endl;

    std::string outfile;
    std::cout << "Enter output file name: ";
    std::cin >> outfile;
    std::cout << std::endl;


		std::vector<std::string> lines(0);
		std::vector<std::string> paths(0);
		readLines(lines, infile);
		extractPath(lines,paths);

    std::ofstream output;
    output.open(outfile.c_str());

		std::set<std::string> uniquePaths;
		for(unsigned int i = 0; i < paths.size(); ++i)
			uniquePaths.insert(paths[i]);

		std::set<std::string>::iterator it;
		for (it = uniquePaths.begin(); it != uniquePaths.end(); ++it) 
			output << *it << std::endl;
		
		output.close();
}
