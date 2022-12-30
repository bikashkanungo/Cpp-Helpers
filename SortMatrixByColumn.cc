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
#include <algorithm>
#include <utility>

bool compareByFirst(std::pair<int,double> p1, 
        std::pair<int,double> p2){
    return p1.first < p2.first;
}
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

void readDataFromFile(std::vector<std::vector<double> > & data, 
        const std::string fileName)
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
                    std::cout << "Undefined behavior: Value read in data file is not a number" << std::endl;
            }
            data.push_back(rowVals);
        }
        else
            std::cout << "Empty line found in data file: " << fileName << std::endl;
    }

    readFile.close();

}
    void
readNodalFieldFromFile(std::vector<double> & nodalVals,
        const std::string fileName,
        const int numLocalNodes,
        const std::vector<int> & globalNodeIds)

{

    nodalVals.resize(numLocalNodes);
    std::vector<std::vector<double> > data(0);
    readDataFromFile(data, fileName);
    const int numGlobalNodes = data.size();
    const int numCols = data[0].size();
    std::vector<std::pair<int,double> > allGlobalNodeIdsAndNodalVals(numGlobalNodes);
    for(unsigned int i = 0; i < numGlobalNodes; ++i)
    {
        const int globalNodeId = (int)(data[i][0]);
        const double nodalVal = data[i][numCols-1];
        std::pair<int,double> a(globalNodeId,nodalVal);
        allGlobalNodeIdsAndNodalVals[i] = a;
    }
    std::sort(allGlobalNodeIdsAndNodalVals.begin(), 
            allGlobalNodeIdsAndNodalVals.end(),
            compareByFirst);

    std::cout << "\n \n Printing sorted vals" << std::endl;
    for(unsigned int i = 0; i < numGlobalNodes; ++i)
    {
        std::cout << allGlobalNodeIdsAndNodalVals[i].first << "\t"  << 
            allGlobalNodeIdsAndNodalVals[i].second << std::endl;
    }
    for(unsigned int iLoc = 0; iLoc < numLocalNodes; ++iLoc)
    {
        const int globalId = globalNodeIds[iLoc];
        nodalVals[iLoc] = allGlobalNodeIdsAndNodalVals[globalId].second;
    }
}

void writeTestMatrix(const std::string fileName,
                     const int numRows,
                     const int numCols)
{

  std::ofstream fileToWrite;
  fileToWrite.open(fileName.c_str());
  srand(0);
  std::set<int> idsSet;
  std::vector<int> ids(0);
  while(idsSet.size() < numRows)
  {
    int val = (rand() % numRows);
    int sizePrev = idsSet.size();
    idsSet.insert(val);
    int sizeNext = idsSet.size();
    if(sizeNext > sizePrev)
      ids.push_back(val);
  }
  
  for(unsigned int i = 0; i < numRows; ++i)
  {
    fileToWrite <<  ids[i] << "\t";
    for(unsigned int j = 0; j < numCols-1; ++j)
    {
      fileToWrite << (rand()/(RAND_MAX + 1.0)) << "\t";
    }
    fileToWrite << std::endl;
  }
  fileToWrite.close();
}

int main()
{

    std::string fileName("TestMatrix");
    int numGlobalNodes = 20;
    //writeTestMatrix(fileName,
    //                numGlobalNodes,
    //                6);

    int numLocalNodes = 5;
    std::vector<int> globalNodeIdsToRead(numLocalNodes);
    srand(0);
    for(unsigned int i = 0; i < numLocalNodes; ++i)
    {
        globalNodeIdsToRead[i] = (rand() % numGlobalNodes);
        std::cout << globalNodeIdsToRead[i] << std::endl;
    }

    std::vector<double> nodalVals(0);
    readNodalFieldFromFile(nodalVals, fileName, numLocalNodes, globalNodeIdsToRead);

    std::cout << "\n \n Printing Vals" << std::endl;
    for(unsigned int i = 0; i < numLocalNodes; ++i)
    {
        std::cout << globalNodeIdsToRead[i] << "\t" << nodalVals[i] << std::endl;
    }

    return 0;
}
