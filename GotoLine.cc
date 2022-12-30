#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <vector>
#include <limits>
#include <stdexcept>
std::ifstream & gotoLine(std::ifstream& file, unsigned int num)
{
  file.seekg(std::ios::beg);
  for(int i=0; i < num; ++i)
  {
	file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
  return file;
}

int main()
{
  unsigned int startLineNum = 6;
  unsigned int numLinesToRead = 3;
  unsigned int numColumns = 2;
  std::vector<std::vector<double>> data(0);
  std::vector<double> rowData(numColumns, 0.0);
  std::ifstream readFile;
  readFile.open("TestData.txt");
  std::string readLine;
  std::string word;
  int columnCount;
  if(readFile.is_open())
  {
	gotoLine(readFile, startLineNum);
	int lineCount = 0;
	while(std::getline(readFile, readLine) && lineCount < numLinesToRead)
	{
	  std::istringstream iss(readLine);
	  columnCount = 0; 
	  while(iss >> word && columnCount < numColumns)
	  {
		rowData[columnCount] = atof(word.c_str());
		columnCount++;
	  }
	  data.push_back(rowData);
	  lineCount++;
	}
  }
  else
  {
	throw std::invalid_argument("File not found");
  }

  readFile.close();

  for(int i = 0; i < numLinesToRead; ++i)
  {
	for(unsigned int j = 0; j < numColumns; ++j)
	  std::cout << data[i][j] << " ";
	std::cout << std::endl;
  } 

  return 0;
}

