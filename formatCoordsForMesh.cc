#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>

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

    std::ifstream readFile;
    readFile.open(infile.c_str());

    std::ofstream output;
    output.open(outfile.c_str());

    if(readFile.is_open())
    {

        std::string readLine;
        while(!readFile.eof())
        {

            //
            // get the next line 
            //
            std::getline(readFile, readLine);
            std::istringstream lineString(readLine);
            std::string dummyString;
            int count = 0;
            while(lineString >> dummyString)
            {

                if(count == 0)
                    output << "{";

                output << dummyString;

                if(count < 3)
                    output << ", ";
                if(count == 3)
                    output << "}" << std::endl;

                count++;

            }
        }
    }

    readFile.close();
    output.close();
}



