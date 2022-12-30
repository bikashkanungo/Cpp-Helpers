#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>

/***
 * It creates a non-periodic diamond lattice.
 * The idea is to first build a simple cubic, 
 * add the face-centered atoms and then add the
 * four interior atoms. This aavoid duplication
 * of atoms.
 *
 ***/

int main()
{

  double lc;
  int nx,ny,nz;
   
  std::ofstream file;
  std::string filename;
  std::string element;

  std::cout << "Enter the element symbol: ";
  std::cin >> element;
  std::cout << "\nEnter lattice constant in bohr: ";
  std::cin >> lc;
  std::cout << "\nEnter NX: ";
  std::cin >> nx;
  std::cout << "\nEnter NY: ";
  std::cin >> ny;
  std::cout << "\nEnter NZ: ";
  std::cin >> nz;

  filename = element;
  std::stringstream s;
  s << nx;
  filename += "_" + s.str() + "x";
  s.str(std::string());
  s << ny;
  filename += s.str() + "x";
  s.str(std::string());
  s << nz;
  filename += s.str();
  std::cout << filename << std::endl;

  file.open(filename.c_str());

  std::vector<double> origin(3);
  origin[0] = nx/2.0;
  origin[1] = ny/2.0;
  origin[2] = nz/2.0;

  std::vector<std::vector<double> > atoms(0);

  //
  // add the simple cubic atoms
  // 
  for(unsigned int i = 0; i <= nx; ++i)
  {

    for(unsigned int j = 0; j <= ny;  ++j)
    { 
  
      for(unsigned int k = 0; k <= nz; ++k)
      {

        std::vector<double> atom(3);
        atom[0] = i;
        atom[1] = j;
        atom[2] = k;
        atoms.push_back(atom);

      }

    }

  }
  
  //
  // add the face-centered atoms
  //
  for(unsigned int i = 0; i <= nx; ++i)
  {

    for(unsigned int j = 0; j < ny; ++j)
    {
	
      for(unsigned int k = 0; k < nz; ++k)
      {

        std::vector<double> atom(3);
    	atom[0] = i;
    	atom[1] = (0.5+j);
	atom[2] = (0.5+k);
	atoms.push_back(atom);

      }

    }

  }

  for(unsigned int i = 0; i < nx; ++i)
  {

    for(unsigned int j = 0; j <= ny; ++j)
    {
	
      for(unsigned int k = 0; k < nz; ++k)
      {

        std::vector<double> atom(3);
    	atom[0] = (0.5+i);
    	atom[1] = j;
	atom[2] = (0.5+k);
	atoms.push_back(atom);

      }

    }

  }

  for(unsigned int i = 0; i < nx; ++i)
  {

    for(unsigned int j = 0; j < ny; ++j)
    {
	
      for(unsigned int k = 0; k <= nz; ++k)
      {

        std::vector<double> atom(3);
    	atom[0] = (0.5+i);
    	atom[1] = (0.5+j);
	atom[2] = k;
	atoms.push_back(atom);

      }

    }

  }

  //
  // add the four interior atoms
  //
    for(unsigned int i = 0; i < nx; ++i)
  {

    for(unsigned int j = 0; j < ny; ++j)
    {
	
      for(unsigned int k = 0; k < nz; ++k)
      {

	std::vector<double> atom(3);

	//
	// interior atom 1
	//
        atom[0] = (0.75+i);
    	atom[1] = (0.75+j);
	atom[2] = (0.75+k);
	atoms.push_back(atom);

	//
	// interior atom 2
	//
        atom[0] = (0.75+i);
    	atom[1] = (0.25+j);
	atom[2] = (0.25+k);
	atoms.push_back(atom);	

 	//
	// interior atom 2
	//
        atom[0] = (0.25+i);
    	atom[1] = (0.75+j);
	atom[2] = (0.25+k);
	atoms.push_back(atom);		

	//
	// interior atom 2
	//
        atom[0] = (0.25+i);
    	atom[1] = (0.25+j);
	atom[2] = (0.75+k);
	atoms.push_back(atom);	

      }

    }

  }

  std::cout << "\nNumber of atoms: " << atoms.size() << std::endl;;  
 
  //
  // shift the origin to (nx/2, ny/2, nz/2) and
  // multiply the fractional coordinates with lattice constant
  //
  for(unsigned int i = 0; i < atoms.size(); ++i)
  {

    for(unsigned int j = 0; j < 3; ++j)
    {
 
      atoms[i][j] -= origin[j];      
      atoms[i][j] *= lc;

    }

    file << "{" << element << ", " << atoms[i][0] << ", " << atoms[i][1] << ", " << atoms[i][2] << " }" << std::endl;

  }

}

  
