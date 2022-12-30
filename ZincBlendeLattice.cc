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


	void
shiftAndScale(std::vector<double> & coord,
		std::vector<double> & origin,
		const double scale)
{

	for(unsigned int i = 0; i < coord.size(); ++i)
	{
		coord[i] -= origin[i];
		coord[i] *= scale;
	}

}

// y = A*x
void matVecMult(std::vector<std::vector<double> > & A,
		const std::vector<double>  & x,
		std::vector<double>  & y)
{
	int M = A.size();
	int N = A[0].size();
	y.resize(M, 0.0);
	for(unsigned int i = 0; i < M; ++i)
	{
		for(unsigned int j = 0; j < N; ++j)
		{
			y[i] += A[i][j]*x[j];
		}
	}

}

void
addVecs(const std::vector<double> & x,
		std::vector<double> & y)
{
	for(unsigned int i = 0; i < x.size(); ++i)
		y[i] += x[i];
}
	void
printAtomicCoord(const std::vector<double> & coord,
		std::string symbol,
		std::ofstream & file)
{
	file << symbol << "\t" << coord[0] << "\t" << coord[1] << "\t" << coord[2] << std::endl;
}

int main()
{

	double lc;
	int nx,ny,nz;

	std::ofstream file;
	std::string filename;
	std::string element1;
	std::string element2;

	std::cout << "Enter the element symbol for FCC atoms: ";
	std::cin >> element1;
	std::cout << "Enter the element symbol for FCC atoms: ";
	std::cin >> element2;
	std::cout << "\nEnter lattice constant: ";
	std::cin >> lc;
	std::cout << "\nEnter NX: ";
	std::cin >> nx;
	std::cout << "\nEnter NY: ";
	std::cin >> ny;
	std::cout << "\nEnter NZ: ";
	std::cin >> nz;

	filename = element1+element2;
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
	origin[0] = 0.0;//nx/2.0;
	origin[1] = 0.0;//ny/2.0;
	origin[2] = 0.0;//nz/2.0;

	std::vector<std::vector<double> > primitiveVectorsMat(3, std::vector<double>(3,0.0));
	primitiveVectorsMat[0][0] = 0.0; 
	primitiveVectorsMat[1][0] = 0.5; 
	primitiveVectorsMat[2][0] = 0.5;

	primitiveVectorsMat[0][1] = 0.5; 
	primitiveVectorsMat[1][1] = 0.0; 
	primitiveVectorsMat[2][1] = 0.5;

	primitiveVectorsMat[0][2] = 0.5; 
	primitiveVectorsMat[1][2] = 0.5; 
	primitiveVectorsMat[2][2] = 0.0;

	int numBasisAtoms = 2;
	std::vector<std::string> basisAtomStrings(numBasisAtoms);
	basisAtomStrings[0] = element1;
	basisAtomStrings[1] = element2;
	std::vector<std::vector<double> > basisAtomPositions(numBasisAtoms, std::vector<double>(3,0.0));
	basisAtomPositions[0][0] = 0.0;
	basisAtomPositions[0][1] = 0.0;
	basisAtomPositions[0][2] = 0.0;

	basisAtomPositions[1][0] = 0.25;
	basisAtomPositions[1][1] = 0.25;
	basisAtomPositions[1][2] = 0.25;

	std::vector<std::vector<double> > atoms(0);
	std::vector<std::string> atomStrings(0);
	for(unsigned int i = 0; i < nx; ++i)
	{

		for(unsigned int j = 0; j < ny;  ++j)
		{ 

			for(unsigned int k = 0; k < nz; ++k)
			{

				std::vector<double> indices(3);
				std::vector<double> coord(3);
				indices[0] = i;
				indices[1] = j;
				indices[2] = k;
				matVecMult(primitiveVectorsMat, indices, coord);
				for(unsigned int iBasisAtom = 0; iBasisAtom < numBasisAtoms; ++iBasisAtom)
				{
					std::vector<double> atom = coord;
					addVecs(basisAtomPositions[iBasisAtom], atom);	
					atoms.push_back(atom);
					atomStrings.push_back(basisAtomStrings[iBasisAtom]);
				}

			}

		}

	}

	const int numAtoms = atoms.size();
	for(unsigned int i = 0; i < numAtoms; ++i)
	{

	  for(unsigned int j = 0; j < 3; ++j)
	  {

	    atoms[i][j] -= origin[j];      
	    atoms[i][j] *= lc;

	  }
	}
	std::cout << "\nNumber of atoms: " << numAtoms << std::endl;;  
	for(unsigned int i = 0; i < numAtoms; ++i)
		file << atomStrings[i] << "\t" << atoms[i][0] << "\t" << atoms[i][1] << "\t" << atoms[i][2] << std::endl;
	////
	//// shift the origin to (nx/2, ny/2, nz/2) and
	//// multiply the fractional coordinates with lattice constant
	////

	//  file << "{" << element << ", " << atoms[i][0] << ", " << atoms[i][1] << ", " << atoms[i][2] << " }" << std::endl;

	//}

		}


