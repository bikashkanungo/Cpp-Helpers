
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

	std::cout << "\nEnter NX: ";
	std::cin >> nx;
	std::cout << "\nEnter NY: ";
	std::cin >> ny;
	std::cout << "\nEnter NZ: ";
	std::cin >> nz;

	std::cout << "\nEnter output filename: ";
	std::cin >> filename;
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

	std::vector<double> origin(3,0.0);
	//origin[0] = nx/2.0;
	//origin[1] = ny/2.0;
	//origin[2] = nz/2.0;

	std::vector<std::vector<double> > primitiveVectorsMat(3, std::vector<double>(3));
	primitiveVectorsMat[0][0] = 2.0; 
	primitiveVectorsMat[1][0] = 0.0; 
	primitiveVectorsMat[2][0] = 0.0;

	primitiveVectorsMat[0][1] = 0.0; 
	primitiveVectorsMat[1][1] = 2.0; 
	primitiveVectorsMat[2][1] = 0.0;

	primitiveVectorsMat[0][2] = 0.0; 
	primitiveVectorsMat[1][2] = 0.0; 
	primitiveVectorsMat[2][2] = 2.0;

	int numBasisAtoms = 3;
	std::vector<std::string> basisAtomStrings(numBasisAtoms);
	basisAtomStrings[0] = "H";
	basisAtomStrings[1] = "H";
	basisAtomStrings[2] = "O";

	std::vector<int> basisAtomTypes(numBasisAtoms);
	basisAtomTypes[0] = 1;
	basisAtomTypes[1] = 1;
	basisAtomTypes[2] = 2;

	std::vector<double> basisAtomCharges(numBasisAtoms);
	basisAtomCharges[0] = 0.4238;
	basisAtomCharges[1] = 0.4238;
	basisAtomCharges[2] = -0.8476;


	std::vector<std::vector<double> > basisAtomPositions(numBasisAtoms, std::vector<double>(3,0.0));
	basisAtomPositions[0][0] = 0.0;
	basisAtomPositions[0][1] = 0.0;
	basisAtomPositions[0][2] = 0.0;

	basisAtomPositions[1][0] = 1.633;
	basisAtomPositions[1][1] = 0.0;
	basisAtomPositions[1][2] = 0.0;

	basisAtomPositions[2][0] = 0.8165;
	basisAtomPositions[2][1] = 0.557735;
	basisAtomPositions[2][2] = 0.0;

	int numBasisBonds = 2;
	std::vector<int> basisBondTypes(numBasisBonds);
	basisBondTypes[0] = 1;
	basisBondTypes[1] = 1;

	std::vector<std::vector<int> > basisBondAtomType(numBasisBonds, std::vector<int>(2));
	basisBondAtomType[0][0] = 0;
	basisBondAtomType[0][1] = 2;

	basisBondAtomType[1][0] = 1;
	basisBondAtomType[1][1] = 2;

	int numBasisAngles = 1;
	std::vector<int> basisAngleTypes(numBasisAngles);
	basisAngleTypes[0] = 1;
	std::vector<std::vector<int> > basisAngleAtomType(numBasisAngles, std::vector<int>(3));
	basisAngleAtomType[0][0] = 0;
	basisAngleAtomType[0][1] = 2;
	basisAngleAtomType[0][2] = 1;


	std::vector<std::vector<double> > atoms(0);
	std::vector<std::string> atomStrings(0);
	std::vector<int> atomIds(0);
	std::vector<int> atomMolIds(0);
	std::vector<double> atomCharges(0);
	std::vector<int> atomTypes(0);

	std::vector<int> bondIds(0);
	std::vector<int> bondTypes(0);
	std::vector<std::vector<int> > bondAtomIds(0);

	std::vector<int> angleIds(0);
	std::vector<int> angleTypes(0);
	std::vector<std::vector<int> > angleAtomIds(0);
	//
	// add the simple cubic atoms
	//
	int molIdCount = 0; 
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
					//atomStrings.push_back(basisAtomStrings[iBasisAtom]);
					atomIds.push_back(molIdCount*numBasisAtoms + iBasisAtom + 1);
					atomMolIds.push_back(molIdCount + 1);
					atomTypes.push_back(basisAtomTypes[iBasisAtom]);
					atomCharges.push_back(basisAtomCharges[iBasisAtom]);
					atoms.push_back(atom);
				}

				// bond info
				for(unsigned int iBasisBond = 0; iBasisBond < numBasisBonds; ++iBasisBond)
				{
					bondIds.push_back(molIdCount*numBasisBonds + iBasisBond + 1);
					bondTypes.push_back(basisBondTypes[iBasisBond]);
					int atomIdOffset = molIdCount*numBasisAtoms;
					std::vector<int> atomIds(2);
					atomIds[0] = atomIdOffset + basisBondAtomType[iBasisBond][0] + 1; 
					atomIds[1] = atomIdOffset + basisBondAtomType[iBasisBond][1] + 1;
					bondAtomIds.push_back(atomIds); 
				}

				// angle info
				//
				for(unsigned int iBasisAngle = 0; iBasisAngle < numBasisAngles; ++iBasisAngle)
				{
					angleIds.push_back(molIdCount*numBasisAngles + iBasisAngle + 1);
					angleTypes.push_back(basisAngleTypes[iBasisAngle]);
					int atomIdOffset = molIdCount*numBasisAtoms;
					std::vector<int> atomIds(3);
					atomIds[0] = atomIdOffset + basisAngleAtomType[iBasisAngle][0] + 1;
					atomIds[1] = atomIdOffset + basisAngleAtomType[iBasisAngle][1] + 1;
					atomIds[2] = atomIdOffset + basisAngleAtomType[iBasisAngle][2] + 1;
					angleAtomIds.push_back(atomIds);
				}	

				molIdCount++;

			}

		}

	}

	unsigned int numAtoms = atomIds.size();
	unsigned int numBonds = bondIds.size();
	unsigned int numAngles = angleIds.size();
	file << "\n Atoms\n" << std::endl;
	for(unsigned int i = 0; i < numAtoms; ++i)
	{
		file << atomIds[i] << "\t" << atomMolIds[i] << "\t" << atomTypes[i] << "\t" << 
			atomCharges[i] << "\t" << atoms[i][0] << "\t" << atoms[i][1] << "\t" << atoms[i][2] << std::endl;
	}
	
	file << "\n Bonds\n" << std::endl;
	for(unsigned int i = 0; i < numBonds; ++i)
	{
		file << bondIds[i] << "\t" << bondTypes[i] << "\t" << bondAtomIds[i][0] << "\t" << bondAtomIds[i][1] << std::endl;
	}

	file << "\n Angles\n" << std::endl;
	for(unsigned int i = 0; i < numAngles; ++i)
	{
		file << angleIds[i] << "\t" << angleTypes[i] << "\t" << angleAtomIds[i][0] << "\t" << 
			angleAtomIds[i][1] << "\t" << angleAtomIds[i][2] << std::endl;
	}

	file.close();
}
