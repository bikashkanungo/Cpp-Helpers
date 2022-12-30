#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>

#define NUM_CHILDREN 2

class Cell {

  public:
	Cell()
	  : d_name("") 
	{
	
	}

	Cell(const std::string name)
	  : d_name(name) 
	{
	  std::cout << "Creating cell: " << d_name << std::endl;
	}
	
	void setName(const std::string name) 
	{
	  d_name=name;
	}
	
	std::string getName() const  
	{
	  return d_name;
	}

	~Cell() {std::cout << "Deleting cell: " << d_name << std::endl;}
  private:
	std::string d_name;
};

class Triangulation {

  public:
	typedef std::vector<std::shared_ptr<Cell>>::iterator cellIterator;
	typedef std::vector<std::shared_ptr<Cell>>::const_iterator const_cellIterator;

	Triangulation(const std::vector<std::string> cellNames)
	  : d_cellPtrVector(0)
	{
	  for(unsigned int i = 0; i < cellNames.size(); ++i)
		d_cellPtrVector.push_back(std::make_shared<Cell>(cellNames[i]));
	}

	~Triangulation() 
	{
	}

	void printUseCount()
	{
	  for(unsigned int iCell = 0; iCell < d_cellPtrVector.size(); ++iCell)
	  {
		std::cout << "Use Count: " << d_cellPtrVector[iCell]->getName() << " "  
		  << d_cellPtrVector[iCell].use_count() << std::endl;
	  }
	}
	cellIterator begin() {return d_cellPtrVector.begin();}
	cellIterator end() {return d_cellPtrVector.end();}
	const_cellIterator begin() const  {return d_cellPtrVector.begin();}
	const_cellIterator end() const {return d_cellPtrVector.end();}

  private:
	std::vector<std::shared_ptr<Cell>> d_cellPtrVector;
};

class ParentToChildCellsManager {

  public:

	ParentToChildCellsManager()
	  : d_triangulationVector(0)
	{
	}

	std::vector<std::shared_ptr<const Cell>>
	  createChildCells(const std::string parentCellName)
	  {

		const unsigned int numChildren = NUM_CHILDREN;
		std::vector<std::string> childCellNames(numChildren);
		for(unsigned int i = 0; i < numChildren; ++i)
		  childCellNames[i] = parentCellName + "-" + std::to_string(i);

		auto triangulation = std::make_shared<Triangulation>(childCellNames);
		d_triangulationVector.push_back(triangulation);

		std::vector<std::shared_ptr<const Cell>> returnValue(0);

		Triangulation::cellIterator cellIter = triangulation->begin();
		unsigned int iCell = 0;
		for(; cellIter != triangulation->end(); ++cellIter)
		{
		  returnValue.push_back(*cellIter);
		}

		return returnValue;
	  }

	void popLast()
	{
	  d_triangulationVector[d_triangulationVector.size()-1]->printUseCount();
	  d_triangulationVector.pop_back();
	}

  private:
	std::vector<std::shared_ptr<Triangulation>> d_triangulationVector;
};


void recursive(const std::string parentName,
	const std::string globalCellName,
	const int  recursionLevel,
	const unsigned int maxRecursion,
	ParentToChildCellsManager  & parentToChildCellsManager)
{
  std::cout << "ParentName: " << parentName << std::endl; 
  std::string parentNameCopy(parentName);
  // remove the hyphens from the name
  parentNameCopy.erase(std::remove(parentNameCopy.begin(), parentNameCopy.end(), '-'), 
	  parentNameCopy.end());
  // convert the parentName string to int
  const unsigned int parentNameInt = std::stoi(parentNameCopy);

  if(recursionLevel < maxRecursion)
  {
	std::vector<std::shared_ptr<const Cell> > childCells =
	  parentToChildCellsManager.createChildCells(parentName);

	const unsigned int numChildren = childCells.size();
	// sub-divide if the parentName is even
	if(parentNameInt % 2 == 0)
	{
	  for(unsigned int iChild = 0; iChild < numChildren; ++iChild)
	  {
		recursive(childCells[iChild]->getName(),
			globalCellName,
			recursionLevel+1,
			maxRecursion,
			parentToChildCellsManager);

	  }
	}
	parentToChildCellsManager.popLast();
	for(unsigned int iCell = 0; iCell < childCells.size(); ++iCell)
	{
	  std::cout << "Use Count: " << childCells[iCell]->getName() << " "  
		<< childCells[iCell].use_count() << std::endl;
	}

  }

}

int main()
{

  ParentToChildCellsManager parentToChildCellsManager;
  int recursionLevel = 0;
  recursive("0","0", recursionLevel, 3, parentToChildCellsManager);
}
