#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <mm_malloc.h>
#define dim 3
#define p 6
#define MIN_Q_1D 30
#define MAX_Q_1D 60

typedef unsigned int size_type;

  size_type
sizeTypePow(size_type base, size_type e)
{
  size_type result = 1;
  for (;;)
  {
	if (e & 1)
	  result *= base;
	e >>= 1;
	if (!e)
	  break;
	base *= base;
  }
  return result;
}

void
get1DIds(const size_type I,
	std::vector<size_type> & ids)
{
  size_type index = I;
  for(int iDim = dim-1; iDim >= 0; --iDim)
  {
	ids[iDim] = index/sizeTypePow(p+1,iDim);
	index = index - ids[iDim]*(sizeTypePow(p+1,iDim));
  }
}


int main()
{

  size_type NQ3D = 100000;
  size_type nCells = 20;
  //std::vector<size_type> QMap(NQ3D*dim);
  size_type * QMap = (size_type *)_mm_malloc(
						NQ3D*dim*sizeof(size_type), 64);
  std::vector<size_type> NQ1Ds(dim);
  size_type nNode = sizeTypePow(p+1,dim);
  size_type shapeFunction1DStorageSize = 0;
  std::vector<size_type> shapeFunction1DStartIds(dim,0);
  for(size_type iDim = 0; iDim < dim; ++iDim)
  {
	NQ1Ds[iDim] = (size_type)(MIN_Q_1D + (MAX_Q_1D-MIN_Q_1D)*(rand()/(RAND_MAX+0.0)));
	shapeFunction1DStartIds[iDim] = shapeFunction1DStorageSize;
	shapeFunction1DStorageSize += NQ1Ds[iDim]*(p+1);
  }

  for(size_type iQ = 0; iQ < NQ3D; ++iQ)
  {
    for(size_type iDim = 0; iDim < dim; ++iDim)
	{
	  QMap[iQ*dim + iDim] = (size_type)((NQ1Ds[iDim]-1)*(rand()/(RAND_MAX+0.0))); 
	}
  }

  std::cout << "shapeFunction1DStorageSize: " <<  shapeFunction1DStorageSize << std::endl; 
  double * shapeFunctionVals1D = (double *)_mm_malloc(shapeFunction1DStorageSize*sizeof(double), 64);
  double * shapeVals = (double *)_mm_malloc(nNode*NQ3D*sizeof(double), 64);
  std::vector<std::vector<double>> shapeVals2DVec(nNode, std::vector<double>(NQ3D,1.0));
  std::fill(&shapeVals[0], &shapeVals[nNode*NQ3D], 1.0);
  std::vector<size_type> ids1D(dim);
  auto start1 = std::chrono::high_resolution_clock::now();  
#pragma omp simd
  for(size_type iCell = 0; iCell < nCells; ++iCell)
  {
	for(size_type iNode = 0; iNode < nNode; ++iNode)
	{
	  get1DIds(iNode, ids1D);
	  //std::vector<double> & shapeValsINode = shapeVals2DVec[iNode];
	  for(size_type iQ = 0; iQ < NQ3D; ++iQ)
	  {
		const size_type IQIndex = iNode*NQ3D + iQ;
		double dummy = 1.0;
		for(size_type iDim = 0; iDim < dim; ++iDim)
		{ 
		  const size_type index = shapeFunction1DStartIds[iDim] + ids1D[iDim]* QMap[iQ*dim + iDim];
		  shapeVals[IQIndex] *= shapeFunctionVals1D[index]; 
		  //dummy *= shapeFunctionVals1D[index];
		  //shapeValsINode[iQ] = shapeFunctionVals1D[index]; 
		}
	  }
	}
  }

  auto time1 = std::chrono::duration_cast<std::chrono::nanoseconds>
	(std::chrono::high_resolution_clock::now() - start1).count();  

  auto start2 = std::chrono::high_resolution_clock::now();  

  for(size_type iQ = 0; iQ < NQ3D; ++iQ)
  {
	for(size_type iNode = 0; iNode < nNode; ++iNode)
	{
	  get1DIds(iNode, ids1D);
	  const size_type IQIndex = iQ*nNode + iNode;
	  for(size_type iDim = 0; iDim < dim; ++iDim)
	  { 
		const size_type index = shapeFunction1DStartIds[iDim] + ids1D[iDim]* QMap[iQ*dim + iDim];
		shapeVals[IQIndex] *= shapeFunctionVals1D[index]; 
	  }
	}
  }

  auto time2 = std::chrono::duration_cast<std::chrono::nanoseconds>
	(std::chrono::high_resolution_clock::now() - start2).count();  

  auto start3 = std::chrono::high_resolution_clock::now();  

  for(size_type iDim = 0; iDim < dim; ++iDim)
  {
	for(size_type iNode = 0; iNode < nNode; ++iNode)
	{
	  get1DIds(iNode, ids1D);
	  for(size_type iQ = 0; iQ < NQ3D; ++iQ)
	  {
		const size_type index = shapeFunction1DStartIds[iDim] + ids1D[iDim]* QMap[iQ*dim + iDim];
		const size_type IQIndex = iNode*NQ3D + iQ;
		shapeVals[IQIndex] *= shapeFunctionVals1D[index]; 
	  }
	}
  }

  auto time3 = std::chrono::duration_cast<std::chrono::nanoseconds>
	(std::chrono::high_resolution_clock::now() - start3).count();  

  auto start4 = std::chrono::high_resolution_clock::now();  
  for(size_type iDim = 0; iDim < dim; ++iDim)
  {
	for(size_type iQ = 0; iQ < NQ3D; ++iQ)
	{
	  const size_type index = shapeFunction1DStartIds[iDim] + ids1D[iDim]* QMap[iQ*dim + iDim];
	  for(size_type iNode = 0; iNode < nNode; ++iNode)
	  {
		get1DIds(iNode, ids1D);
		const size_type IQIndex = iNode*NQ3D + iQ;
		shapeVals[IQIndex] *= shapeFunctionVals1D[index]; 
	  }
	}
  }

  auto time4 = std::chrono::duration_cast<std::chrono::nanoseconds>
	(std::chrono::high_resolution_clock::now() - start4).count();  

  _mm_free(QMap);
  _mm_free(shapeFunctionVals1D);
  _mm_free(shapeVals);

  std::cout << "Time1: " << time1 << std::endl;
  std::cout << "Time2: " << time2 << std::endl;
  std::cout << "Time3: " << time3 << std::endl;
  std::cout << "Time4: " << time4 << std::endl;
}
