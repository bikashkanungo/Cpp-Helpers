#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <ctime>

/***
 * This code performs a coordinate transformation between two coordinate 
 * systems with common origin. 
 *
 * Let u1, u2, u3 be the unit vectors of the first coordinate system.
 * Let v1, v2, v3 be the unit vectors of the second coordinate system.
 * A position vector x in the first coordinate system is given as
 * x = a1*u1 + a2*u2 + a3*u3, where a={a1,a2,a3} 
 * are the coordinates in the u basis.
 * Similarly, for the second coordinate system, x = b1*v1 + b2*v2 + b3*v3, 
 * where b={b1,b2,b3} are the coordinates in the v basis. 
 * Thus, x = U*a = V*b, where U={u1,u2,u3} and V={v1,v2,v3} are matrices 
 * containing the u and v basis vectors, respectively. 
 * Thus, given the coordinates "a" in u basis, the coordinates "b" in the v basis
 * are given by, b = V^{-1}*U*a  
 ***/

// Function to get cofactor of A[p][q] in temp[][]. n is current 
// dimension of A[][] 
void getCofactor(std::vector<std::vector<double> >  & A, 
		std::vector<std::vector<double> > & temp, 
		int p, 
		int q, 
		int n) 
{ 
	int i = 0, j = 0; 

	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++) 
	{ 
		for (int col = 0; col < n; col++) 
		{ 
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != p && col != q) 
			{ 
				temp[i][j++] = A[row][col]; 

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1) 
				{ 
					j = 0; 
					i++; 
				} 
			} 
		} 
	} 
} 

/* Recursive function for finding determinant of matrix. 
   n is current dimension of A[][]. */
double determinant(std::vector<std::vector<double> >  & A, int n) 
{ 
	double D = 0.0; // Initialize result 

	//  Base case : if matrix contains single element 
	if (n == 1) 
		return A[0][0]; 

	int N = A.size();
	std::vector<std::vector<double> > temp(N, std::vector<double>(N,0.0)); // To store cofactors 

	double sign = 1;  // To store sign multiplier 

	// Iterate for each element of first row 
	for (int f = 0; f < n; f++) 
	{ 
		// Getting Cofactor of A[0][f] 
		getCofactor(A, temp, 0, f, n); 
		D += sign * A[0][f] * determinant(temp, n - 1); 

		// terms are to be added with alternate sign 
		sign = -sign; 
	} 

	return D; 
}

// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(std::vector<std::vector<double> > & A,
	     std::vector<std::vector<double> > & adj) 
{ 
	int N = A.size();
	if (N == 1) 
	{ 
		adj[0][0] = 1; 
		return; 
	} 

	// temp is used to store cofactors of A[][] 
	double sign = 1.0; 
	std::vector<std::vector<double> > temp(N, std::vector<double>(N,0.0)); 

	for (int i=0; i<N; i++) 
	{ 
		for (int j=0; j<N; j++) 
		{ 
			// Get cofactor of A[i][j] 
			getCofactor(A, temp, i, j, N); 

			// sign of adj[j][i] positive if sum of row 
			// and column indexes is even. 
			sign = ((i+j)%2==0)? 1.0: -1.0; 

			// Interchanging rows and columns to get the 
			// transpose of the cofactor matrix 
			adj[j][i] = (sign)*(determinant(temp, N-1)); 
		} 
	} 
} 

// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool inverse(std::vector<std::vector<double> > & A, 
	     std::vector<std::vector<double> > & inverse) 
{ 
	// Find determinant of A[][] 
	int N = A.size();
	double det = determinant(A, N); 
	if (det == 0) 
	{ 
		std::cout << "Singular matrix, can't find its inverse"; 
		return false; 
	} 

	// Find adjoint 
	std::vector<std::vector<double> >  adj(N, std::vector<double>(N,0.0)); 
	adjoint(A, adj); 

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
	for (int i=0; i<N; i++) 
		for (int j=0; j<N; j++) 
			inverse[i][j] = adj[i][j]/det; 

	return true; 
} 

void matMatMult(std::vector<std::vector<double> > & A,
		std::vector<std::vector<double> > & B,
		std::vector<std::vector<double> > & C)
{
	int M = A.size();
	int K = B.size();
	int N = B[0].size();
	C.resize(M, std::vector<double>(N,0.0));
	for(unsigned int i = 0; i < M; ++i)
	{
		for(unsigned int j = 0; j < N; ++j)
		{
			for(unsigned int k = 0; k < K; ++k)
				C[i][j] += A[i][k]*B[k][j];
		}
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

void printMat(const std::vector<std::vector<double> > & A)
{

	int M = A.size();
	int N = A[0].size();
	for(unsigned int i = 0; i < M; ++i)
	{
		for(unsigned int j = 0; j < N; ++j)
		{
			std::cout << A[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}

// y = B^{-1}*A*x
void rotateCoordinates(const std::vector<double> & x,
		       const std::vector<std::vector<double> > & A,
		       const std::vector<std::vector<double> > & BInv,
		       std::vector<double> & y)	 
{
	y.resize(BInv.size(),0.0);
	std::vector<double> z(A.size(),0.0);
	matVecMult(A, x, z);
	matVecMult(BInv, z, y);
}
			

int main()
{

	int N = 3;
	std::string fileName1, fileName2;
	std::cout << "Enter filename for the basis vectors for first coordinate system: ";
	std::cin >> fileName1;
	std::cout << "Enter filename for the basis vectors for second coordinate system: ";
	std::cin >> fileName2; 
	std::vector<std::vector<double> > A(N, std::vector<double>(N,0.0));
	std::vector<std::vector<double> > B(N, std::vector<double>(N,0.0));
	for(unsigned int i = 0; i < N; ++i)
	{	
		for(unsigned int j = 0; j < N; ++j)
		{ 
			A[i][j] = (rand()/(RAND_MAX + 0.0));
		}
	}

	inverse(A, AInv);
	readBasisVectors(fileName, A);	
	readBasisVectors(fileName, B);	
} 
