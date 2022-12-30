#include <cstdlib>
#include <iostream>
#include <vector>

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

