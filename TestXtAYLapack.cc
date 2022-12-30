#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <complex>
#include <valarray>
#include <cassert>
#include <time.h>

extern "C" {
	
  void dgemm_(char* transA,
                    char* transB,
		    int *m,
		    int *n,
		    int *k,
		    double *alpha,
		    double *A,
		    int *lda,
		    double *B,
		    int *ldb,
		    double *beta,
		    double *C,
		    int *ldc);
      
  void zgemm_(char* transA,
                    char* transB,
                    int *m,
                    int *n,
                    int *k,
		    std::complex<double> * alpha,
		    std::complex<double> * A,
		    int *lda,
		    std::complex<double> * B,
		    int *ldb,
		    std::complex<double> * beta,
		    std::complex<double> * C,
		    int *ldc);


  void dsyevd_(char* jobz, 
  		     char* uplo, 
  		     int* n, 
  		     double* A, 
  		     int *lda, 
  		     double* w, 
                     double* work, 
  		     int* lwork, 
  		     int* iwork, 
  		     int* liwork, 
  		     int* info);

  void dcopy_(int *n,
	      double *x,
	      int *incx,
	      double *y,
	      int *incy);

}

struct gen_rand { 
    double range;
  public:
    gen_rand(double r=1.0) : range(r) {}
    double operator()() 
    { 
      return (rand()/(double)RAND_MAX) * range;
    }
};
void initVectorRandom(std::vector<double> & x)
{
  const int size = x.size();
  std::generate_n(x.begin(), size, gen_rand());
  
}

void symmetrize(std::vector<double> & X, 
		const int size)
{

  for(unsigned int i = 0; i < size; ++i)
  {
    for(unsigned int j = i; j < size; ++j)
    {
      double & Xij = X[i*size+j];
      double & Xji = X[j*size+i];
      const double val = 0.5*(Xij+Xji);
      Xij = val;
      Xji = val;
    }
  }
}


int main()
{

  int Asize, Xsize, Ysize;
  std::cout << "Enter A size: ";
  std::cin >> Asize;
  std::cout << std::endl;
  std::cout << "Enter X size: ";
  std::cin >> Xsize;
  std::cout << std::endl;
  std::cout << "Enter Y size: ";
  std::cin >> Ysize;
  std::cout << std::endl;

  std::vector<double> A(Asize*Asize);
  std::vector<double> X(Asize*Xsize);
  std::vector<double> Y(Asize*Ysize);
  std::vector<double> AY(Asize*Ysize,0.0);
  std::vector<double> Z(Xsize*Ysize,0.0);
  std::vector<double> ZLP(Xsize*Ysize,0.0);
  

  initVectorRandom(A);
  initVectorRandom(X);
  initVectorRandom(Y);
  symmetrize(A,Asize);

  std::cout << "Init done" << std::endl;
  std::cout << "Xsize: " <<  Xsize << std::endl;
  std::cout << "Ysize: " <<  Ysize << std::endl;
  std::cout << "Asize: " <<  Asize << std::endl;

  char transA = 'N';
  char transB = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int inc = 1;
  int Ysize1 = Ysize;

  dgemm_(&transA,
	   &transB,
	   &Asize,
	   &Ysize,
	   &Asize,
	   &alpha,
	   &A[0],
	   &Asize,
	   &Y[0],
	   &Asize,
	   &beta,
	   &AY[0],
	   &Asize);


}


