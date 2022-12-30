#if !defined(dft_GaussianFunctionManagerMO_h)
#define fem_GaussianFunctionManagerMO_h

#include <string>
#include <vector>
#include <set>
//
//
//

namespace dft {

  //
  // forward declarations
  //

  /**
   * @brief Singleton class to read and store the Gaussian basis parameters to construct the input density
   * 
   */
  class GaussianFunctionManagerMO {

    //
    // types
    //
    public:

      struct contractedGaussian 
      {
	int L;
	double * alpha;
	double * c;
	char lquantum;
      };

      struct basis {
	int L; // Number of Gaussians used for contraction
	const double * alpha; // exponents of the L Gaussians used for contraction
	const double * c; // Coefficients of for each of the 
	int l; // l quantum number
	int m; // m quantum number
	double * normConsts; //Normalization constant for each primitive Gaussian
	const double * origin; // Origin or position of the atom about which the Gaussian is defined
	double basisNormConst; // Normalization constant for the resulting contracted Gaussian
      };

      //
      // methods
      //
    public:

      /**
       * @brief Constructor 
       */
      GaussianFunctionManagerMO();

      /**
       * @brief Destructor 
       *
       * @param enrichedFunctionId Id of the enrichedFunction
       * @param mollifierRadius The mollifier radius value to set it to
       */
      ~GaussianFunctionManagerMO();

      double getMOValue(const double * point, const int MOId);

      int getNumberMOs();
      double getBasisFunctionValue(const int basisId, const double * point);

      /**
       * @brief get the gradient of a basis function at a given point
       *
       * @param basisId Id of the basis function
       * @param point Point at which the density is to be computed
       *
       * @return vector containing the gradient of the basis function at the point
       */
      std::vector<double>
	getBasisFunctionGradient(const int basisId, const double * point);

      /**
       * @brief get the double derivatives of the basis function at a point
       *
       * @param point Point at which the basis function is to be computed
       *
       * @return double derivatives of the basis function at the point
       */
      std::vector<double>
	getBasisFunctionDoubleDerivatives(const int basisId, const double * point);

      /**
       * @brief get the laplacian of a basis function at a given point
       *
       * @param basisId Id of the basis function
       * @param point Point at which the density is to be computed
       *
       * @return laplacian of the basis function at the point
       */
      double
	getBasisFunctionLaplacian(const int basisId, const double * point);

      /**
       * @brief get the number of basis functions
       *
       * @return Number of basis functions
       */
      int getNumberBasisFunctions();

    private:

      //
      // store atomic coordinates
      // 
      std::vector<std::vector<double> > d_atomicCoords;

      //
      // store the MOs 
      // 
      std::vector<std::vector<double> > d_MOMat;

      //
      // store basis file names for each atom
      //
      std::vector<std::string> d_basisFileNames;

      //
      //store the unique basis file names
      //
      std::set<std::string> d_uniqueBasisFileNames;

      //
      //store the contracted gaussians for each atom type
      //
      std::vector<std::vector<dft::GaussianFunctionManagerMO::contractedGaussian*> > d_contractedGaussians;

      //
      // store basis function paramters
      //
      std::vector<dft::GaussianFunctionManagerMO::basis*> d_basisFunctions;


  };

}

#endif // dft_GaussianFunctionManagerMO_h

