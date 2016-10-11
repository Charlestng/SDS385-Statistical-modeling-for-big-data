#include <RcppEigen.h>
#include <algorithm>    // std::max

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseVector;
typedef Eigen::MappedSparseMatrix<double>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
typedef Map<VectorXi>  MapVeci;


template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

inline double invSqrt( const double& x ) {
    double y = x;
    double xhalf = ( double )0.5 * y;
    long long i = *( long long* )( &y );
    i = 0x5fe6ec85e7de30daLL - ( i >> 1 );//LL suffix for (long long) type for GCC
    y = *( double* )( &i );
    y = y * ( ( double )1.5 - xhalf * y * y );
    
    return y;
}


list sparsesgd_logit(MapMatd X, VectorXd Y, VectorXd M, double eta, int npass, VectorXd beta0, double lambda=1.0, double discount = 0.01){
	
	// Initialise the number of observations and the number of features
	int n_obs = X.cols();
	int n_features = X.rows();
	
	// Variables we will need to keep values during runtime
	SparseVector<double> x(n_features);
		
	// LOOP 1: General loop that counts the number of times we go over the whole data set
	for(int pass = 0, pass < npass, pass++){
		
		// LOOP 2: going through each observation in the data set for stochastic gradient descent
		for(int obs, obs < n_obs, obs++){
			
			// The goal here is to update the intercept "alpha" at each iteration based on:
			// The value of Y_hat
			// The value of beta
			
			x = X.innerVector(obs);
			
		}
	}
}


