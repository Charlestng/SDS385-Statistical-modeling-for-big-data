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


List sparsesgd_logit(MapMatd X, VectorXd Y, VectorXd M, int npass, VectorXd beta0, double lambda=1.0, double discount = 0.01){
	// X is the design matrix stored in column-major format
  	// i.e. with features for case i stores in column i
  	// Y is the vector of counts
	// M is the vector of sample sizes
	  // Thus Y[i] ~ Binomial( M[i], w[i])  )
	  // w[i] = 1/{1+exp(- x[i] dot beta)}
  	// where Beta is the regression vector we want to estimate
  	// lambda is the l1 regularization parameter
  
  
	// Initialise the number of observations and the number of features
	int n_obs = X.cols();
	int n_features = X.rows();
	
	// Initialise parameters
		// w_hat is intialised this way because we want to avoid it being equal to 1.0 or 0.0
		// It is the expectation of Y
	double w_hat = (Y.sum() + 1.0) / (M.sum() * 2.0);
		// alpha is the inverse logistic of the expectation of Y
	double alpha = log(w_hat / (1.0 - w_hat));
	
	double beta = beta0;
	
	// Variables we will need in order to work with observations during runtime
	SparseVector<double> x(n_features);
	
	// Bookkeeping variables for the log likelihood and the values of beta
	NumericVector n_LL_track(npass * n_obs, 0.0);
	double n_LL_avg = 0.0;
	double e_psi;
	double psi;
	int global_iterator;
	int feat;
	double y_hat;
	double delta;
	double g0squared;	
	double weight;
	double gammatilde;
	double mu;
	double grad;
	
	// Initialize Gsquared and beta
	VectorXd beta(numfeatures);
	VectorXd Gsquared(numfeatures);
	for(int feat=0; feat < numfeatures; feat++) {
		Gsquared(feat) = 1e-3;
    	beta(feat) = beta0(feat);
		}	
	
	// Keeping track of the last time a feature's parameter has been updated for the purpose
	// of lazy update
	NumericVector last_update(n_features, 0.0);
		
	// LOOP 1: General loop that counts the number of times we go over the whole data set
	for(int pass = 0, pass < npass, pass++){
		// Initialising the number of iterations 
		global_iterator = 0;
		
		// LOOP 2: going through each observation in the data set for stochastic gradient descent
		for(int obs, obs < n_obs, obs++){
			
			// The goal here is to update the intercept "alpha" at each iteration based on:
			// The value of Y_hat
			// The value of beta
			
			// First we calculate the estimate of Y at this time of the loop
			
			x = X.innerVector(obs);
			psi = alpha + x.dot(beta);
			e_psi = exp(psi);
			y_hat = M[obs] * e_psi / (1.0 + e_psi);
			
			// Then we calculate the negative log likelihood using y_hat
				// Using a decaying moving average of the contribution to the n_LL of each observation
			n_LL_avg = n_LL_avg * (1.0 - discount) + (M[obs] * log(1.0 + e_psi) + (M[obs]- Y[obs]) * psi);
			n_LL_track[global_iterator] = n_LL_avg;
			
			// Update the intercept
			// Update intercept 
			delta = Y[obs] - y_hat;
			g0squared += delta*delta;
			alpha += (1.0 / sqrt(g0squared)) * delta;
			
			// LOOP 3: going through each feature for each observation, except the sparsity of the feature vector
			// for each observation allows us to do lazy updating thanks to the "last_update" tracker
			for (SparseVector<double>::InnerIterator it(x); it; ++it) {
				
				// which feature are we working with ?
				feat = it.index();
					
				// For a better fit, it's useful to put a lower penalty on the high values of beta
        		weight = 1.0 / (1.0 + fabs(beta(feat))); 
        		
        		// Step 1: Since the features are sparse, we don't update them all at each iteration
        		// there fore this step aggregates the penalty for all the times we did not touch this feature
        		// using the general iteration counter global_iterator and the last update tracker
        		double skip = global_iterator - last_update(feat);
        		gammatilde = skip / sqrt(Gsquared(feat));
        		beta(feat) = sgn(beta(feat)) * fmax(0.0, fabs(beta(feat)) - gammatilde * weight * lambda);
        		
        		// Update the last-update vector
        		last_update(feat) = global_iterator;
        		
        		// Step 2: Now we compute the update for this observation and this feature.

        		// gradient of negative log likelihood
        		grad = -delta * it.value();
        		
        		// update adaGrad scaling for this feature
        		Gsquared(feat) += grad * grad;
        		
        		// scaled stepsize
				gammatilde = 1.0 / sqrt(Gsquared(feat));       			
        		mu = beta(feat) - gammatilde * grad;
        		beta(feat) = sgn(mu) * fmax(0.0, fabs(mu) - gammatilde * weight * lambda);
			}
			global_iterator++; // increment global counter
		}
	}
	// At the very end, apply the accumulated penalty for the variables we haven't touched recently
	for(int feat=0; feat < n_features; feat++) {
		double skip = global_iterator - last_update(feat);
        gammatilde = skip / sqrt(Gsquared(feat));
        beta(feat) = sgn(beta(feat)) * fmax(0.0, fabs(beta(feat)) - gammatilde * weight * lambda);
	}
	return List::create(Named("alpha") = alpha,
                      Named("beta") = beta,
                      Named("n_LL_track") = n_LL_track);
}
