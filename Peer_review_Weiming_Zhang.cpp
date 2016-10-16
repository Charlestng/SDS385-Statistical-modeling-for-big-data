// Peer review:
// Reviewer: Charles TANGUY, Reviewee: Weiming Zhang
// All the entries of the peer review will start by the sequence:
//
//review
//
// For more readability

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

typedef MapMatd::InnerIterator InIterMat;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// Use this function to calculate inverse square root.
inline double invSqrt( const double& x ) {
double y = x;
double xhalf = ( double )0.5 * y;
long long i = *( long long* )( &y );
i = 0x5fe6ec85e7de30daLL - ( i >> 1 ); //LL suffix for (long long) type for GCC
y = *( double* )( &i );
y = y * ( ( double )1.5 - xhalf * y * y );

return y;
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

List improvingSGD(MapMatd X, VectorXd Y, VectorXd m, double step, VectorXd beta0, double lambda=0.0, int npass=1){
  // X is the design matrix stored in column-major format
  // i.e. with features for case i stores in column i
  // Y is the vector of counts
  // M is the vector of sample sizes
  // Thus Y[i] ~ Binomial( M[i], w[i])  )
  // w[i] = 1/{1+exp(- x[i] dot beta)}
  // where Beta is the regression vector we want to estimate
  // lambda is the l1 regularization parameter

  //INITIALIZE VALUES:
  int numobs = X.cols(); 
  int numfearures = X.rows();

  int iter = 0;
  int j = 0;
  double epsi = 1e-4;
  double Yi = Y(0);
  double mi = m(0);
  double grad = 0;
  double wi = .5;
  double wi_exponent = 0;

  VectorXd beta(numfearures);			//Initialize Beta_hat vector
  VectorXd hist_grad(numfearures);		//Initialize vector to store running hist_grad
  VectorXd adj_grad(numfearures);		//Initialize vector to store adj_grad
  for (int i=0;i<numfearures;i++){
    hist_grad(i) = 1e-3;
    beta(i)=beta0(i);
    adj_grad(i) = 0;
  }

  //Initialize elements to hold X, Y and m for a single observation (column).
  SparseVector<double> Xi(numobs);
  Xi=X.innerVector(0);

  // Bookkeeping: how long has it been since the last update of each feature?
  NumericVector last_updated(numfearures,0.0);
  
  //Initialize vectors to hold log-likelihood and running avg neg log-likelihood.
  double neglik = 0;
  NumericVector neglik_ra(numobs*npass,0.0);
  
  //Initialize variable to calculate penalty.
  double skip = 1;
  double accum_l2_penalty = 0;
  double gam = 0;
  
  // Outer loop: number of passes over data set
  double k = 0; // global interation counter
  for (int pass=0; pass < npass; pass++){

    //Loop over each observation (columns of x)
    for (int i=0; i<numobs; i++){

      //For linear predictor and E(Y[i]) from features
      Xi = X.innerVector(i);
      Yi = Y(i);
      mi = m(i);			

      //Update wi probability value
      wi_exponent = Xi.dot(beta);
      wi = 1 / (1 + exp(-wi_exponent));

      //Update neglik
      neglik = -(Yi * log(wi + epsi) + (mi - Yi) * log(1 - wi + epsi));
      if(iter > 0) {
        neglik_ra(iter) = (neglik_ra(iter-1) * (iter-1) + neglik) / iter;
      }	
	  
//
// Review
//
// I am a little confused with line 117, you run the statement if and only if the integer variable iter
// is strictly bigger than 0, the initial value being 0, you don't update the negative log likelihood
// before the first update of beta. What confuses me is that I could not find a place in the code where
// the value of iter is updated. This is a problem because the variable iter helps you update the L2
// penalty at line 160. My intuition is that maybe iter plays an equivalent role as k (the global iteration
// counter), which will very easily fix the problem.

      //Iterate over the active features for this instance
      for (SparseVector<double>::InnerIterator it(Xi); it; ++it){
        j = it.index();
        skip = k - last_updated(j);
        if (skip > 5){  skip = 5;} //Cap skip at 5
        last_updated(j) = k;

        //Update penalty and beta
        gam = step * adj_grad(j);
        accum_l2_penalty = beta(j) * ((1 - pow(1 + lambda * gam,skip)) / (1 - lambda * gam) );
        beta(j) -= accum_l2_penalty; 

//
// Review
//
// Here at line 140, 141, there might be cases where the penalty ((1 - pow(1 + lambda * gam,skip)) / (1 - lambda * gam) )
// will be greater than 1, in that case, the beta update at line 141 will cause beta(j) to change sign.
// The penalty is only supposed to shrink the non important features to zero.
// Thankfully to things in your code compensate for this:
//     - the fact that you limit the number of skips to 5
//     - the fact that the penalty is always of beta's opposite sign
// Worst case scenario, the penalty will make beta be alternatively negative and positive every 5 iteration. There is
// little chance that it will take crazy values.

        //Calculate l2 norm penalty
        double l2penalty = 2*lambda*beta(j);

        //Update hist_grad and gradient
        grad = (mi * wi - Yi) * it.value() + l2penalty;  
        hist_grad(j) += grad * grad;

        //Calculate adj_grad and beta
        adj_grad(j) = grad * invSqrt(hist_grad(j) + epsi);
//
// Review 
//
// Adding an epsilon to the historical sum of square gradients was very smart here at line 163, I did
// not think about that at the beginning :)
        beta(j) -= step*adj_grad(j);	
      }
      k++; //increment global counter
    }	
  }
  
  // At the very end, apply the accumulated penalty for the variables we haven't touched recently
  for (int j=0; j< numfearures; ++j){
    skip = (iter - 1) - last_updated(j); 
    if (skip > 5){  skip = 5;} //Cap skip at 5
    gam = step * adj_grad(j);	
    accum_l2_penalty = beta(j) * ((1 - pow(1 + lambda * gam,skip)) / (1 - lambda * gam) );
    beta(j) -= accum_l2_penalty; 
  }	

  return List::create(Named("beta") = beta,
                      Named("numobs") = numobs,
                      Named("numfearures") = numfearures,
                      Named("iter") = iter,
                      Named("loglik") = neglik_ra);
  
}

//
// Review
//
// Overall, you did a very good job of making your code very readable and accessible for an
// outsider to look at it, and that was one of the main point of having these peer reviews!
// You might want to look into the "iter" variable issue I mentioned, and maybe add a 
// soft_threshold to counter the cumulative penalty if it gets too big, I am not sure if this
// would have a major impact on the speed.
// I do not know if you use a data set for which you add a column of 1s or if you scale your
// data beforehand, but you might want to add a few lines of code to update the intercept 
// independently of beta so it is not penalized. The fit of the model should be better overall.
// Great job Weiming!
