#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calc_sa_above_cpp(NumericVector sa, NumericVector parentrow) {
    int n = sa.size();
    NumericVector cum;
    cum = clone(sa);

    // accumulate all internodes above downwards from tips.  relies on row order being correct.
    for(int i = n - 1; i > 0; i--) {
        cum[ (parentrow[i] - 1) ] += cum[ i ];  // -1 since C is zero-based
    }

    return cum;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
calc_sa_above_cpp(1:10, 0:9)
*/
