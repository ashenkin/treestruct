#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calc_total_x_above_internode_cpp(NumericVector x, NumericVector parentrow) {
    // this is a generic algorithm to compute the total [anything (e.g. vol, sa)] above each internode.
    int n = x.size();
    NumericVector cum;
    cum = clone(x);

    // accumulate all internodes above downwards from tips.  relies on row order being correct.
    //for(int i = n - 1; i > 0; i--) {
    for(int i = 0; i < n - 1; i++) {
        cum[ (parentrow[i] - 1) ] += cum[ i ];  // -1 since C is zero-based
    }

    return cum;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
calc_total_x_above_internode_cpp(1:10, c(2:10,0))
*/
