#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector assign_branchnum_cpp(NumericVector furcations, LogicalVector is_tip){
    // furcations is the vector of furcations at each internode.  1 = no furcations.
    // is_tip is a logical vector indicating whether that row is a branch tip
    int n = furcations.size();
    NumericVector branchnum(n);
    int currbranch=1;
    branchnum[0]=currbranch;
    for(int i = 1; i < n; i++){
        if (furcations[i] > 1 || is_tip[i]) {
            currbranch++;
        }
        branchnum[i] = currbranch;
    }
    return branchnum;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
assign_branchnum_cpp(rep(c(rep(1, 10), 2), 10), rep(c(T, F, F, F, F), 22))
*/
