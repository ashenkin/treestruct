#include <Rcpp.h>
using namespace Rcpp;

// from https://stackoverflow.com/questions/46368188/how-to-efficiently-calculate-path-lengths-in-tree-topology/46369457#46369457

// [[Rcpp::export]]
NumericVector calc_pathlen_cpp(NumericVector len, NumericVector idx){
    int n = len.size();
    NumericVector res(n);
    double cumsum=0;
    int j;
    res[0] = len[0];

    for(int i = 1; i < n; i++){
        j = idx[ i ] - 1;
        cumsum = res[ j ] + len[ i ];
        res[ i ] = cumsum;
    }
    return res;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
calc_pathlen_cpp(runif(100), c(1:100))
*/
