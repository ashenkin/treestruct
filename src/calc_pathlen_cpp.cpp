#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector calc_pathlen_cpp(NumericVector len, NumericVector parent_idx){
    // idx is a pointer to the row of the parent of the current internode
    // idx must be ordered from tip to base
    int n = len.size();
    NumericVector pathlen;
    pathlen = clone(len);

    for(int currrow = n - 1; currrow >= 0; currrow--){
        if (!R_IsNA(parent_idx[ currrow ])) {
            pathlen[ currrow ] += pathlen[ parent_idx[ currrow ] - 1 ]; // minus one since c is zero-indexed
        }
    }

    return pathlen;
}


/*** R
calc_pathlen_cpp(1:10, c(2:10,NA))
*/
