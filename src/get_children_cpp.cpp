#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::set<int> get_children(int thisrow, IntegerVector parentrow) {
    // get all the subsequent children internode rows given a particular internode row
    int n = parentrow.size();
    std::set<int> daughterrows; // sets are ordered, so faster to look up values for comparison than vectors
    daughterrows.insert(thisrow);

    for(int i = n-1; i > 0; i--) { // traverse from base to tips.  relies on row order being correct.
        if (daughterrows.find(parentrow[i]) != daughterrows.end()) {
            daughterrows.insert(i+1); // beware 0-indices in c
        }
    }

    return daughterrows;
}

/*** R
get_children(4, c(2,3,4,5,NA))
*/
