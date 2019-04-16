#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
List assign_branchnum_cpp(NumericVector furcations, LogicalVector is_tip){
    // furcations is the vector of furcations at each internode.  1 = no furcations.
    // is_tip is a logical vector indicating whether that row is a branch tip
    // this routine relies on internodes being in the correct order (which reorder_internodes does).
    int n = furcations.size();
    NumericVector branchnum(n);
    NumericVector order_in_branch(n);
    int currbranch=1;
    int currorderinbranch=1;
    branchnum[0]=1;
    order_in_branch[0]=1;
    for(int i = 1; i < n; i++){
        if (furcations[i] > 1 || is_tip[i]) {
            currbranch++;
            currorderinbranch = 1;
        }
        branchnum[i] = currbranch;
        order_in_branch[i] = currorderinbranch++;
    }
    List ret;
    ret["branchnum"] = branchnum;
    ret["order_in_branch"] = order_in_branch;
    return ret;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
assign_branchnum_cpp(rep(c(rep(1, 10), 2), 10), rep(c(T, F, F, F, F), 22))
*/
