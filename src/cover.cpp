#include <Rcpp.h>
using namespace Rcpp;

// Given a logical vector, returns an integer vector (1-based) of the positions in the vector which are true
IntegerVector which_true(LogicalVector x) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; i++) { if (x[i]) y.push_back(i+1); }
  return wrap(y);
}

// Given an (n x 2d) matrix of min/max bounds of any 'iso-oriented' rectangles in the plane, i.e. whose edges are parallel to the coordinate axes, 
// and a matrix of point cloud data 'x', return a list of indices of points in x which fall in the level set bounds given by 'bnds'.
// [[Rcpp::export]]
List constructIsoAlignedLevelSets(const NumericMatrix& x, const NumericMatrix& bnds, bool save_bounds=true){
  if (x.ncol() != (bnds.ncol()/2)){ Rcpp::stop("dimension of points != dimension of bounds matrix / 2."); }
  const int n_level_sets = bnds.nrow(), d = bnds.ncol()/2;
  List level_sets = List(n_level_sets);
  LogicalVector level_set_test = LogicalVector(x.nrow(), true); // which pts lie in the set; use logical vector to shorten code and use vectorized &
  for (int i = 0; i < n_level_sets; ++i){
    NumericMatrix::ConstRow ls_bnds = bnds.row(i); // Update level set bounds
    std::fill(level_set_test.begin(), level_set_test.end(), true);// Reset to all true prior to doing logical range checks
    for (int d_i = 0; d_i < d; ++d_i){
      level_set_test = level_set_test & ((x.column(d_i) >= ls_bnds[d_i]) & (x.column(d_i) <= ls_bnds[d + d_i]));
    }

    // Save the level set
    IntegerVector ls = which_true(level_set_test);
    if (save_bounds){ ls.attr("bounds") = ls_bnds; }
    level_sets[i] = ls;
  }
  return(level_sets);
}
