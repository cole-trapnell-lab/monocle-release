#include <Rcpp.h>
using namespace Rcpp;

// Compute jaccard coefficient between nearest-neighbor sets
//
// Weights of both i->j and j->i are recorded if they have intersection. In this case
// w(i->j) should be equal to w(j->i). In some case i->j has weights while j<-i has no
// intersections, only w(i->j) is recorded. This is determinded in code `if(u>0)`. 
// In this way, the undirected graph is symmetrized by halfing the weight 
// in code `weights(r, 2) = u/(2.0*ncol - u)/2`.
//
// Author: Chen Hao, Date: 25/09/2015; updated by Xiaojie Qiu Nov. 12, 2017


// [[Rcpp::export]]
NumericMatrix jaccard_coeff(NumericMatrix idx, bool weight) {
  int nrow = idx.nrow(), ncol = idx.ncol(), r = 0;
  NumericMatrix weights(nrow*ncol, 3);
  
  for(int i = 0; i < nrow; i ++) {
    for(int j = 0; j < ncol; j ++) {
      int k = idx(i,j) - 1;
      
      weights(r, 0) = i + 1;
      weights(r, 1) = k + 1;
      weights(r, 2) = 1;
      
      if(weight == TRUE) {
        
        NumericVector nodei = idx(i, _);
        NumericVector nodej = idx(k, _);
        
        int u = intersect(nodei, nodej).size();  // count intersection number
        int v = 2 * ncol - u;  // count union number
        
        if(u>0) { 
          // weights(r, 0) = i + 1;
          // weights(r, 1) = k + 1;
          // weights(r, 2) = u / (2.0 * ncol - u) / 2;  // symmetrize the graph
          
          weights(r, 2) = (double) u / (double) v;  // normalize the values
        }
      }
      
      r ++;
      
    }
  }
  
  weights(_, 2) = weights(_, 2) / max(weights(_, 2));
  
  return weights;
}

/***
 edges$C = jaccard_dist
 edges = subset(edges, C != 0)
 edges$C = edges$C/max(edges$C)
 */