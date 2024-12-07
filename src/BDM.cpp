#include <Rcpp.h>
using namespace Rcpp;

//' @title BDM
//' @description this function computes the backward degree metric of a graph
//' @param graph A numeric matrix representing the graph adjacency matrix
//' @return A numeric vector containing the backward degree metric for each node
//' @export
// [[Rcpp::export]]
NumericVector BDM(NumericMatrix graph) {
   int m = graph.ncol();  // number of columns (nodes in the graph)
   NumericVector backward_degree(m, 0.0);  // Initialize vector to store backward degree

   // Loop through each node
   for (int i = 0; i < m; ++i) {
     // Calculate the backward degree for each node by summing the adjacency values
     for (int j = 0; j <= i; ++j) {  // Iterating over the upper triangular part (or full matrix)
       backward_degree[i] += graph(i, j);
     }
   }
   return backward_degree;
}
