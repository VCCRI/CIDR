#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;
using namespace RcppParallel;

// helper inline code to calculate individual dissimilarity entries - given
// particular count and truth "column" vectors
template <typename InputIterator1, typename InputIterator2, typename InputIterator3, typename InputIterator4>
inline double
    my_accumulate(InputIterator1 count_col1_it,
                  InputIterator1 count_col1_end,
                  InputIterator2 count_col2_it,
                  InputIterator3 truth_col1_it,
                  InputIterator4 truth_col2_it,
                  double threshold){
        double part_dist=0;
        for (; count_col1_it != count_col1_end ; ++count_col1_it) {
            if ((!*truth_col1_it & !*truth_col2_it) | (std::max(*count_col1_it,*count_col2_it) > threshold)) {
                part_dist += pow(std::abs(*count_col1_it - *count_col2_it),2);
            }
            ++count_col2_it;
            ++truth_col1_it;
            ++truth_col2_it;
        }
        return part_dist;
    }

// worker object for use in parallel code - calculates dissimilarity matrix entries
// for a given begin/end range
struct CidrDistance : public Worker {
    
    // Distance matrix - input & output
    RMatrix<double> dist;
    
    // truth matrix - input
    const RMatrix<int> truth;
    
    // Counts matrix - input
    const RMatrix<double> counts;
    
    const double threshold;
    
    // cidrDistance constructor - convert input matrices to RMatrix type
    CidrDistance(NumericMatrix dist, const IntegerMatrix truth, const NumericMatrix counts, const double threshold) 
        : dist(dist), truth(truth), counts(counts), threshold(threshold) {}
    
    // TODO - only expecting begin & end  - as parallelFor only works with begin/end 
    // maybe don't need all the others - just counts_begin & counts_end
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i<end; ++i) {
            // hold first columns in memory
            RMatrix<double>::Column count_col_1 = counts.column(i);
            RMatrix<int>::Column truth_col_1 = truth.column(i);
            for (std::size_t j=0; j<i; ++j) {
                RMatrix<double>::Column count_col_2 = counts.column(j);
                RMatrix<int>::Column truth_col_2 = truth.column(j);
                // this needs changing to new namenclature
                dist(j,i) = my_accumulate(count_col_1.begin(), count_col_1.end(), count_col_2.begin(), truth_col_1.begin(), truth_col_2.begin(), threshold);
            }
        }
    }
};

// [[Rcpp::export]]
NumericMatrix cpp_dist(NumericMatrix dist, IntegerMatrix truth, NumericMatrix counts, int ncol, double threshold) {
    // create worker object for parallel code
    CidrDistance cidrDistance(dist, truth, counts, threshold); 
    
    // do the work in parallel
    parallelFor(0, ncol, cidrDistance);
    
    return dist;
}



struct ProbDropoutFunc {
  double a;
  double b;
  
  inline double getProb(double x) {
    return 1 / ( 1 + exp(a*(x - b)) );
  }
  
  inline double distance(double count1, int isDropoutCandidate1, double count2, int isDropoutCandidate2) {
    if (isDropoutCandidate1 && isDropoutCandidate2) {
      return 0;
    } else {
      if (isDropoutCandidate1 /*&& !isDropoutCandidate2*/) {
        double p = getProb(count2);
        count1 = p*count2 + (1 - p)*count1;
      }
      if (isDropoutCandidate2 /*&& !isDropoutCandidate1*/) {
        double p = getProb(count1);
        count2 = p*count1 + (1 - p)*count2;
      }
      double diff = count1 - count2;
      return diff*diff;
    }
  }
};

// helper inline code to calculate individual dissimilarity entries - given
// particular count and truth "column" vectors
template <typename InputIterator1, typename InputIterator2, typename InputIterator3, typename InputIterator4>
inline double
  weighted_accumulate(InputIterator1 count_col1_it,
                      InputIterator1 count_col1_end,
                      InputIterator2 count_col2_it,
                      InputIterator3 truth_col1_it,
                      InputIterator4 truth_col2_it,
                      ProbDropoutFunc func) {
    double part_dist=0;
    for (; count_col1_it != count_col1_end ; ++count_col1_it) {
      part_dist += func.distance(*count_col1_it, *truth_col1_it, *count_col2_it, *truth_col2_it);
      ++count_col2_it;
      ++truth_col1_it;
      ++truth_col2_it;
    }
    return part_dist;
  }

// worker object for use in parallel code - calculates dissimilarity matrix entries
// for a given begin/end range
struct WeightedCidrDistance : public Worker {
  
  // Distance matrix - input & output
  RMatrix<double> dist;
  
  // truth matrix - input
  const RMatrix<int> truth;
  
  // Counts matrix - input
  const RMatrix<double> counts;
  
  const ProbDropoutFunc coefs;
  
  // cidrDistance constructor - convert input matrices to RMatrix type
  WeightedCidrDistance(NumericMatrix dist, const IntegerMatrix truth, const NumericMatrix counts, const ProbDropoutFunc coefs) 
    : dist(dist), truth(truth), counts(counts), coefs(coefs) {}
  
  // TODO - only expecting begin & end  - as parallelFor only works with begin/end 
  // maybe don't need all the others - just counts_begin & counts_end
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i<end; ++i) {
      // hold first columns in memory
      RMatrix<double>::Column count_col_1 = counts.column(i);
      RMatrix<int>::Column truth_col_1 = truth.column(i);
      for (std::size_t j=0; j<i; ++j) {
        RMatrix<double>::Column count_col_2 = counts.column(j);
        RMatrix<int>::Column truth_col_2 = truth.column(j);
        // this needs changing to new namenclature
        dist(j,i) = weighted_accumulate(count_col_1.begin(), count_col_1.end(), count_col_2.begin(), truth_col_1.begin(), truth_col_2.begin(), coefs);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix cpp_dist_weighted(NumericMatrix dist, IntegerMatrix truth, NumericMatrix counts, int ncol, double a, double b) {
  // create worker object for parallel code
  ProbDropoutFunc coefs = { a, b }; 
  WeightedCidrDistance cidrDistance(dist, truth, counts, coefs); 
  
  // do the work in parallel
  parallelFor(0, ncol, cidrDistance);
  
  return dist;
}
