#include <Rcpp.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;
inline double
my_accumulate(NumericMatrix::Column::iterator __ci_first,
              NumericMatrix::Column::iterator __ci_last,
              NumericMatrix::Column::iterator __cj_first,
              IntegerMatrix::Column::iterator __ti_first,
              IntegerMatrix::Column::iterator __ti_last,
              IntegerMatrix::Column::iterator __tj_first,
              double threshold){
    double __part_dist=0;
    for (; __ci_first != __ci_last; ++__ci_first) {
        if ((!*__ti_first & !*__tj_first) | (std::max(*__ci_first,*__cj_first) > threshold)) {
                __part_dist += pow(std::abs(*__ci_first - *__cj_first),2);
            }
            ++__cj_first;
            ++__ti_first;
            ++__tj_first;
    }
    return __part_dist;
}

// [[Rcpp::export]]
NumericMatrix cpp_dist(NumericMatrix D, IntegerMatrix T, NumericMatrix C, int ncol, double threshold) {
    int i, j;
    i = j = 0;
    for (i=0; i < ncol; i++) {
        NumericMatrix::Column c1 = C(_, i);
        IntegerMatrix::Column t1 = T(_, i);
        for (j=0; j<i; j++) {
            NumericMatrix::Column c2 = C(_, j);
            IntegerMatrix::Column t2 = T(_, j);
            D(j,i) = my_accumulate(c1.begin(), c1.end(), c2.begin(), t1.begin(), t1.end(), t2.begin(), threshold);
        }
    }
    return D;
}
