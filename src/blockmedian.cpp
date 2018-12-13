// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List blockmedian(NumericMatrix Xr, int w) {
    int n = Xr.nrow(), p = Xr.ncol();
    arma::mat X(Xr.begin(), n, p, false);
    arma::mat Xnew(n, p, arma::fill::zeros);
    arma::vec temp_v(w * w);
    arma::vec v = arma::vectorise(X);
    int rmin = 0, rmax = 0, cmin = 0, cmax = 0;
    double med;
    
    // median of all non-NA values in matrix
    med = arma::median(v.elem(arma::find_finite(v)));

    // median values in w x w windows
    for (int i = 0; i < n; ++i) {
        rmin = std::max(0, (int) (i - floor(w / 2.0)));
        rmax = std::min(rmin + w - 1, n - 1);
        rmin = std::max(rmax - w + 1, 0);
        for (int j = 0; j < p; ++j) {
            cmin = std::max(0, (int) (j - floor(w / 2.0)));
            cmax = std::min(cmin + w - 1, p - 1);
            cmin = std::max(cmax - w + 1, 0);
            temp_v = arma::vectorise(X.submat(rmin, cmin, rmax, cmax));
            temp_v = temp_v.elem(arma::find_finite(temp_v));
            if (temp_v.n_elem > 0) {
                try {
                    Xnew(i, j) = arma::median(temp_v);
                } catch (...) {
                    cout << "Error at i = " << i << ", j = " << j;
                    Xnew(i, j) = arma::datum::nan;
                };
            } else {
                Xnew(i, j) = arma::datum::nan;
            };
        };
    };
    
    return Rcpp::List::create(Named("local") = Xnew,
                              Named("global") = med);
}
