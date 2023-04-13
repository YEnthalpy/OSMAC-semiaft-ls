#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

NumericMatrix center(NumericMatrix M, NumericVector w, int N)
{
  int p = M.ncol();
  int n = M.nrow();
  NumericVector cmean(p);
  for (int i = 0; i < n; i++)
  {
    double wi = w[i];
    for (int j = 0; j < p; j++)
    {
      cmean[j] += (M(i, j) / wi);
    }
  }
  cmean = cmean / n / N;
  NumericMatrix out(n, p);
  for (int i = 0; i < n; i++)
  {
    out(i, _) = M(i, _) - cmean;
  }
  return out;
}

// [[Rcpp::export]]

List km(NumericVector e, NumericVector d,
         NumericVector p)
{
  int t = e.length();
  // int l = 0;
  // while (l < t)
  // {
  //   if (d[l] == 1)
  //   {
  //     break;
  //   }
  //   l++;
  // }
  NumericVector s(t, 1.0);
  // NumericVector eu(t);
  double den = 0;
  int i = t - 1;
  // int k = t - 1;
  while (i >= 0)
  {
    double ei = e[i];
    int j = i;
    double nut = 0;
    while (j >= 0)
    {
      if (e[j] < ei)
      {
        break;
      }
      else
      {
        double pj = p[j];
        den += 1 / pj;
        if (d[j] == 1)
        {
          nut += 1 / pj;
        }
      }
      j--;
    }
    // i = j;
    // if (nut == 0)
    // {
    //   continue;
    // }
    double tmp = 1 - nut / den;
    s[j + 1] = tmp;
    i = j;
    // s[k] = tmp;
    // eu[k] = e[i + 1];
    // k--;
  }
  // s.erase(0, k + 1);
  // eu.erase(0, k + 1);
  // NumericVector out(t - k - 1);
  NumericVector out(t);
  double start = 1;
  for (int i = 0; i < t; i++)
  {
    start = start * s[i];
    out[i] = start;
  }
  // if (l != 0)
  // {
  //   out.push_front(1);
  //   eu.push_front(e[l - 1]);
  // }
  List out1 = List::create(e, out);
  return out1;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_smth(arma::mat x, arma::vec y, arma::uvec d,
                     arma::vec beta, arma::vec p, int n)
{
    arma::vec e = y - x * beta;
    int r = y.n_elem;
    int m = x.n_cols;
    arma::mat out(r, m);
    for (int i = 0; i < r; i++)
    {
        if (d[i] == 0)
        {
            continue;
        }
        arma::rowvec nut(m);
        for (int j = 0; j < r; j++)
        {
            arma::rowvec dx = x.row(i) - x.row(j);
            double rij = arma::norm(dx);
            if (rij == 0)
            {
                continue;
            }
            double tmp = arma::normcdf((e[j] - e[i]) / rij);
            nut += dx * tmp / p[j];
        }
        out.row(i) = nut / p[i] / r / r / n;
    }
    return out;
}


// NumericVector intpol(NumericVector x, NumericVector y1,
//                     NumericVector y2, NumericVector z)
// {
//   int n = x.length();
//   int nout = z.length();
//   NumericVector out(nout);
//   double tmp = x[0] * (1 - y2[0]);
//   for (int i = 0; i < nout; i++)
//   {
//     double v = z[i];
//     int k = 0;
//     int l = n - 1;
//     if (v < x[k])
//     {
//       out[i] = y1[0] + tmp;
//     }
//     else if (v > x[l])
//     {
//       out[i] = 0;
//     }
//     else
//     {
//       while (k < l - 1)
//       {
//         int kl = (k + l) / 2;
//         if (v < x[kl])
//         {
//           l = kl;
//         }
//         else
//         {
//           k = kl;
//         }
//       }
//       if (v == x[l])
//       {
//         out[i] = y1[l] / y2[l];
//       }
//       else
//       {
//         out[i] = y1[k] / y2[k];
//       }
//     }
//   }
//   return out;
// }

// NumericVector multi(NumericMatrix M, NumericVector v)
// {
//   NumericVector out(M.nrow());
//   for (int j = 0; j < M.ncol(); j++)
//   {
//     double tmp = v[j];
//     for (int i = 0; i < M.nrow(); i++)
//     {
//       out[i] += M(i, j) * tmp;
//     }
//   }
//   return out;
// }

// NumericMatrix multi2(NumericMatrix M, NumericVector v)
// {
//   int p = M.ncol();
//   NumericMatrix out(M.nrow(), p);
//   for (int i = 0; i < p; i++)
//   {
//     out(_, i) = M(_, i) * v;
//   }
//   return out;
// }

// NumericMatrix multi3(NumericMatrix M, NumericVector v, int N)
// {
//   int p = M.ncol();
//   NumericMatrix out(M.nrow(), p);
//   for (int i = 0; i < p; i++)
//   {
//     out(_, i) = M(_, i) / v / N;
//   }
//   return out;
// }

// List eRes(NumericVector e, NumericVector d, NumericVector pi,
//           IntegerVector ind_km)
// {
//   Function od("order");
//   Function approx("approx");
//   NumericVector e_sub = e[ind_km - 1];
//   IntegerVector ord_sub = od(e_sub);
//   NumericVector ei_sub = e_sub[ord_sub - 1];
//   NumericVector d_sub = (as<NumericVector>(d[ind_km - 1]))[ord_sub - 1];
//   NumericVector pi_sub = (as<NumericVector>(pi[ind_km - 1]))[ord_sub - 1];
//   NumericVector tmp = km(ei_sub, d_sub, pi_sub);
//   IntegerVector ord = od(e);
//   NumericVector ei = e[ord - 1];
//   NumericVector ei_su = sort_unique(ei_sub);
//   NumericVector s = as<List>(approx(Named("x") = ei_su, Named("y") = tmp,
//                                     Named("xout") = ei, Named("yleft") = 1,
//                                     Named("yright") = 0))["y"];
//   NumericVector ediff = diff(ei);
//   ediff.push_back(0);
//   NumericVector ehat_tmp = cumsum(rev(ediff * s));
//   NumericVector ehat = rev(ehat_tmp);
//   NumericVector ehat2 = ehat + s * ei;
//   NumericVector ehatnew = ehat / s + ei;
//   int p = ehatnew.length();
//   for (int i = 0; i < p; i++)
//   {
//     if (ISNA(ehatnew[i]))
//     {
//       ehatnew[i] = ei[i];
//     }
//   }
//   ehatnew[ord - 1] = clone(ehatnew);
//   List out = List::create(ehatnew, ehat2);
//   return out;
// }

// NumericMatrix uls(NumericMatrix x, NumericVector y, NumericVector delta,
//                   NumericVector pi, NumericVector beta,
//                   IntegerVector ind_km, int N)
// {
//   int r = x.nrow();
//   int p = x.ncol();
//   NumericMatrix x_c = clone(x);
//   x_c.erase(0, r);
//   x_c.attr("dim") = Dimension(r, p - 1);
//   NumericMatrix xdif = center(x_c, pi, N);
//   NumericVector py = multi(x, beta);
//   NumericVector tmp = eRes(y - py, delta, pi, ind_km)[0];
//   NumericVector hy = delta * y + (1 - delta) * (tmp + py);
//   NumericVector ydif = hy - mean(hy / pi) / N;
//   NumericVector beta_c = clone(beta);
//   beta_c.erase(0);
//   NumericMatrix out = multi3(multi2(xdif, ydif - multi(xdif, beta_c)), pi, N);
//   return out;
// }

// List semi_ls_est(NumericMatrix x, NumericVector y,
//                  NumericVector delta, NumericVector pi, int N)
// {
//   Function solve("solve");
//   Function crossprod("crossprod");
//   Function colSums("colSums");
//   int r = x.nrow();
//   int p = x.ncol();
//   NumericVector beta(p);
//   NumericMatrix x_c = clone(x);
//   x_c.erase(0, r);
//   x_c.attr("dim") = Dimension(r, p - 1);
//   NumericMatrix xdif = center(x_c, pi, N);
//   IntegerVector ind_km = seq_len(r);
//   for (int i = 0; i < 100; i++)
//   {
//     NumericVector py = multi(x, beta);
//     NumericVector tmp = eRes(y - py, delta, pi, ind_km)[0];
//     NumericVector hy = delta * y + (1 - delta) * (tmp + py);
//     NumericVector ydif = hy - mean(hy / pi) / N;
//     NumericVector beta_new = solve(crossprod(xdif, multi3(xdif, pi, 1)),
//                                    colSums(multi2(xdif, ydif / pi)));
//     NumericVector beta_newt = clone(beta_new);
//     beta_new.push_front(max(as<NumericVector>(eRes(y - multi(x_c, beta_newt),
//                                                    delta, pi, ind_km)[1])));
//     double e = sqrt(sum(pow(beta_new - beta, 2)));
//     if (e < 1e-5)
//     {
//       break;
//     }
//     else
//     {
//       beta = beta_new;
//     }
//   }
//   List out = List::create(beta);
//   return out;
// }

// NumericMatrix resp(NumericMatrix x, NumericVector y, NumericVector delta,
//                    NumericVector beta, NumericVector pi, int N, int B)
// {
//   Function crossprod("crossprod");
//   Function Var("var");
//   int p = x.ncol() - 1;
//   int r = x.nrow();
//   IntegerVector ind_km = seq_len(r);
//   NumericMatrix Z(r, B, rexp(r * B, 1));
//   NumericMatrix g = uls(x, y, delta, pi, beta, ind_km, N);
//   NumericMatrix V = Var(crossprod(Z, g));
//   NumericMatrix zbs(B, p, ifelse(as<logicalVector>(rbinom(B * p, 1, 0.5)), 1, -1));
//   NumericMatrix zb = zbs % * % (V % ^ % (-0.5));
//   NumericMatrix beta1 = zb / sqrt(r) +
// }

// int BinarySearch(NumericVector v, int low, int high, double x)
// {
//   if (x >= v[high - 1])
//   {
//     return high;
//   }
//   int mid = (low + high) / 2;
//   if (v[mid] == x)
//   {
//     return mid;
//   }
//   if (mid > 0 && v[mid - 1] <= x && x < v[mid])
//   {
//     return mid - 1;
//   }
//   if (x < v[mid])
//   {
//     return BinarySearch(v, low, mid - 1, x);
//   }
//   else
//   {
//     return BinarySearch(v, mid + 1, high, x);
//   }
// }

// // [[Rcpp::export]]

// NumericVector pred(NumericVector vf, NumericVector vs,
//                    NumericVector km)
// {
//   int t = vf.length();
//   int p = vs.length();
//   NumericVector out(t);
//   for (int i = 0; i < t; i++)
//   {
//     double vfi = vf[i];
//     if (vfi < vs[0])
//     {
//       out[i] = 0;
//     }
//     else if (vfi > vs[p - 1])
//     {
//       out[i] = 1;
//     }
//     else
//     {
//       int m = BinarySearch(vs, 0, p, vfi);
//       out[i] = km[m];
//     }
//   }
//   return out;
// }
