#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// log-pdf kum
inline double logpdf_kum(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = log(0.1e1 - tau) / log(-pow(mu, theta) + 0.1e1);
  double t1 = pow(x, theta);
  double t3 = log(0.1e1 - t1);
  double t5 = log(alpha);
  double t6 = lnx;
  double t8 = log(theta);
  return(alpha * t3 + t5 + theta * t6 + t8 - t6 - t3);
}

// [[Rcpp::export]]
NumericVector cpp_dkum(const NumericVector x,
                       const NumericVector mu,
                       const NumericVector theta,
                       const double tau,
                       const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_kum(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf kum
inline double cdf_kum(double x, double mu, double theta, double tau)
{
  double alpha = log(0.1e1 - tau) / log(-pow(mu, theta) + 0.1e1);
  double t1 = pow(x, theta);
  double t3 = pow(0.1e1 - t1, alpha);
  return(0.1e1 - t3);
}

// [[Rcpp::export]]
NumericVector cpp_pkum(const NumericVector x,
                       const NumericVector mu,
                       const NumericVector theta,
                       const double tau,
                       const bool lowertail = true,
                       const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_kum(x[i], mu[i % nmu], theta[i % nth], tau);

    if (!lowertail) out = 0.1e1 - out;
    if (logprob) out = Rcpp::log(out);

    return(out);
}

// inv-cdf kum
inline double invcdf_kum(double x, double mu, double theta, double tau)
{
  double alpha = log(0.1e1 - tau) / log(-pow(mu, theta) + 0.1e1);
  double t3 = pow(0.1e1 - x, 0.1e1 / alpha);
  return(pow(0.1e1 - t3, 0.1e1 / theta));
}

// [[Rcpp::export]]
NumericVector cpp_qkum(const NumericVector x,
                       const NumericVector mu,
                       const NumericVector theta,
                       const double tau,
                       const bool lowertail = true,
                       const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  if(lowertail)
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_kum(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_kum(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood kum

// [[Rcpp::export]]
double cpp_loglikekum(NumericVector x,
                      NumericVector lnx,
                      int n,
                      NumericVector mu,
                      NumericVector theta,
                      double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_kum(x[i], lnx[i], mu[i], theta[i], tau);
    return(-ll);
}

// kum score vector

// [[Rcpp::export]]
NumericMatrix cpp_gradientkum(int n,
                              NumericVector x,
                              NumericMatrix U,
                              NumericVector dmu_deta,
                              NumericVector dtheta_dzeta,
                              NumericVector mu,
                              NumericVector theta,
                              double tau)
{
  for(int i = 0; i < n; i++)
  {
    double t2 = log(1 - tau);
    double t3 = pow(mu[i], theta[i]);
    double t4 = -t3 + 0.1e1;
    double t5 = log(t4);
    double t6 = t5 * t5;
    double t9 = pow(x[i], theta[i]);
    double t10 = 0.1e1 - t9;
    double t11 = log(t10);
    double t12 = t2 / t6 * t11;
    double t14 = 0.1e1 / mu[i];
    double t15 = 0.1e1 / t4;
    double t19 = 0.1e1 / t5;
    double t20 = t19 * t3;
    double t25 = log(mu[i]);
    double t30 = log(x[i]);
    double t33 = t9 * t30 / t10;
    double t34 = t12 * t3 * theta[i] * t14 * t15 + t20 * theta[i] * t14 * t15;
    double t35 = t12 * t3 * t25 * t15 - t2 * t19 * t33 + t20 * t25 * t15 + t30 + 0.1e1 / theta[i] + t33;
    U(i,0) = t34 * dmu_deta[i];
    U(i,1) = t35 * dtheta_dzeta[i];
  }
  return(U);
}

// kum hessian matrix

// [[Rcpp::export]]
NumericMatrix cpp_hessiankum(int n,
                             NumericVector x,
                             NumericMatrix H,
                             NumericVector mu,
                             NumericVector theta,
                             double tau)
{
  for(int i = 0; i < n; i++)
  {
    double t2 = log(0.1e1 - tau);
    double t3 = pow(mu[i], theta[i]);
    double t4 = -t3 + 0.1e1;
    double t5 = log(t4);
    double t6 = t5 * t5;
    double t9 = t2 / t6 / t5;
    double t10 = pow(x[i], theta[i]);
    double t11 = 0.1e1 - t10;
    double t12 = log(t11);
    double t13 = t9 * t12;
    double t14 = t3 * t3;
    double t15 = theta[i] * theta[i];
    double t17 = mu[i] * mu[i];
    double t18 = 0.1e1 / t17;
    double t19 = t4 * t4;
    double t20 = 0.1e1 / t19;
    double t22 = t14 * t15 * t18 * t20;
    double t25 = 0.1e1 / t6;
    double t26 = t2 * t25;
    double t27 = t26 * t12;
    double t29 = 0.1e1 / t4;
    double t30 = t18 * t29;
    double t37 = t25 * t14;
    double t38 = t15 * t18;
    double t39 = t38 * t20;
    double t41 = 0.1e1 / t5;
    double t42 = t41 * t3;
    double t48 = t41 * t14;
    double t51 = t12 * t14;
    double t53 = log(mu[i]);
    double t55 = 0.1e1 / mu[i];
    double t56 = theta[i] * t55;
    double t57 = t53 * t20 * t56;
    double t62 = t53 * t29;
    double t70 = log(x[i]);
    double t72 = t26 * t10 * t70;
    double t73 = 0.1e1 / t11;
    double t74 = t73 * t3;
    double t80 = t20 * theta[i] * t55;
    double t90 = 0.2e1 * t9 * t51 * t57 + t26 * t12 * t3 * t56 * t62 + t27 * t3 * t55 * t29 + t26 * t51 * t57 - t72 * t74 * t56 * t29 + t37 * t53 * t80 + t42 * theta[i] * t55 * t53 * t29 + t42 * t55 * t29 + t48 * t53 * t80;
    double t91 = t53 * t53;
    double t93 = t14 * t91 * t20;
    double t103 = t2 * t41;
    double t104 = t70 * t70;
    double t106 = t10 * t104 * t73;
    double t108 = t10 * t10;
    double t110 = t11 * t11;
    double t112 = t108 * t104 / t110;
    double t114 = t91 * t20;
    double t120 = 0.2e1 * t13 * t93 - 0.2e1 * t72 * t74 * t62 + t27 * t3 * t91 * t29 + t27 * t93 - t103 * t106 - t103 * t112 + t37 * t114 + t42 * t91 * t29 + t48 * t114 - 0.1e1 / t15 + t106 + t112;
    H(i,0) = 0.2e1 * t13 * t22 + t27 * t3 * t15 * t30 - t27 * t3 * theta[i] * t30 + t27 * t22 + t37 * t39 + t42 * t38 * t29 - t42 * theta[i] * t18 * t29 + t48 * t39;
    H(i,1) = t90;
    H(i,2) = t120;
  }
  return(H);
}
