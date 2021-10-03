#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf unit-gumbel
inline double logpdf_ugumbel(double x, double lnx, double mu, double theta, double tau)
{
  double alpha =  log(-(mu - 0.1e1) / mu) * theta + log(-0.1e1 / log(tau));
  double t1 = lnx;
  double t3 = 0.1e1 - x;
  double t4 = log(t3);
  double t6 = exp(alpha);
  double t8 = pow(x, theta);
  double t11 = pow(t3, -theta);
  double t14 = log(theta);
  return(-alpha - theta * t1 + theta * t4 - 0.1e1 / t6 / t8 / t11 + t14 - t1 - t4);
}

// [[Rcpp::export]]
NumericVector cpp_dugumbel(const NumericVector x,
                            const NumericVector mu,
                            const NumericVector theta,
                            const double tau,
                            const bool logprob = false)
{
  const int n = x.length();
  NumericVector out(n);
  const int nmu = mu.length();
  const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_ugumbel(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out);
  else return(Rcpp::exp(out));
}

// cdf unit-gumbel
inline double cdf_ugumbel(double x, double mu, double theta, double tau)
{
  double alpha =  log(-(mu - 0.1e1) / mu) * theta + log(-0.1e1 / log(tau));
  double t1 = exp(-alpha);
  double t2 = pow(x, -theta);
  double t5 = pow(0.1e1 - x, theta);
  return(exp(-t1 * t2 * t5));
}

// [[Rcpp::export]]
NumericVector cpp_pugumbel(const NumericVector x,
                            const NumericVector mu,
                            const NumericVector theta,
                            const double tau,
                            const bool lowertail = true,
                            const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_ugumbel(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf unit-gumbel
inline double invcdf_ugumbel(double x, double mu, double theta, double tau)
{
  double alpha =  log(-(mu - 0.1e1) / mu) * theta + log(-0.1e1 / log(tau));
  double t1 = 0.1e1 / theta;
  double t3 = exp(t1 * alpha);
  double t4 = log(x);
  double t6 = pow(-0.1e1 / t4, t1);
  return(0.1e1 / (t3 + t6) * t6);
}

// [[Rcpp::export]]
NumericVector cpp_qugumbel(const NumericVector x,
                            const NumericVector mu,
                            const NumericVector theta,
                            const double tau,
                            const bool lowertail = true,
                            const bool logprob = false)
{
  const int n = x.length();
  NumericVector out(n);
  const int nmu = mu.length();
  const int nth = theta.length();

  if(lowertail)
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ugumbel(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ugumbel(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood unit-gumbel

// [[Rcpp::export]]
double cpp_loglikeugumbel(NumericVector x,
                           NumericVector lnx,
                           int n,
                           NumericVector mu,
                           NumericVector theta,
                           double tau)
{
  double ll = 0;

  for (int i = 0; i < n; i++)
    ll += logpdf_ugumbel(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// unit-gumbel score vector

// [[Rcpp::export]]
NumericMatrix cpp_gradientugumbel(int n,
                                   NumericVector x,
                                   NumericMatrix U,
                                   NumericVector dmu_deta,
                                   NumericVector dtheta_dzeta,
                                   NumericVector mu,
                                   NumericVector theta,
                                   double tau)
{
  double t1 = log(tau);
  for (int i = 0; i < n; i++)
  {
    double t1 = -mu[i] + 0.1e1;
    double t2 = 0.1e1 / t1;
    double t4 = 0.1e1 / mu[i];
    double t8 = log(t1);
    double t9 = log(mu[i]);
    double t10 = t8 - t9;
    double t12 = exp(-t10 * theta[i]);
    double t14 = pow(x[i], -theta[i]);
    double t15 = x[i] - 0.1e1;
    double t16 = 0.1e1 / t15;
    double t17 = pow(-t16, -theta[i]);
    double t18 = t14 * t17;
    double t19 = log(tau);
    double t24 = log(x[i]);
    double t25 = log(-t15);
    double t28 = t12 * t14;
    double t31 = log(-t16);
    double t48 = theta[i] * t2 + theta[i] * t4 - (-t2 - t4) * theta[i] * t12 * t18 * t19;
    double t49 = 0.1e1 / theta[i] - t8 + t9 - t24 + t25 + (-t10 * t12 * t18 - t28 * t24 * t17 - t28 * t17 * t31) * t19;
    U(i,0) = t48 * dmu_deta[i];
    U(i,1) = t49 * dtheta_dzeta[i];
  }
  return(U);
}

// unit-gumbel hessian matrix

// [[Rcpp::export]]
NumericMatrix cpp_hessianugumbel(int n,
                                  NumericVector x,
                                  NumericMatrix H,
                                  NumericVector mu,
                                  NumericVector theta,
                                  double tau)
{
  double t10 = log(tau);
  for (int i = 0; i < n; i++)
  {
    double t1 = -mu[i] + 0.1e1;
    double t2 = t1 * t1;
    double t3 = 0.1e1 / t2;
    double t5 = mu[i] * mu[i];
    double t6 = 0.1e1 / t5;
    double t10 = log(t1);
    double t11 = log(mu[i]);
    double t12 = t10 - t11;
    double t14 = exp(-t12 * theta[i]);
    double t16 = pow(x[i], -theta[i]);
    double t18 = 0.1e1 / (x[i] - 0.1e1);
    double t19 = pow(-t18, -theta[i]);
    double t20 = t16 * t19;
    double t21 = log(tau);
    double t22 = t20 * t21;
    double t24 = 0.1e1 / t1;
    double t25 = 0.1e1 / mu[i];
    double t26 = -t24 - t25;
    double t27 = t26 * t26;
    double t28 = theta[i] * theta[i];
    double t36 = -t12 * t26 * theta[i];
    double t37 = t14 * t16;
    double t41 = t26 * theta[i] * t14;
    double t42 = log(x[i]);
    double t43 = t16 * t42;
    double t44 = t43 * t19;
    double t46 = log(-t18);
    double t47 = t20 * t46;
    double t54 = t19 * t21;
    double t64 = t12 * t12;
    double t67 = -t12 * t14;
    double t72 = t42 * t42;
    double t79 = t46 * t46;
    H(i,0) = theta[i] * t3 - theta[i] * t6 - (-t3 + t6) * theta[i] * t14 * t22 + t27 * t28 * t14 * t22;
    H(i,1) = t24 + t25 + (-t26 * t14 * t20 - t36 * t37 * t19 + t41 * t44 + t41 * t47) * t21;
    H(i,2) = -0.1e1 / t28 + (t64 * t14 * t20 - 0.2e1 * t67 * t44 - 0.2e1 * t67 * t47 + t37 * t72 * t19 + 0.2e1 * t37 * t42 * t19 * t46 + t37 * t19 * t79) * t21;
  }
  return(H);
}

