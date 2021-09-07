#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf unit-Gompertz

inline double logpdf_ugompertz(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (0.1e1 / pow(mu, theta) - 0.1e1);
  double t1 = log(alpha);
  double t2 = log(theta);
  double t3 = lnx;
  double t5 = pow(x, theta);
  return(t1 + t2 - theta * t3 + alpha - alpha / t5 - t3);
}

// [[Rcpp::export]]
NumericVector cpp_dugompertz(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_ugompertz(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf unit-gompertz
inline double cdf_ugompertz(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (0.1e1 / pow(mu, theta) - 0.1e1);
  double t1 = pow(x, -theta);
  return(exp(alpha - alpha * t1));
}

// [[Rcpp::export]]
NumericVector cpp_pugompertz(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_ugompertz(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);
  return(out);
}

// inv-cdf unit-gompertz
inline double invcdf_ugompertz(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (0.1e1 / pow(mu, theta) - 0.1e1);
  double t1 = log(x);
  double t6 = pow((alpha - t1) / alpha, 0.1e1 / theta);
  return(0.1e1 / t6);
}

// [[Rcpp::export]]
NumericVector cpp_qugompertz(const NumericVector x,
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
      out[i] = invcdf_ugompertz(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ugompertz(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood unit-gompertz

// [[Rcpp::export]]
double cpp_loglikeugompertz(NumericVector x,
                            NumericVector lnx,
                            int n,
                            NumericVector mu,
                            NumericVector theta,
                            double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_ugompertz(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// unit-gompertz score vector

// [[Rcpp::export]]
NumericMatrix cpp_gradientugompertz(int n,
                                    NumericVector x,
                                    NumericMatrix U,
                                    NumericVector dmu_deta,
                                    NumericVector dtheta_dzeta,
                                    NumericVector mu,
                                    NumericVector theta,
                                    double tau)
{
  double t8 = log(tau);
  for(int i = 0; i < n; i++)
  {
    double t1 = pow(mu[i], -theta[i]);
    double t2 = t1 - 0.1e1;
    double t3 = 0.1e1 / t2;
    double t4 = t3 * t1;
    double t5 = 0.1e1 / mu[i];
    double t9 = t2 * t2;
    double t11 = t8 / t9;
    double t13 = t1 * theta[i] * t5;
    double t15 = pow(x[i], theta[i]);
    double t16 = 0.1e1 / t15;
    double t20 = log(mu[i]);
    double t23 = log(x[i]);
    U(i,0) = (t4 * theta[i] * t5 - t11 * t13 + t11 * t16 * t13) * dmu_deta[i];
    U(i,1) = (t4 * t20 + 0.1e1 / theta[i] - t23 - t11 * t1 * t20 + t11 * t16 * t1 * t20 - t8 * t3 * t16 * t23) * dtheta_dzeta[i];
  }
  return(U);
}

// unit-gompertz hessian matrix

// [[Rcpp::export]]
NumericMatrix cpp_hessianugompertz(int n,
                                   NumericVector x,
                                   NumericMatrix H,
                                   NumericVector mu,
                                   NumericVector theta,
                                   double tau)
{
  double t17 = log(tau);
  for(int i = 0; i < n; i++)
  {
    double t1 = pow(mu[i], -theta[i]);
    double t2 = t1 - 0.1e1;
    double t3 = t2 * t2;
    double t4 = 0.1e1 / t3;
    double t5 = t1 * t1;
    double t6 = t4 * t5;
    double t7 = theta[i] * theta[i];
    double t8 = mu[i] * mu[i];
    double t9 = 0.1e1 / t8;
    double t10 = t7 * t9;
    double t12 = 0.1e1 / t2;
    double t13 = t12 * t1;
    double t20 = t17 / t3 / t2;
    double t22 = t5 * t7 * t9;
    double t25 = t17 * t4;
    double t27 = t1 * t7 * t9;
    double t29 = t1 * theta[i];
    double t30 = t29 * t9;
    double t32 = pow(x[i], theta[i]);
    double t33 = 0.1e1 / t32;
    double t34 = t20 * t33;
    double t37 = t25 * t33;
    double t41 = log(mu[i]);
    double t43 = 0.1e1 / mu[i];
    double t44 = t41 * theta[i] * t43;
    double t56 = theta[i] * t43;
    double t63 = t33 * t1;
    double t66 = log(x[i]);
    double t70 = t6 * t44 - t13 * t44 + t13 * t43 - 0.2e1 * t20 * t5 * t44 + t25 * t1 * t44 - t25 * t1 * t43 + 0.2e1 * t34 * t5 * t41 * t56 - t37 * t29 * t43 * t41 + t25 * t63 * t43 - t37 * t66 * t1 * t56;
    double t71 = t41 * t41;
    double t91 = t66 * t66;
    H(i,0) = t6 * t10 - t13 * t10 - t13 * theta[i] * t9 - 0.2e1 * t20 * t22 + t25 * t27 + t25 * t30 + 0.2e1 * t34 * t22 - t37 * t27 - t37 * t30;
    H(i,1) = t70;
    H(i,2) = t6 * t71 - t13 * t71 - 0.1e1 / t7 - 0.2e1 * t20 * t5 * t71 + t25 * t1 * t71 + 0.2e1 * t20 * t33 * t5 * t71 - 0.2e1 * t37 * t1 * t41 * t66 - t25 * t63 * t71 + t17 * t12 * t33 * t91;
  }
  return(H);
}

