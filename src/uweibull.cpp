#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf unit-weibull
inline double logpdf_uweibull(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -log(tau) / pow(-log(mu), theta);
  double t1 = log(alpha);
  double t2 = lnx;
  double t3 = log(-t2);
  double t5 = log(theta);
  double t7 = log(-0.1e1 / t2);
  double t8 = pow(-t2, theta);
  return(t1 + theta * t3 + t5 - t2 + t7 - alpha * t8);
}

// [[Rcpp::export]]
NumericVector cpp_duweibull(const NumericVector x,
                            const NumericVector mu,
                            const NumericVector theta,
                            const double tau,
                            const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_uweibull(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out);
  else return(Rcpp::exp(out));
}

// cdf unit-weibull
inline double cdf_uweibull(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / pow(-log(mu), theta);
  double t1 = log(x);
  double t2 = pow(-t1, theta);
  return(exp(-alpha * t2));
}

// [[Rcpp::export]]
NumericVector cpp_puweibull(const NumericVector x,
                            const NumericVector mu,
                            const NumericVector theta,
                            const double tau,
                            const bool lowertail = true,
                            const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_uweibull(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf unit-weibull
inline double invcdf_uweibull(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / pow(-log(mu), theta);
  double t1 = 0.1e1 / theta;
  double t2 = pow(alpha, t1);
  double t4 = log(x);
  double t5 = pow(-t4, t1);
  double t7 = exp(0.1e1 / t2 * t5);
  return(0.1e1 / t7);
}

// [[Rcpp::export]]
NumericVector cpp_quweibull(const NumericVector x,
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
      out[i] = invcdf_uweibull(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_uweibull(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood unit-weibull

// [[Rcpp::export]]
double cpp_loglikeuweibull(NumericVector x,
                           NumericVector lnx,
                           int n,
                           NumericVector mu,
                           NumericVector theta,
                           double tau)
{
  double ll = 0;

  for (int i = 0; i < n; i++)
    ll += logpdf_uweibull(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// unit-weibull score vector

// [[Rcpp::export]]
NumericMatrix cpp_gradientuweibull(int n,
                                   NumericVector x,
                                   NumericMatrix U,
                                   NumericVector dmu_deta,
                                   NumericVector dtheta_dzeta,
                                   NumericVector mu,
                                   NumericVector theta,
                                   double tau)
{
  double t8 = log(tau);
  for (int i = 0; i < n; i++)
  {
    double t3 = log(mu[i]);
    double t5 = theta[i] / mu[i] / t3;
    double t6 = pow(-t3, theta[i]);
    double t9 = 0.1e1 / t6 * t8;
    double t10 = log(x[i]);
    double t11 = pow(-t10, theta[i]);
    double t15 = log(-t3);
    double t16 = log(-t10);
    double t17 = -t5 - t9 * t11 * t5;
    double t18 = -t15 + t16 + 0.1e1 / theta[i] - t9 * t11 * t15 + t9 * t11 * t16;
    U(i,0) = t17 * dmu_deta[i];
    U(i,1) = t18 * dtheta_dzeta[i];
  }
  return(U);
}

// unit-weibull hessian matrix

// [[Rcpp::export]]
NumericMatrix cpp_hessianuweibull(int n,
                                  NumericVector x,
                                  NumericMatrix H,
                                  NumericVector mu,
                                  NumericVector theta,
                                  double tau)
{
  double t12 = log(tau);
  for (int i = 0; i < n; i++)
  {
    double t1 = mu[i] * mu[i];
    double t2 = 0.1e1 / t1;
    double t3 = theta[i] * t2;
    double t4 = log(mu[i]);
    double t5 = 0.1e1 / t4;
    double t6 = t3 * t5;
    double t7 = t4 * t4;
    double t8 = 0.1e1 / t7;
    double t9 = t3 * t8;
    double t10 = pow(-t4, theta[i]);
    double t13 = 0.1e1 / t10 * t12;
    double t14 = log(x[i]);
    double t15 = pow(-t14, theta[i]);
    double t16 = t13 * t15;
    double t17 = theta[i] * theta[i];
    double t24 = 0.1e1 / mu[i];
    double t25 = t24 * t5;
    double t26 = log(-t4);
    double t33 = log(-t14);
    double t37 = -t25 + t16 * t26 * theta[i] * t25 - t13 * t15 * t24 * t5 - t16 * t33 * theta[i] * t25;
    double t39 = t26 * t26;
    double t46 = t33 * t33;
    H(i,0) = (t6 + t9 + t16 * t17 * t2 * t8 + t16 * t6 + t16 * t9);
    H(i,1) = t37;
    H(i,2) = (-0.1e1 / t17 + t13 * t15 * t39 - 0.2e1 * t13 * t15 * t33 * t26 + t13 * t15 * t46);
  }
  return(H);
}

