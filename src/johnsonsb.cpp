#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf johnsonsb

inline double logpdf_johnsonsb(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE) - theta * log(mu / (0.1e1 - mu));
  double t1 = log(theta);
  double t2 = lnx;
  double t5 = 0.1e1 - x;
  double t6 = log(t5);
  double t10 = log(x / t5);
  double t13 = pow(alpha + theta * t10, 0.2e1);
  double t15 = exp(-0.5000000000e0 * t13);
  double t16 = log(t15);
  return(-0.918938533204672e0 + t1 - t2 - t6 + t16);
}

// [[Rcpp::export]]
NumericVector cpp_djohnsonsb(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_johnsonsb(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf johnsonsb
inline double cdf_johnsonsb(double x, double mu, double theta, double tau)
{
  double alpha = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE) - log(mu / (0.1e1 - mu)) * theta;
  return(R::pnorm(alpha + theta * log(x / (0.1e1 - x)), 0.0, 1.0, TRUE, FALSE));
}

// [[Rcpp::export]]
NumericVector cpp_pjohnsonsb(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_johnsonsb(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf johnsonsb
inline double invcdf_johnsonsb(double x, double mu, double theta, double tau)
{
  double alpha = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE) - log(mu / (0.1e1 - mu)) * theta;
  double rqnorm = exp((R::qnorm(x, 0.0, 1.0, TRUE, FALSE) - alpha) / theta);
  return(rqnorm / (0.1e1 + rqnorm));
}

// [[Rcpp::export]]
NumericVector cpp_qjohnsonsb(const NumericVector x,
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
      out[i] = invcdf_johnsonsb(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_johnsonsb(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }

  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood johnsonsb

// [[Rcpp::export]]
double cpp_loglikejohnsonsb(NumericVector x, NumericVector lnx, int n, NumericVector mu, NumericVector theta, double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_johnsonsb(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// score vector johnsonsb

// [[Rcpp::export]]
NumericMatrix cpp_gradientjohnsonsb(int n,
                                    NumericVector x,
                                    NumericMatrix U,
                                    NumericVector dmu_deta,
                                    NumericVector dtheta_dzeta,
                                    NumericVector mu,
                                    NumericVector theta,
                                    double tau)
{
  double rqnorm = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t1 = theta[i] * theta[i];
    double t2 = 0.1e1 / mu[i];
    double t3 = t1 * t2;
    double t4 = log(x[i]);
    double t9 = log(0.1e1 - x[i]);
    double t13 = 0.1e1 / (mu[i] - 0.1e1);
    double t14 = t1 * t13;
    double t19 = rqnorm * theta[i];
    double t25 = log(-t13);
    double t28 = log(mu[i]);
    double t29 = t1 * t28;
    double t39 = theta[i] * t28;
    double t44 = theta[i] * t25;
    double t62 = t28 * t28;
    double t65 = t25 * t25;
    double t68 = t4 * t4;
    double t71 = t9 * t9;
    double t74 = 0.1e1 / theta[i] + 0.2000000000e1 * t39 * t4 - 0.2000000000e1 * t39 * t9 + 0.2000000000e1 * t44 * t4 - 0.2000000000e1 * t44 * t9 + 0.2000000000e1 * theta[i] * t4 * t9 + 0.100000000000000e1 * rqnorm * t28 + 0.100000000000000e1 * rqnorm * t25 - 0.1000000000e1 * rqnorm * t4 + 0.100000000000000e1 * rqnorm * t9 - 0.2000000000e1 * t39 * t25 - 0.1000000000e1 * theta[i] * t62 - 0.1000000000e1 * theta[i] * t65 - 0.1000000000e1 * theta[i] * t68 - 0.1000000000e1 * theta[i] * t71;
    U(i,0) = (0.100000000000000e1 * t3 * t4 - 0.1000000000e1 * t3 * t9 - 0.1000000000e1 * t14 * t4 + 0.1000000000e1 * t14 * t9 + 0.100000000000000e1 * t19 * t2 - 0.1000000000e1 * t19 * t13 - 0.1000000000e1 * t3 * t25 + 0.1000000000e1 * t29 * t13 - 0.1000000000e1 * t29 * t2 + 0.1000000000e1 * t1 * t25 * t13) * dmu_deta[i];
    U(i,1) = t74 * dtheta_dzeta[i];
  }
  return(U);
}

// hessian matrix johnsonsb

// [[Rcpp::export]]
NumericMatrix cpp_hessianjohnsonsb(int n,
                                   NumericVector x,
                                   NumericMatrix H,
                                   NumericVector mu,
                                   NumericVector theta,
                                   double tau)
{
  double rqnorm = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t1 = theta[i] * theta[i];
    double t2 = mu[i] * mu[i];
    double t3 = 0.1e1 / t2;
    double t4 = t1 * t3;
    double t5 = log(x[i]);
    double t10 = log(0.1e1 - x[i]);
    double t13 = mu[i] - 0.1e1;
    double t14 = t13 * t13;
    double t15 = 0.1e1 / t14;
    double t16 = t1 * t15;
    double t21 = rqnorm * theta[i];
    double t26 = 0.1e1 / t13;
    double t28 = log(-t26);
    double t31 = 0.1e1 / mu[i];
    double t35 = log(mu[i]);
    double t36 = t1 * t35;
    double t46 = -0.1000000000e1 * t4 * t5 + 0.100000000000000e1 * t4 * t10 + 0.1000000000e1 * t16 * t5 - 0.1000000000e1 * t16 * t10 - 0.1000000000e1 * t21 * t3 + 0.1000000000e1 * t21 * t15 + 0.100000000000000e1 * t4 * t28 + 0.2000000000e1 * t1 * t31 * t26 - 0.1000000000e1 * t36 * t15 - 0.1000000000e1 * t4 + 0.1000000000e1 * t36 * t3 - 0.1000000000e1 * t16 - 0.1000000000e1 * t1 * t28 * t15;
    double t47 = theta[i] * t31;
    double t52 = theta[i] * t26;
    double t63 = theta[i] * t35;
    double t71 = 0.2000000000e1 * t47 * t5 - 0.2000000000e1 * t47 * t10 - 0.2000000000e1 * t52 * t5 + 0.2000000000e1 * t52 * t10 + 0.100000000000000e1 * rqnorm * t31 - 0.1000000000e1 * rqnorm * t26 - 0.2000000000e1 * t47 * t28 + 0.2000000000e1 * t63 * t26 - 0.2000000000e1 * t63 * t31 + 0.2000000000e1 * theta[i] * t28 * t26;
    double t85 = t35 * t35;
    double t87 = t28 * t28;
    double t89 = t5 * t5;
    double t91 = t10 * t10;
    double t93 = -0.1e1 / t1 + 0.2000000000e1 * t35 * t5 - 0.2000000000e1 * t35 * t10 + 0.2000000000e1 * t28 * t5 - 0.2000000000e1 * t28 * t10 + 0.2000000000e1 * t5 * t10 - 0.2000000000e1 * t35 * t28 - 0.1000000000e1 * t85 - 0.1000000000e1 * t87 - 0.1000000000e1 * t89 - 0.1000000000e1 * t91;
    H(i,0) = t46;
    H(i,1) = t71;
    H(i,2) = t93;
  }
  return(H);
}
