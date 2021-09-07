#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf leeg

inline double logpdf_leeg(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -0.1e1 / (-0.1e1 + tau) * (-0.1e1 + pow(mu, -theta) * tau);
  double t1 = pow(x, theta);
  double t4 = log(0.1e1 + alpha * t1);
  double t7 = log(0.1e1 + alpha);
  double t8 = log(theta);
  double t9 = lnx;
  return(-0.2e1 * t4 + t7 + t8 + theta * t9 - t9);
}

// [[Rcpp::export]]
NumericVector cpp_dleeg(const NumericVector x,
                        const NumericVector mu,
                        const NumericVector theta,
                        const double tau,
                        const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_leeg(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf leeg
inline double cdf_leeg(double x, double mu, double theta, double tau)
{
  double alpha = -0.1e1 / (-0.1e1 + tau) * (-0.1e1 + pow(mu, -theta) * tau);
  double t2 = pow(x, theta);

  return((0.1e1 + alpha) * t2 / (0.1e1 + alpha * t2));
}

// [[Rcpp::export]]
NumericVector cpp_pleeg(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_leeg(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf leeg
inline double invcdf_leeg(double x, double mu, double theta, double tau)
{
  double alpha = -0.1e1 / (-0.1e1 + tau) * (-0.1e1 + pow(mu, -theta) * tau);
  double t1 = 0.1e1 / theta;
  double t2 = pow(x, t1);
  double t6 = pow(-0.1e1 / (-0.1e1 - alpha + x * alpha), t1);
  return(t2 * t6);
}

// [[Rcpp::export]]
NumericVector cpp_qleeg(const NumericVector x,
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
      out[i] = invcdf_leeg(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_leeg(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }

  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood leeg

// [[Rcpp::export]]
double cpp_loglikeleeg(NumericVector x, NumericVector lnx, int n, NumericVector mu, NumericVector theta, double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_leeg(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// score vector leeg

// [[Rcpp::export]]
NumericMatrix cpp_gradientleeg(int n,
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
    double t1 = pow(mu[i], -theta[i]);
    double t4 = t1 * theta[i] / mu[i];
    double t6 = 0.1e1 / (-0.1e1 + tau);
    double t7 = tau * t6;
    double t8 = pow(x[i], theta[i]);
    double t11 = (t1 * tau - 0.1e1) * t6;
    double t14 = 0.1e1 / (0.1e1 - t11 * t8);
    double t21 = t7 / (0.1e1 - t11);
    double t24 = log(mu[i]);
    double t25 = t1 * t24;
    double t28 = log(x[i]);
    double t29 = -0.2e1 * t4 * t7 * t8 * t14 + t4 * t21;
    double t30 = -0.2e1 * (t25 * t7 * t8 - t11 * t8 * t28) * t14 + t25 * t21 + 0.1e1 / theta[i] + t28;
    U(i,0) = t29 * dmu_deta[i];
    U(i,1) = t30 * dtheta_dzeta[i];
  }
  return(U);
}

// hessian matrix leeg

// [[Rcpp::export]]
NumericMatrix cpp_hessianleeg(int n,
                              NumericVector x,
                              NumericMatrix H,
                              NumericVector mu,
                              NumericVector theta,
                              double tau)
{
  for(int i = 0; i < n; i++)
  {
    double t2 = pow(mu[i], -theta[i]);
    double t3 = theta[i] * theta[i];
    double t5 = mu[i] * mu[i];
    double t6 = 0.1e1 / t5;
    double t7 = t2 * t3 * t6;
    double t8 = -0.1e1 + tau;
    double t9 = 0.1e1 / t8;
    double t10 = tau * t9;
    double t11 = pow(x[i], theta[i]);
    double t14 = (t2 * tau - 0.1e1) * t9;
    double t17 = 0.1e1 - t14 * t11;
    double t18 = 0.1e1 / t17;
    double t20 = t10 * t11 * t18;
    double t23 = t2 * theta[i];
    double t24 = t23 * t6;
    double t27 = t2 * t2;
    double t29 = t27 * t3 * t6;
    double t30 = tau * tau;
    double t31 = t8 * t8;
    double t32 = 0.1e1 / t31;
    double t33 = t30 * t32;
    double t34 = t11 * t11;
    double t35 = t17 * t17;
    double t36 = 0.1e1 / t35;
    double t42 = 0.1e1 - t14;
    double t43 = 0.1e1 / t42;
    double t44 = t10 * t43;
    double t49 = t42 * t42;
    double t50 = 0.1e1 / t49;
    double t51 = t33 * t50;
    double t55 = 0.1e1 / mu[i];
    double t56 = t23 * t55;
    double t57 = log(mu[i]);
    double t58 = t57 * tau;
    double t59 = t9 * t11;
    double t63 = t2 * t55;
    double t64 = t10 * t11;
    double t67 = log(x[i]);
    double t68 = t11 * t67;
    double t75 = t2 * t57;
    double t80 = t75 * t64 - t14 * t68;
    double t83 = t55 * tau;
    double t86 = 0.2e1 * t80 * t36 * t23 * t83 * t59;
    double t90 = t56 * t58 * t9 * t43;
    double t92 = t63 * t44;
    double t96 = theta[i] * t55;
    double t99 = t27 * t57 * t30 * t32 * t50 * t96;
    double t114 = t57 * t57;
    double t115 = t2 * t114;
    double t122 = t67 * t67;
    double t129 = t80 * t80;
    H(i,0) = 0.2e1 * t7 * t20 + 0.2e1 * t24 * t20 + 0.2e1 * t29 * t33 * t34 * t36 - t7 * t44 - t24 * t44 - t29 * t51;
    H(i,1) = 0.2e1 * t75 * t96 * t20 - 0.2e1 * t63 * tau * t59 * t18 - 0.2e1 * t23 * t83 * t59 * t67 * t18 + t86 - t90 + t92 - t99;
    H(i,2) = -0.2e1 * (-t115 * t64 + 0.2e1 * t75 * tau * t59 * t67 - t14 * t11 * t122) * t18 + 0.2e1 * t129 * t36 - t115 * t44 - t27 * t114 * t51 - 0.1e1 / t3;
  }
  return(H);
}
