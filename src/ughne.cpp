#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf ughne

inline double logpdf_ughne(double x, double lnx, double mu, double theta, double tau)
{
  double rqnorm = R::qnorm(tau / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  double alpha = -log(mu) / pow(-rqnorm, 0.1e1 / theta);
  double t1 = log(theta);
  double t2 = lnx;
  double t5 = log(-0.1e1 / t2);
  double t6 = t2;
  double t7 = log(-t6);
  double t8 = log(alpha);
  double t13 = pow(alpha, -theta);
  double t14 = t13 * t13;
  double t15 = pow(-t6, theta);
  double t16 = t15 * t15;
  return(-0.2257913526e0 + t1 + t5 + theta * (t7 - t8) - 0.5e0 * t14 * t16 - t6);
}

// [[Rcpp::export]]
NumericVector cpp_dughne(const NumericVector x,
                         const NumericVector mu,
                         const NumericVector theta,
                         const double tau,
                         const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_ughne(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf ughne
inline double cdf_ughne(double x, double mu, double theta, double tau)
{
  double rqnorm = R::qnorm(tau / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  double alpha = -log(mu) / pow(-rqnorm, 0.1e1 / theta);
  double rpnorm = R::pnorm(-pow(-log(x) / alpha, theta) , 0.0, 1.0, TRUE, FALSE);
  return(0.2e1 * rpnorm);
}

// [[Rcpp::export]]
NumericVector cpp_pughne(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_ughne(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf ughne
inline double invcdf_ughne(double x, double mu, double theta, double tau)
{
  double rqnorm1 = R::qnorm(tau / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  double alpha = -log(mu) / pow(-rqnorm1, 0.1e1 / theta);
  double rqnorm2 = pow(-R::qnorm(x / 0.2e1, 0.0, 1.0, TRUE, FALSE), 0.1e1 / theta);
  return(exp(-alpha * rqnorm2));
}

// [[Rcpp::export]]
NumericVector cpp_qughne(const NumericVector x,
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
      out[i] = invcdf_ughne(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ughne(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood ughne

// [[Rcpp::export]]
double cpp_loglikeughne(NumericVector x,
                            NumericVector lnx,
                            int n,
                            NumericVector mu,
                            NumericVector theta,
                            double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_ughne(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// score-vector ughne

// [[Rcpp::export]]
NumericMatrix cpp_gradientughne(int n,
                                NumericVector x,
                                NumericMatrix U,
                                NumericVector dmu_deta,
                                NumericVector dtheta_dzeta,
                                NumericVector mu,
                                NumericVector theta,
                                double tau)
{
  double rqnorm = R::qnorm(tau / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t1 = 0.1e1 / mu[i];
    double t3 = log(mu[i]);
    double t4 = 0.1e1 / t3;
    double t7 = 0.1e1 / theta[i];
    double t8 = pow(-rqnorm, t7);
    double t11 = t3 / t8;
    double t12 = 0.2e1 * theta[i];
    double t13 = pow(-t11, -t12);
    double t16 = log(x[i]);
    double t17 = t16;
    double t18 = pow(-t17, t12);
    double t23 = log(-t16);
    double t24 = log(-t11);
    double t25 = log(-rqnorm);
    double t26 = t7 * t25;
    double t35 = log(-t17);
    double t36 = -0.1000000000e1 * theta[i] * t1 * t4 + 0.1000000000e1 * t13 * theta[i] * t1 * t4 * t18;
    double t37 = t7 + t23 - t24 - 0.1000000000e1 * t26 - 0.5e0 * t13 * (-0.2e1 * t24 - 0.2000000000e1 * t26) * t18 - 0.10e1 * t13 * t18 * t35;
    U(i,0) = t36 * dmu_deta[i];
    U(i,1) = t37 * dtheta_dzeta[i];
  }
  return(U);
}

// hessian-matrix ughne

// [[Rcpp::export]]
NumericMatrix cpp_hessianughne(int n,
                               NumericVector x,
                               NumericMatrix H,
                               NumericVector mu,
                               NumericVector theta,
                               double tau)
{
  double rqnorm = R::qnorm(tau / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t1 = mu[i] * mu[i];
    double t2 = 0.1e1 / t1;
    double t3 = theta[i] * t2;
    double t4 = log(mu[i]);
    double t5 = 0.1e1 / t4;
    double t8 = t4 * t4;
    double t9 = 0.1e1 / t8;
    double t12 = 0.1e1 / theta[i];
    double t13 = pow(-rqnorm, t12);
    double t16 = t4 / t13;
    double t17 = 0.2e1 * theta[i];
    double t18 = pow(-t16, -t17);
    double t19 = theta[i] * theta[i];
    double t22 = log(x[i]);
    double t23 = t22;
    double t24 = pow(-t23, t17);
    double t25 = t2 * t9 * t24;
    double t28 = t18 * theta[i];
    double t36 = 0.1e1 / mu[i];
    double t39 = t28 * t36;
    double t40 = log(-t16);
    double t42 = log(-rqnorm);
    double t45 = -0.2e1 * t40 - 0.2000000000e1 * t12 * t42;
    double t51 = t5 * t24;
    double t54 = log(-t23);
    double t58 = -0.1000000000e1 * t36 * t5 + 0.1000000000e1 * t39 * t5 * t45 * t24 + 0.1000000000e1 * t18 * t36 * t51 + 0.2000000000e1 * t39 * t51 * t54;
    double t60 = t45 * t45;
    double t69 = t54 * t54;
    H(i,0) = 0.1000000000e1 * t3 * t5 + 0.1000000000e1 * t3 * t9 - 0.2000000000e1 * t18 * t19 * t25 - 0.1000000000e1 * t28 * t2 * t5 * t24 - 0.1000000000e1 * t28 * t25;
    H(i,1) = t58;
    H(i,2) = -0.1e1 / t19 - 0.5e0 * t18 * t60 * t24 - 0.20e1 * t18 * t45 * t24 * t54 - 0.20e1 * t18 * t24 * t69;
  }
  return(H);
}
