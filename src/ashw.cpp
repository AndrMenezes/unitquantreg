#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// log-pdf ashw
inline double logpdf_ashw(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -log(tau) * pow(log((0.1e1 + sqrt( (1 - mu * mu))) / mu), -theta);
  double t1 = log(alpha);
  double t3 = log(0.1e1 / x + sqrt(pow(0.1e1 / x, 0.2e1) - 0.1e1));
  double t4 = log(t3);
  double t6 = log(theta);
  double t8 = log(0.1e1 - x);
  double t10 = lnx;
  double t12 = log(0.1e1 + x);
  double t14 = pow(t3, theta);
  return(t1 + t4 * (theta - 0.1e1) + t6 - t8 / 0.2e1 - t10 - t12 / 0.2e1 - t14 * alpha);
}

// [[Rcpp::export]]
NumericVector cpp_dashw(const NumericVector x,
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
    out[i] = logpdf_ashw(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out);
  else return(Rcpp::exp(out));
}

// cdf ashw
inline double cdf_ashw(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) * pow(log((0.1e1 + sqrt( (1 - mu * mu))) / mu), -theta);
  double t1 = log(0.1e1 / x + sqrt(pow(0.1e1 / x, 0.2e1) - 0.1e1));
  double t2 = pow(t1, theta);
  return(exp(-alpha * t2));
}

// [[Rcpp::export]]
NumericVector cpp_pashw(const NumericVector x,
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

  for(int i = 0; i < n; i++)
    out[i] = cdf_ashw(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf ashw
inline double invcdf_ashw(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) * pow(log((0.1e1 + sqrt( (1 - mu * mu))) / mu), -theta);
  double t1 = log(x);
  double t4 = log(-t1 / alpha);
  double t7 = exp(t4 / theta);
  return(0.1e1 / cosh(t7));
}

// [[Rcpp::export]]
NumericVector cpp_qashw(const NumericVector x,
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
      out[i] = invcdf_ashw(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ashw(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood ashw

// [[Rcpp::export]]
double cpp_loglikeashw(NumericVector x,
                           NumericVector lnx,
                           int n,
                           NumericVector mu,
                           NumericVector theta,
                           double tau)
{
  double ll = 0;

  for (int i = 0; i < n; i++)
    ll += logpdf_ashw(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// ashw score vector

// [[Rcpp::export]]
NumericMatrix cpp_gradientashw(int n,
                                   NumericVector x,
                                   NumericMatrix U,
                                   NumericVector dmu_deta,
                                   NumericVector dtheta_dzeta,
                                   NumericVector mu,
                                   NumericVector theta,
                                   double tau)
{
  double t18 = log(tau);
  for (int i = 0; i < n; i++)
  {
    double t1 = mu[i] * mu[i];
    double t3 = sqrt(0.1e1 - t1);
    double t5 = 0.1e1 + t3;
    double t9 = theta[i] * (-0.1e1 / t3 - t5 / t1);
    double t11 = 0.1e1 / t5 * mu[i];
    double t14 = log(t5 / mu[i]);
    double t15 = 0.1e1 / t14;
    double t19 = pow(t14, -theta[i]);
    double t20 = t18 * t19;
    double t22 = log(0.1e1 / x[i] + sqrt(pow(0.1e1 / x[i], 0.2e1) - 0.1e1));
    double t23 = pow(t22, theta[i]);
    double t28 = log(t14);
    double t29 = log(t22);
    double t30 = -t9 * t11 * t15 - t20 * t9 * t11 * t15 * t23;
    double t31 = -t28 + t29 + 0.1e1 / theta[i] - t20 * t28 * t23 + t20 * t23 * t29;
    U(i,0) = dmu_deta[i] * t30;
    U(i,1) = dtheta_dzeta[i] * t31;
  }
  return(U);
}

// ashw hessian matrix

// [[Rcpp::export]]
NumericMatrix cpp_hessianashw(int n,
                                  NumericVector x,
                                  NumericMatrix H,
                                  NumericVector mu,
                                  NumericVector theta,
                                  double tau)
{
  double t43 = log(tau);
  for (int i = 0; i < n; i++)
  {
    double t1 = mu[i] * mu[i];
    double t2 = 0.1e1 - t1;
    double t3 = sqrt(t2);
    double t7 = 0.1e1 / t3;
    double t8 = 0.1e1 / mu[i];
    double t10 = 0.1e1 + t3;
    double t16 = theta[i] * (-0.1e1 / t3 / t2 * mu[i] + t7 * t8 + 0.2e1 * t10 / t1 / mu[i]);
    double t17 = 0.1e1 / t10;
    double t18 = t17 * mu[i];
    double t20 = log(t10 * t8);
    double t21 = 0.1e1 / t20;
    double t26 = -t7 - t10 / t1;
    double t27 = theta[i] * t26;
    double t28 = t10 * t10;
    double t29 = 0.1e1 / t28;
    double t36 = t26 * t26;
    double t37 = theta[i] * t36;
    double t38 = t29 * t1;
    double t39 = t20 * t20;
    double t40 = 0.1e1 / t39;
    double t44 = pow(t20, -theta[i]);
    double t45 = t43 * t44;
    double t46 = theta[i] * theta[i];
    double t49 = log(0.1e1 / x[i] + sqrt(pow(0.1e1 / x[i], 0.2e1) - 0.1e1));
    double t50 = pow(t49, theta[i]);
    double t52 = t38 * t40 * t50;
    double t55 = t21 * t50;
    double t56 = t18 * t55;
    double t58 = t45 * t27;
    double t63 = t26 * t17;
    double t71 = log(t20);
    double t78 = log(t49);
    double t82 = -t63 * mu[i] * t21 + t58 * t18 * t21 * t71 * t50 - t45 * t26 * t56 - t58 * t18 * t55 * t78;
    double t84 = t71 * t71;
    double t91 = t78 * t78;
    H(i,0) = -t16 * t18 * t21 - t27 * t29 * t1 * t21 * t7 - t27 * t17 * t21 + t37 * t38 * t40 + t45 * t46 * t36 * t52 - t45 * t16 * t56 - t58 * t38 * t55 * t7 - t45 * theta[i] * t63 * t55 + t45 * t37 * t52;
    H(i,1) = t82;
    H(i,2) = -0.1e1 / t46 + t45 * t84 * t50 - 0.2e1 * t45 * t71 * t50 * t78 + t45 * t50 * t91;
  }
  return(H);
}

