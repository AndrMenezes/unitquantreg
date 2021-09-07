#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf unit-logistic
inline double logpdf_ulogistic(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = log(tau) - log(0.1e1 - tau) - log(-(mu / (-0.1e1 + mu))) * theta;
  double t1 = lnx;
  double t3 = 0.1e1 / (-0.1e1 + x);
  double t4 = log(-t3);
  double t7 = log(theta);
  double t8 = exp(alpha);
  double t10 = pow(-x * t3, theta);
  double t13 = log(0.1e1 + t8 * t10);
  return(alpha + theta * (t1 + t4) + t7 - 0.2e1 * t13 - t1 + t4);
}

// [[Rcpp::export]]
NumericVector cpp_dulogistic(const NumericVector x,
                            const NumericVector mu,
                            const NumericVector theta,
                            const double tau,
                            const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_ulogistic(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out);
  else return(Rcpp::exp(out));
}

// cdf unit-logistic
inline double cdf_ulogistic(double x, double mu, double theta, double tau)
{
  double alpha = log(tau) - log(0.1e1 - tau) - log(-(mu / (-0.1e1 + mu))) * theta;
  double t1 = exp(alpha);
  double t5 = pow(-(x / (-1 + x)), theta);
  double t6 = t1 * t5;
  return(t6 / (0.1e1 + t6));
}

// [[Rcpp::export]]
NumericVector cpp_pulogistic(const NumericVector x,
                            const NumericVector mu,
                            const NumericVector theta,
                            const double tau,
                            const bool lowertail = true,
                            const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_ulogistic(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf unit-logistic
inline double invcdf_ulogistic(double x, double mu, double theta, double tau)
{
  double alpha = -(-0.1e1 + tau) * pow(-(-0.1e1 + mu) / mu, -theta) / tau;
  double t7 = pow((0.1e1 - x) / x / alpha, 0.1e1 / theta);
  return(0.1e1 / (t7 + 0.1e1));
}

// [[Rcpp::export]]
NumericVector cpp_qulogistic(const NumericVector x,
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
      out[i] = invcdf_ulogistic(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ulogistic(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood unit-logistic

// [[Rcpp::export]]
double cpp_loglikeulogistic(NumericVector x,
                           NumericVector lnx,
                           int n,
                           NumericVector mu,
                           NumericVector theta,
                           double tau)
{
  double ll = 0;

  for (int i = 0; i < n; i++)
    ll += logpdf_ulogistic(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// unit-logistic score vector

// [[Rcpp::export]]
NumericMatrix cpp_gradientulogistic(int n,
                                   NumericVector x,
                                   NumericMatrix U,
                                   NumericVector dmu_deta,
                                   NumericVector dtheta_dzeta,
                                   NumericVector mu,
                                   NumericVector theta,
                                   double tau)
{
  for (int i = 0; i < n; i++)
  {
    double t2 = -0.1e1 + mu[i];
    double t3 = 0.1e1 / mu[i];
    double t5 = t2 * t3;
    double t7 = pow(-t5, -theta[i]);
    double t8 = (-0.1e1 + tau) * t7;
    double t10 = mu[i] * mu[i];
    double t15 = theta[i] * (-t3 + t2 / t10);
    double t18 = 0.1e1 / t2 * mu[i];
    double t19 = 0.1e1 / tau;
    double t21 = 0.1e1 - x[i];
    double t23 = t21 / x[i];
    double t24 = pow(t23, theta[i]);
    double t25 = t19 * t24;
    double t29 = 0.1e1 / (0.1e1 - t8 * t25);
    double t37 = log(-t5);
    double t42 = log(t23);
    double t50 = log(t21);
    double t51 = log(x[i]);
    U(i,0) = (0.2000000000e1 * t8 * t15 * t18 * t25 * t29 + 0.1000000000e1 * t15 * t18) * dmu_deta[i];
    U(i,1) = (-0.2e1 * (t8 * t37 * t19 * t24 - t8 * t25 * t42) * t29 - 0.1000000000e1 * t37 + t50 - t51 + 0.1e1 / theta[i]) * dtheta_dzeta[i];
  }
  return(U);
}

// unit-logistic hessian matrix

// [[Rcpp::export]]
NumericMatrix cpp_hessianulogistic(int n,
                                  NumericVector x,
                                  NumericMatrix H,
                                  NumericVector mu,
                                  NumericVector theta,
                                  double tau)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = -0.1e1 + tau;
    double t2 = -0.1e1 + mu[i];
    double t3 = 0.1e1 / mu[i];
    double t4 = t2 * t3;
    double t5 = pow(-t4, -theta[i]);
    double t6 = t1 * t5;
    double t7 = theta[i] * theta[i];
    double t8 = mu[i] * mu[i];
    double t9 = 0.1e1 / t8;
    double t11 = -t3 + t2 * t9;
    double t12 = t11 * t11;
    double t13 = t7 * t12;
    double t15 = t2 * t2;
    double t16 = 0.1e1 / t15;
    double t17 = t16 * t8;
    double t18 = 0.1e1 / tau;
    double t21 = (1 - x[i]) / x[i];
    double t22 = pow((double) t21, theta[i]);
    double t23 = t18 * t22;
    double t25 = 0.1e1 - t6 * t23;
    double t26 = 0.1e1 / t25;
    double t27 = t23 * t26;
    double t35 = 0.2e1 * theta[i] * (t9 - t2 / t8 / mu[i]);
    double t37 = 0.1e1 / t2;
    double t38 = t37 * mu[i];
    double t39 = t38 * t27;
    double t42 = theta[i] * t11;
    double t43 = t6 * t42;
    double t44 = t16 * mu[i];
    double t49 = t22 * t26;
    double t53 = t1 * t1;
    double t54 = t5 * t5;
    double t57 = tau * tau;
    double t59 = t22 * t22;
    double t61 = t25 * t25;
    double t62 = 0.1e1 / t61;
    double t69 = t42 * t37;
    double t71 = log(-t4);
    double t73 = t71 * t18 * t22;
    double t79 = log((double) t21);
    double t80 = t23 * t79;
    double t88 = t6 * t73 - t6 * t80;
    double t92 = t11 * t37;
    double t93 = mu[i] * t18;
    double t97 = 0.2e1 * t88 * t62 * t6 * theta[i] * t92 * t93 * t22;
    double t98 = t92 * mu[i];
    double t116 = t71 * t71;
    double t123 = t79 * t79;
    double t129 = t88 * t88;
    H(i,0) = 0.2e1 * t6 * t13 * t17 * t27 + 0.2e1 * t6 * t35 * t39 - 0.2e1 * t43 * t44 * t27 + 0.2e1 * t43 * t37 * t18 * t49 + 0.2e1 * t53 * t54 * t13 * t17 / t57 * t59 * t62 + t35 * t38 - t42 * t44 + t69;
    H(i,1) = -0.2e1 * (t43 * t38 * t73 - t6 * t11 * t38 * t23 - t43 * t38 * t80) * t26 - t97 + t98;
    H(i,2) = -0.2e1 * (-t6 * t116 * t18 * t22 + 0.2e1 * t6 * t71 * t80 - t6 * t23 * t123) * t26 + 0.2e1 * t129 * t62 - 0.1e1 / t7;
  }
  return(H);
}

