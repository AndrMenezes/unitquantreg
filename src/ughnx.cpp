#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf ughnx

inline double logpdf_ughnx(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -mu / pow(R::qnorm(tau/ 0.2e1 + 0.1e1 / 0.2e1, 0.0, 1.0, TRUE, FALSE), 0.1e1 / theta) / (-0.1e1 + mu);
  double t1 = log(theta);
  double t2 = lnx;
  double t5 = 0.1e1 - x;
  double t6 = log(t5);
  double t7 = t6;
  double t8 = log(alpha);
  double t12 = pow(x, theta);
  double t13 = t12 * t12;
  double t14 = theta;
  double t15 = pow(t5, -t14);
  double t16 = t15 * t15;
  double t18 = pow(alpha, -t14);
  double t19 = t18 * t18;
  return(-0.22579135264472743239e0 + t1 - t2 - t7 + theta * (t2 - t7 - t8) - 0.5e0 * t13 * t16 * t19);
}

// [[Rcpp::export]]
NumericVector cpp_dughnx(const NumericVector x,
                         const NumericVector mu,
                         const NumericVector theta,
                         const double tau,
                         const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_ughnx(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf ughnx
inline double cdf_ughnx(double x, double mu, double theta, double tau)
{
  double alpha = -mu / pow(R::qnorm(tau/ 0.2e1 + 0.1e1 / 0.2e1, 0.0, 1.0, TRUE, FALSE), 0.1e1 / theta) / (-0.1e1 + mu);
  double t6 = pow(x / (0.1e1 - x) / alpha, theta);
  double t7 = R::pnorm(t6, 0.0, 1.0, TRUE, FALSE);
  return(0.2e1 * t7 - 0.1e1);
}

// [[Rcpp::export]]
NumericVector cpp_pughnx(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_ughnx(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf ughnx
inline double invcdf_ughnx(double x, double mu, double theta, double tau)
{
  double alpha = -mu / pow(R::qnorm(tau/ 0.2e1 + 0.1e1 / 0.2e1, 0.0, 1.0, TRUE, FALSE), 0.1e1 / theta) / (-0.1e1 + mu);
  double t2 = R::qnorm(x / 0.2e1 + 0.1e1 / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  double t4 = pow(t2, 0.1e1 / theta);
  double t5 = alpha * t4;
  return(t5 / (t5 + 0.1e1));
}

// [[Rcpp::export]]
NumericVector cpp_qughnx(const NumericVector x,
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
      out[i] = invcdf_ughnx(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ughnx(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood ughnx

// [[Rcpp::export]]
double cpp_loglikeughnx(NumericVector x,
                            NumericVector lnx,
                            int n,
                            NumericVector mu,
                            NumericVector theta,
                            double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_ughnx(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// score-vector ughnx

// [[Rcpp::export]]
NumericMatrix cpp_gradientughnx(int n,
                                NumericVector x,
                                NumericMatrix U,
                                NumericVector dmu_deta,
                                NumericVector dtheta_dzeta,
                                NumericVector mu,
                                NumericVector theta,
                                double tau)
{
  double t2 = R::qnorm(tau / 0.2e1 + 0.1e1 / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t3 = 0.1e1 / theta[i];
    double t4 = pow(t2, t3);
    double t5 = 0.1e1 / t4;
    double t6 = -0.1e1 + mu[i];
    double t7 = 0.1e1 / t6;
    double t10 = mu[i] * t5;
    double t11 = t6 * t6;
    double t15 = -t5 * t7 + t10 / t11;
    double t17 = 0.1e1 / mu[i];
    double t22 = pow(x[i], theta[i]);
    double t23 = t22 * t22;
    double t25 = 0.1e1 - x[i];
    double t26 = theta[i];
    double t27 = pow(t25, -t26);
    double t28 = t27 * t27;
    double t29 = t23 * t28;
    double t31 = t10 * t7;
    double t32 = pow(-t31, -t26);
    double t33 = t32 * t32;
    double t42 = log(x[i]);
    double t43 = log(t25);
    double t45 = log(-t31);
    double t46 = t45;
    double t47 = log(t2);
    double t49 = 0.1000000000e1 * t3 * t47;
    U(i,0) = (theta[i] * t15 * t17 * t4 * t6 - 0.1000000000e1 * t29 * t33 * theta[i] * t15 * t17 * t4 * t6) * dmu_deta[i];
    U(i,1) = (t3 + t42 - t43 - t46 - t49 - 0.10e1 * t29 * t33 * t42 + 0.10e1 * t29 * t33 * t43 - 0.10e1 * t29 * t33 * (-t46 - t49)) * dtheta_dzeta[i];
  }
  return(U);
}

// hessian-matrix ughnx

// [[Rcpp::export]]
NumericMatrix cpp_hessianughnx(int n,
                               NumericVector x,
                               NumericMatrix H,
                               NumericVector mu,
                               NumericVector theta,
                               double tau)
{
  double t2 = R::qnorm(tau / 0.2e1 + 0.1e1 / 0.2e1, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t3 = 0.1e1 / theta[i];
    double t4 = pow(t2, t3);
    double t5 = 0.1e1 / t4;
    double t6 = -0.1e1 + mu[i];
    double t7 = t6 * t6;
    double t8 = 0.1e1 / t7;
    double t11 = mu[i] * t5;
    double t16 = 0.2e1 * t5 * t8 - 0.2e1 * t11 / t7 / t6;
    double t18 = 0.1e1 / mu[i];
    double t19 = t18 * t4;
    double t20 = t19 * t6;
    double t23 = 0.1e1 / t6;
    double t24 = t5 * t23;
    double t28 = -t24 + t11 * t8;
    double t29 = theta[i] * t28;
    double t30 = mu[i] * mu[i];
    double t31 = 0.1e1 / t30;
    double t36 = t29 * t19;
    double t38 = pow(x[i], theta[i]);
    double t39 = t38 * t38;
    double t41 = 0.1e1 - x[i];
    double t42 = theta[i];
    double t43 = pow(t41, -t42);
    double t44 = t43 * t43;
    double t45 = t39 * t44;
    double t47 = t11 * t23;
    double t48 = pow(-t47, -t42);
    double t49 = t48 * t48;
    double t50 = theta[i] * theta[i];
    double t53 = t28 * t28;
    double t55 = t4 * t4;
    double t61 = t45 * t49 * theta[i];
    double t63 = t4 * t6;
    double t71 = t45 * t49;
    double t75 = t28 * t18;
    double t76 = t75 * t63;
    double t77 = t76;
    double t78 = log(x[i]);
    double t79 = t49 * t78;
    double t81 = t29 * t20;
    double t83 = 0.2000000000e1 * t45 * t79 * t81;
    double t84 = log(t41);
    double t85 = t49 * t84;
    double t88 = 0.2000000000e1 * t45 * t85 * t81;
    double t89 = log(-t47);
    double t91 = log(t2);
    double t94 = -t89 - 0.1000000000e1 * t3 * t91;
    double t98 = 0.2000000000e1 * t45 * t49 * t94 * t81;
    double t99 = t71 * t76;
    double t102 = 0.1e1 / t50;
    double t110 = -t24 * t102 * t91 + t11 * t8 * t102 * t91;
    double t117 = t4 * t91 * t6;
    double t131 = t78 * t78;
    double t141 = t84 * t84;
    double t148 = t94 * t94;
    H(i,0) = theta[i] * t16 * t20 - t29 * t31 * t4 * t6 + t36 - 0.2000000000e1 * t45 * t49 * t50 * t53 * t31 * t55 * t7 - 0.1000000000e1 * t61 * t16 * t18 * t63 + 0.1000000000e1 * t61 * t28 * t31 * t63 - 0.1000000000e1 * t71 * t36;
    H(i,1) = t77 + theta[i] * t110 * t20 - t3 * t28 * t18 * t117 - t83 + t88 - t98 - 0.1000000000e1 * t99 - 0.1000000000e1 * t61 * t110 * t18 * t63 + 0.1000000000e1 * t45 * t49 * t3 * t75 * t117;
    H(i,2) = -t102 - 0.20e1 * t45 * t49 * t131 + 0.40e1 * t45 * t79 * t84 - 0.40e1 * t45 * t79 * t94 - 0.20e1 * t45 * t49 * t141 + 0.40e1 * t45 * t85 * t94 - 0.20e1 * t45 * t49 * t148;
  }
  return(H);
}
