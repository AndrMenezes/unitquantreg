#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#define GETV(x, i) x[i % x.length()]

// log-pdf burr-xii
inline double logpdf_uburrxii(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -log(tau) / log(pow(-log(mu), theta) + 0.1e1);
  double t1 = lnx;
  double t2 = pow(-t1, theta);
  double t4 = log(0.1e1 + t2);
  double t6 = log(alpha);
  double t7 = log(-t1);
  double t9 = log(theta);
  double t11 = log(-0.1e1 / t1);
  return(-alpha * t4 + t6 + theta * t7 + t9 - t1 + t11 - t4);
}

// [[Rcpp::export]]
NumericVector cpp_duburrxii(const NumericVector x,
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
    out[i] = logpdf_uburrxii(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out);
  else return(Rcpp::exp(out));
}

// cdf burr-xii
inline double cdf_uburrxii(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / log(pow(-log(mu), theta) + 0.1e1);
  double t1 = log(x);
  double t2 = pow(-t1, theta);
  return(pow(0.1e1 + t2, -alpha));

}

// [[Rcpp::export]]
NumericVector cpp_puburrxii(const NumericVector x,
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
    out[i] = cdf_uburrxii(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf burr-xii
inline double invcdf_uburrxii(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / log(pow(-log(mu), theta) + 0.1e1);
  double t2 = pow(x, -0.1e1 / alpha);
  double t5 = pow(t2 - 0.1e1, 0.1e1 / theta);
  return(exp(-t5));
}

// [[Rcpp::export]]
NumericVector cpp_quburrxii(const NumericVector x,
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
      out[i] = invcdf_uburrxii(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_uburrxii(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }
  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood burr-xii

// [[Rcpp::export]]
double cpp_loglikeuburrxii(NumericVector x,
                           NumericVector lnx,
                           int n,
                           NumericVector mu,
                           NumericVector theta,
                           double tau)
{
  double ll = 0;

  for (int i = 0; i < n; i++)
    ll += logpdf_uburrxii(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// burr-xii score vector

// [[Rcpp::export]]
NumericMatrix cpp_gradientuburrxii(int n,
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
    double t2 = log(mu[i]);
    double t3 = log(-t2);
    double t5 = exp(t3 * theta[i]);
    double t6 = t5 + 0.1e1;
    double t7 = log(t6);
    double t8 = t7 * t7;
    double t10 = t1 / t8;
    double t11 = log(x[i]);
    double t12 = pow(-t11, theta[i]);
    double t13 = 0.1e1 + t12;
    double t14 = log(t13);
    double t15 = 0.1e1 / mu[i];
    double t18 = 0.1e1 / t2;
    double t20 = 0.1e1 / t6;
    double t21 = t5 * t20;
    double t24 = 0.1e1 / t7;
    double t36 = log(-t11);
    double t39 = t12 * t36 / t13;
    double t40 = -t10 * t14 * t15 * t18 * theta[i] * t21 - t24 * t15 * t18 * theta[i] * t5 * t20;
    double t41 = -t10 * t14 * t3 * t5 * t20 + t1 * t24 * t39 - t24 * t3 * t21 + t36 + 0.1e1 / theta[i] - t39;
    U(i,0) = t40 * dmu_deta[i];
    U(i,1) = t41 * dtheta_dzeta[i];
  }
  return(U);
}

// burr-xii hessian matrix

// [[Rcpp::export]]
NumericMatrix cpp_hessianuburrxii(int n,
                                  NumericVector x,
                                  NumericMatrix H,
                                  NumericVector mu,
                                  NumericVector theta,
                                  double tau)
{
  double t1 = log(tau);
  for (int i = 0; i < n; i++)
  {
    double t1 = log(tau);
    double t2 = log(mu[i]);
    double t3 = log(-t2);
    double t5 = exp(t3 * theta[i]);
    double t6 = t5 + 0.1e1;
    double t7 = log(t6);
    double t8 = t7 * t7;
    double t11 = t1 / t8 / t7;
    double t12 = log(x[i]);
    double t13 = pow(-t12, theta[i]);
    double t14 = 0.1e1 + t13;
    double t15 = log(t14);
    double t16 = mu[i] * mu[i];
    double t17 = 0.1e1 / t16;
    double t18 = t15 * t17;
    double t20 = t2 * t2;
    double t21 = 0.1e1 / t20;
    double t22 = theta[i] * theta[i];
    double t23 = t21 * t22;
    double t24 = t5 * t5;
    double t25 = t6 * t6;
    double t26 = 0.1e1 / t25;
    double t27 = t24 * t26;
    double t28 = t23 * t27;
    double t31 = 0.1e1 / t8;
    double t32 = t1 * t31;
    double t33 = t32 * t18;
    double t34 = 0.1e1 / t2;
    double t35 = t34 * theta[i];
    double t36 = 0.1e1 / t6;
    double t37 = t5 * t36;
    double t38 = t35 * t37;
    double t49 = t22 * t24 * t26;
    double t51 = 0.1e1 / t7;
    double t52 = t51 * t17;
    double t55 = theta[i] * t5 * t36;
    double t57 = t52 * t21;
    double t64 = t15 * t3;
    double t66 = 0.1e1 / mu[i];
    double t67 = t66 * t34;
    double t69 = t27 * t67 * theta[i];
    double t72 = t32 * t15;
    double t75 = t32 * t64;
    double t76 = t67 * t55;
    double t79 = log(-t12);
    double t80 = t13 * t79;
    double t81 = 0.1e1 / t14;
    double t88 = t26 * t66 * t35;
    double t94 = t51 * t3;
    double t99 = 0.2e1 * t11 * t64 * t69 - t72 * t67 * t37 - t75 * t76 + t75 * t69 - t32 * t80 * t81 * t76 + t31 * t3 * t24 * t88 - t51 * t66 * t34 * t5 * t36 - t94 * t66 * t38 + t94 * t24 * t88;
    double t101 = t3 * t3;
    double t103 = t101 * t24 * t26;
    double t115 = t1 * t51;
    double t116 = t79 * t79;
    double t118 = t13 * t116 * t81;
    double t120 = t13 * t13;
    double t122 = t14 * t14;
    double t124 = t120 * t116 / t122;
    double t128 = t51 * t101;
    double t132 = 0.2e1 * t11 * t15 * t103 - 0.2e1 * t32 * t80 * t81 * t3 * t37 - t72 * t101 * t5 * t36 + t72 * t103 + t115 * t118 - t115 * t124 + t31 * t101 * t27 - t128 * t37 + t128 * t27 - 0.1e1 / t22 - t118 + t124;
    H(i,0) = 0.2e1 * t11 * t18 * t28 + t33 * t38 + t33 * t21 * theta[i] * t37 - t33 * t23 * t37 + t33 * t28 + t31 * t17 * t21 * t49 + t52 * t34 * t55 + t57 * t55 - t57 * t22 * t5 * t36 + t57 * t49;
    H(i,1) = t99;
    H(i,2) = t132;
  }
  return(H);
}

