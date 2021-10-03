#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// log-pdf ubs

inline double muinv(double mu, double theta, double tau)
{
  double rqnorm = R::qnorm(0.1e1 - tau, 0.0, 1.0, TRUE, FALSE);
  double t1 = log(mu);
  double t3 = rqnorm * rqnorm;
  double t5 = theta * theta;
  double t11 = sqrt(t3 * t5 + 0.4e1);
  return(-t1 - 0.5000000000e0 * t1 * t3 * t5 + 0.5000000000e0 * t1 * rqnorm * theta * t11);
}

inline double logpdf_ubs(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = muinv(mu, theta, tau);
  double t1 = log(theta);
  double t3 = log(alpha);
  double t5 = lnx;
  double t6 = 0.1e1 / t5;
  double t8 = alpha * t6;
  double t9 = sqrt(-t8);
  double t12 = log(t9 - t9 * t8);
  double t13 = theta * theta;
  double t14 = 0.1e1 / t13;
  return(-0.16120857137646180512e1 - t1 - t3 + t12 + 0.50000000000000000000e0 * t14 * t5 / alpha + 0.50000000000000000000e0 * t14 * alpha * t6 + t14 - t5);
}

// [[Rcpp::export]]
NumericVector cpp_dubs(const NumericVector x,
                        const NumericVector mu,
                        const NumericVector theta,
                        const double tau,
                        const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_ubs(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf ubs
inline double cdf_ubs(double x, double mu, double theta, double tau)
{
  double alpha = muinv(mu, theta, tau);
  double t1 = log(x);
  double t2 = -t1 / alpha;
  double t3 =  0.1e1 / t2;
  double t4 = (pow(t2, 0.05e1) - pow(t3, 0.05e1)) / theta;
  double rqnorm = R::pnorm(t4, 0.0, 1.0, TRUE, FALSE);
  return(0.1e1 - rqnorm);
}

// [[Rcpp::export]]
NumericVector cpp_pubs(const NumericVector x,
                       const NumericVector mu,
                       const NumericVector theta,
                       const double tau,
                       const bool lowertail = true,
                       const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_ubs(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf ubs
inline double invcdf_ubs(double x, double mu, double theta, double tau)
{
  double alpha = muinv(mu, theta, tau);
  double rqnorm = R::qnorm(0.1e1 - x, 0.0, 1.0, TRUE, FALSE);
  double t1 = theta * rqnorm;
  double t2 = pow(t1, 0.2e1);
  double t3 = 0.2e1 + t2 - t1 * pow(0.4e1 + t2, 0.05e1);
  return(exp(-0.2e1 * alpha / t3));
}

// [[Rcpp::export]]
NumericVector cpp_qubs(const NumericVector x,
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
      out[i] = invcdf_ubs(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_ubs(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }

  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood ubs

// [[Rcpp::export]]
double cpp_loglikeubs(NumericVector x, NumericVector lnx, int n, NumericVector mu, NumericVector theta, double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_ubs(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// score-vector ubs

// [[Rcpp::export]]
NumericMatrix cpp_gradientubs(int n,
                              NumericVector x,
                              NumericMatrix U,
                              NumericVector dmu_deta,
                              NumericVector dtheta_dzeta,
                              NumericVector mu,
                              NumericVector theta,
                              double tau)
{
  double rqnorm = R::qnorm(0.1e1 - tau, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t1 = 0.1e1 / mu[i];
    double t2 = log(mu[i]);
    double t3 = 0.1e1 / t2;
    double t5 = rqnorm * rqnorm;
    double t6 = theta[i] * theta[i];
    double t7 = t5 * t6;
    double t10 = sqrt(t7 + 0.4e1);
    double t12 = 0.2e1 + t7 - rqnorm * theta[i] * t10;
    double t14 = log(x[i]);
    double t15 = 0.1e1 / t14;
    double t16 = t2 * t12 * t15;
    double t17 = sqrt(t16);
    double t18 = 0.1e1 / t17;
    double t20 = t12 * t15;
    double t31 = 0.1e1 / (0.7071067812e0 * t17 + 0.3535533906e0 * t17 * t16);
    double t33 = 0.1e1 / t6;
    double t34 = t33 * t14;
    double t35 = t2 * t2;
    double t37 = 0.1e1 / t12;
    double t54 = 0.2e1 * t5 * theta[i] - rqnorm * t10 - t5 * rqnorm * t6 / t10;
    double t57 = t54 * t15;
    double t66 = 0.1e1 / t6 / theta[i];
    double t72 = t12 * t12;
    double t73 = -t1 * t3 + (0.3535533906e0 * t18 * t1 * t20 + 0.5303300859e0 * t17 * t1 * t20) * t31 + 0.1000000000e1 * t34 / t35 * t37 * t1 - 0.2500000000e0 * t33 * t1 * t20;
    double t74= -0.1e1 / theta[i] - t54 * t37 + (0.3535533906e0 * t18 * t2 * t57 + 0.5303300859e0 * t17 * t2 * t57) * t31 - 0.2000000000e1 * t66 + 0.2000000000e1 * t66 * t14 * t3 * t37 + 0.1000000000e1 * t34 * t3 / t72 * t54 + 0.5000000000e0 * t66 * t2 * t20 - 0.2500000000e0 * t33 * t2 * t57;
    U(i,0) = t73 * dmu_deta[i];
    U(i,1) = t74 * dtheta_dzeta[i];
  }
  return(U);
}

// hessian-matrix ubs

// [[Rcpp::export]]
NumericMatrix cpp_hessianubs(int n,
                             NumericVector x,
                             NumericMatrix H,
                             NumericVector mu,
                             NumericVector theta,
                             double tau)
{
  double rqnorm = R::qnorm(0.1e1 - tau, 0.0, 1.0, TRUE, FALSE);
  for(int i = 0; i < n; i++)
  {
    double t1 = mu[i] * mu[i];
    double t2 = 0.1e1 / t1;
    double t3 = log(mu[i]);
    double t4 = 0.1e1 / t3;
    double t6 = t3 * t3;
    double t7 = 0.1e1 / t6;
    double t9 = rqnorm * rqnorm;
    double t10 = theta[i] * theta[i];
    double t11 = t9 * t10;
    double t13 = t11 + 0.4e1;
    double t14 = sqrt(t13);
    double t16 = 0.2e1 + t11 - rqnorm * theta[i] * t14;
    double t18 = log(x[i]);
    double t19 = 0.1e1 / t18;
    double t20 = t3 * t16 * t19;
    double t21 = sqrt(t20);
    double t22 = t21 * t20;
    double t23 = 0.1e1 / t22;
    double t25 = t16 * t16;
    double t26 = t18 * t18;
    double t27 = 0.1e1 / t26;
    double t28 = t25 * t27;
    double t31 = 0.1e1 / t21;
    double t32 = t31 * t2;
    double t33 = t16 * t19;
    double t44 = 0.7071067812e0 * t21 + 0.3535533906e0 * t22;
    double t45 = 0.1e1 / t44;
    double t47 = 0.1e1 / mu[i];
    double t48 = t31 * t47;
    double t51 = t21 * t47;
    double t54 = 0.3535533906e0 * t48 * t33 + 0.5303300859e0 * t51 * t33;
    double t55 = t54 * t54;
    double t56 = t44 * t44;
    double t57 = 0.1e1 / t56;
    double t59 = 0.1e1 / t10;
    double t60 = t59 * t18;
    double t63 = 0.1e1 / t16;
    double t68 = t7 * t63;
    double t80 = t9 * rqnorm;
    double t82 = 0.1e1 / t14;
    double t84 = 0.2e1 * t9 * theta[i] - rqnorm * t14 - t80 * t10 * t82;
    double t87 = t27 * t47 * t16;
    double t90 = t84 * t19;
    double t93 = t31 * t3;
    double t103 = t21 * t3;
    double t106 = 0.3535533906e0 * t93 * t90 + 0.5303300859e0 * t103 * t90;
    double t109 = t10 * theta[i];
    double t110 = 0.1e1 / t109;
    double t111 = t110 * t18;
    double t116 = 0.1e1 / t25;
    double t127 = (-0.1767766953e0 * t23 * t3 * t84 * t87 + 0.3535533906e0 * t48 * t90 + 0.2651650430e0 * t93 * t84 * t87 + 0.5303300859e0 * t51 * t90) * t45 - t106 * t57 * t54 - 0.2000000000e1 * t111 * t68 * t47 - 0.1000000000e1 * t60 * t7 * t116 * t84 * t47 + 0.5000000000e0 * t110 * t47 * t33 - 0.2500000000e0 * t59 * t47 * t90;
    double t132 = t9 * t9;
    double t138 = 0.2e1 * t9 - 0.3e1 * t80 * t82 * theta[i] + t132 * rqnorm * t109 / t14 / t13;
    double t140 = t84 * t84;
    double t143 = t140 * t27;
    double t146 = t138 * t19;
    double t156 = t106 * t106;
    double t158 = t10 * t10;
    double t159 = 0.1e1 / t158;
    double t165 = t4 * t116;
    double t187 = t59 - t138 * t63 + t140 * t116 + (-0.1767766953e0 * t23 * t6 * t143 + 0.3535533906e0 * t93 * t146 + 0.2651650430e0 * t31 * t6 * t143 + 0.5303300859e0 * t103 * t146) * t45 - t156 * t57 + 0.6000000000e1 * t159 - 0.6000000000e1 * t159 * t18 * t4 * t63 - 0.4000000000e1 * t111 * t165 * t84 - 0.2000000000e1 * t60 * t4 / t25 / t16 * t140 + 0.1000000000e1 * t60 * t165 * t138 - 0.1500000000e1 * t159 * t3 * t33 + 0.1000000000e1 * t110 * t3 * t90 - 0.2500000000e0 * t59 * t3 * t146;
    H(i,0) = t2 * t4 + t2 * t7 + (-0.1767766953e0 * t23 * t2 * t28 - 0.3535533906e0 * t32 * t33 + 0.2651650430e0 * t32 * t28 - 0.5303300859e0 * t21 * t2 * t33) * t45 - t55 * t57 - 0.2000000000e1 * t60 / t6 / t3 * t63 * t2 - 0.1000000000e1 * t60 * t68 * t2 + 0.2500000000e0 * t59 * t2 * t33;
    H(i,1) = t127;
    H(i,2) = t187;
  }
  return(H);
}
