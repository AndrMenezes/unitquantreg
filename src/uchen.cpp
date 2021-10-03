#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// log-pdf uchen

inline double logpdf_uchen(double x, double lnx, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (exp(pow(-log(mu), theta)) - 0.1e1);
  double t1 = log(alpha);
  double t2 = lnx;
  double t3 = log(-t2);
  double t5 = log(theta);
  double t6 = pow(-t2, theta);
  double t7 = exp(t6);
  return(t1 + theta * t3 - t3 + t5 + t6 + alpha - alpha * t7 - t2);
  }

// [[Rcpp::export]]
NumericVector cpp_duchen(const NumericVector x,
                         const NumericVector mu,
                         const NumericVector theta,
                         const double tau,
                         const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = logpdf_uchen(x[i], log(x[i]), mu[i % nmu], theta[i % nth], tau);

  if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf uchen
inline double cdf_uchen(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (exp(pow(-log(mu), theta)) - 0.1e1);
  double t1 = log(x);
  double t2 = pow(-t1, theta);
  double t3 = exp(t2);
  return(exp(alpha * (0.1e1 - t3)));
}

// [[Rcpp::export]]
NumericVector cpp_puchen(const NumericVector x,
                             const NumericVector mu,
                             const NumericVector theta,
                             const double tau,
                             const bool lowertail = true,
                             const bool logprob = false)
{
  const int n = x.length(); NumericVector out(n);
  const int nmu = mu.length(); const int nth = theta.length();

  for(int i = 0; i < n; i++)
    out[i] = cdf_uchen(x[i], mu[i % nmu], theta[i % nth], tau);

  if (!lowertail) out = 0.1e1 - out;
  if (logprob) out = Rcpp::log(out);

  return(out);
}

// inv-cdf uchen
inline double invcdf_uchen(double x, double mu, double theta, double tau)
{
  double alpha = -log(tau) / (exp(pow(-log(mu), theta)) - 0.1e1);
  double t1 = log(alpha);
  double t2 = log(x);
  double t4 = log(alpha - t2);
  double t7 = pow(-t1 + t4, 0.1e1 / theta);
  double t8 = exp(t7);
  return(0.1e1 / t8);
}

// [[Rcpp::export]]
NumericVector cpp_quchen(const NumericVector x,
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
      out[i] = invcdf_uchen(x[i], mu[i % nmu], theta[i % nth], tau);
  }
  else
  {
    for(int i = 0; i < n; i++)
      out[i] = invcdf_uchen(0.1e1 - x[i], mu[i % nmu], theta[i % nth], tau);
  }

  if(logprob) return(Rcpp::log(out)); else return(out);
}

// log-likelihood uchen

// [[Rcpp::export]]
double cpp_loglikeuchen(NumericVector x, NumericVector lnx, int n, NumericVector mu, NumericVector theta, double tau)
{
  double ll = 0;

  for(int i = 0; i < n; i++)
    ll += logpdf_uchen(x[i], lnx[i], mu[i], theta[i], tau);

  return(-ll);
}

// score vector uchen

// [[Rcpp::export]]
NumericMatrix cpp_gradientuchen(int n,
                                NumericVector x,
                                NumericMatrix U,
                                NumericVector dmu_deta,
                                NumericVector dtheta_dzeta,
                                NumericVector mu,
                                NumericVector theta,
                                double tau)
{
  double t14 = log(tau);
  for(int i = 0; i < n; i++)
  {
    double t1 = log(mu[i]);
    double t2 = pow(-t1, theta[i]);
    double t3 = exp(t2);
    double t4 = t3 - 0.1e1;
    double t5 = 0.1e1 / t4;
    double t6 = t5 * t2;
    double t7 = t6 * theta[i];
    double t8 = 0.1e1 / mu[i];
    double t9 = 0.1e1 / t1;
    double t10 = t8 * t9;
    double t11 = t3 * t3;
    double t16 = log(-t14 * t5);
    double t17 = t16 * t2;
    double t19 = t10 * t3;
    double t22 = theta[i] * theta[i];
    double t23 = log(x[i]);
    double t24 = log(-t23);
    double t28 = t24 * t2;
    double t31 = log(theta[i]);
    double t32 = t31 * t2;
    double t35 = pow(-t23, theta[i]);
    double t36 = t35 * t2;
    double t39 = t23 * t2;
    double t45 = theta[i] * t24;
    double t47 = t24 * t3;
    double t50 = exp(t35);
    double t53 = t16 * t3 - t16 + t45 * t3 - t45 - t47 + t24 + t31 * t3 - t31 + t35 * t3 - t35 - t14 + t14 * t50 - t23 * t3 + t23;
    double t54 = t4 * t4;
    double t56 = t53 / t54;
    double t63 = log(-t1);
    double t66 = t63 * t3;
    double t70 = t2 * t63 * t3;
    double t73 = 0.1e1 / theta[i];
    double t76 = t35 * t24;
    double t83 = -t6 * t63 * t11 + t17 * t66 + t6 * t66 + t47 + t45 * t70 - t24 - t28 * t66 + t73 * t3 + t32 * t66 - t73 + t76 * t3 + t36 * t66 - t76 + t14 * t35 * t24 * t50 - t39 * t66;
    double t84 = (-t7 * t10 * t11 + t17 * theta[i] * t19 + t7 * t19 + t22 * t24 * t2 * t19 - t28 * theta[i] * t19 + t32 * theta[i] * t19 + t36 * theta[i] * t19 - t39 * theta[i] * t19) * t5 - t56 * t2 * theta[i] * t8 * t9 * t3;
    double t85 = t83 * t5 - t56 * t70;
    U(i,0) = t84 * dmu_deta[i];
    U(i,1) = t85 * dtheta_dzeta[i];
  }
  return(U);
}

// hessian matrix uchen

// [[Rcpp::export]]
NumericMatrix cpp_hessianuchen(int n,
                               NumericVector x,
                               NumericMatrix H,
                               NumericVector mu,
                               NumericVector theta,
                               double tau)
{
  double t50 = log(tau);
  for(int i = 0; i < n; i++)
  {
    double t1 = log(mu[i]);
    double t2 = pow(-t1, theta[i]);
    double t3 = exp(t2);
    double t4 = t3 - 0.1e1;
    double t5 = 0.1e1 / t4;
    double t6 = t5 * t2;
    double t7 = t6 * theta[i];
    double t8 = mu[i] * mu[i];
    double t9 = 0.1e1 / t8;
    double t10 = 0.1e1 / t1;
    double t11 = t9 * t10;
    double t12 = t11 * t3;
    double t14 = t1 * t1;
    double t15 = 0.1e1 / t14;
    double t16 = t9 * t15;
    double t17 = t16 * t3;
    double t19 = theta[i] * theta[i];
    double t20 = log(x[i]);
    double t21 = log(-t20);
    double t23 = t19 * t21 * t2;
    double t27 = t21 * t2;
    double t28 = t27 * theta[i];
    double t31 = log(theta[i]);
    double t32 = t31 * t2;
    double t33 = t32 * theta[i];
    double t36 = pow(-t20, theta[i]);
    double t37 = t36 * t2;
    double t38 = t37 * theta[i];
    double t41 = t20 * t2;
    double t42 = t41 * theta[i];
    double t45 = t3 * t3;
    double t48 = t16 * t45;
    double t52 = log(-t50 * t5);
    double t53 = t52 * t2;
    double t54 = t53 * theta[i];
    double t56 = t2 * t2;
    double t57 = t52 * t56;
    double t60 = -t7 * t12 - t7 * t17 - t23 * t12 - 0.2e1 * t23 * t17 + t28 * t12 + t28 * t17 - t33 * t12 - t33 * t17 - t38 * t12 - t38 * t17 + t42 * t12 + t42 * t17 + t7 * t11 * t45 + t7 * t48 - t54 * t12 + t57 * t19 * t17;
    double t61 = t5 * t56;
    double t62 = t61 * t19;
    double t65 = t19 * theta[i] * t21;
    double t68 = t21 * t56;
    double t69 = t68 * t19;
    double t71 = t31 * t56;
    double t74 = t36 * t56;
    double t77 = t20 * t56;
    double t80 = t4 * t4;
    double t81 = 0.1e1 / t80;
    double t82 = t81 * t56;
    double t83 = t82 * t19;
    double t84 = t45 * t3;
    double t90 = t6 * t19;
    double t104 = t62 * t17 + t17 * t65 * t56 - t69 * t17 + t71 * t19 * t17 + t74 * t19 * t17 - t77 * t19 * t17 + t83 * t16 * t84 - 0.3e1 * t62 * t48 - t83 * t48 - t90 * t48 + t53 * t19 * t17 + t90 * t17 + t65 * t2 * t17 + t32 * t19 * t17 + t37 * t19 * t17 - t41 * t19 * t17 - t54 * t17;
    double t107 = 0.1e1 / mu[i];
    double t108 = t107 * t10;
    double t109 = t108 * t45;
    double t111 = t108 * t3;
    double t115 = t28 * t111;
    double t120 = (-t7 * t109 + t54 * t111 + t7 * t111 + t23 * t111 - t115 + t33 * t111 + t38 * t111 - t42 * t111) * t81;
    double t123 = t10 * t3;
    double t124 = theta[i] * t107 * t123;
    double t128 = theta[i] * t21;
    double t130 = t21 * t3;
    double t133 = exp(t36);
    double t136 = t52 * t3 - t52 + t128 * t3 - t128 - t130 + t21 + t31 * t3 - t31 + t36 * t3 - t36 - t50 + t50 * t133 - t20 * t3 + t20;
    double t139 = t136 / t80 / t4;
    double t141 = t19 * t9;
    double t146 = t136 * t81;
    double t147 = t146 * t2;
    double t148 = t15 * t3;
    double t149 = t141 * t148;
    double t151 = theta[i] * t9;
    double t159 = log(-t1);
    double t160 = t159 * t3;
    double t161 = t108 * t160;
    double t163 = t36 * t21;
    double t168 = t82 * t159;
    double t172 = t61 * t159;
    double t174 = t45 * theta[i] * t108;
    double t189 = t33 * t161 + t163 * t2 * t124 + t38 * t161 - t42 * t161 + t168 * t84 * theta[i] * t108 - 0.3e1 * t172 * t174 - t168 * t174 + t69 * t161 + t172 * t124 + t23 * t161 - t68 * t159 * t124 + t71 * t159 * t124 + t74 * t159 * t124 - t77 * t159 * t124;
    double t192 = t159 * t45;
    double t208 = t57 * t159 * t124 - t7 * t108 * t192 + t54 * t161 + t7 * t161 - t28 * t161 - t6 * t109 + t53 * t111 + t6 * t111 - t27 * t111 + t32 * t111 + t37 * t111 - t41 * t111 + t2 * t107 * t123 + 0.2e1 * t115;
    double t215 = t2 * t159 * t3;
    double t217 = t27 * t160;
    double t218 = 0.1e1 / theta[i];
    double t223 = t50 * t36;
    double t227 = -t6 * t192 + t53 * t160 + t6 * t160 + t130 + t128 * t215 - t21 - t217 + t218 * t3 + t32 * t160 - t218 + t163 * t3 + t37 * t160 - t163 + t223 * t21 * t133 - t41 * t160;
    double t228 = t227 * t81;
    double t232 = t56 * t159;
    double t242 = (t189 + t208) * t5 - t228 * t2 * t124 - t120 * t215 + 0.2e1 * t139 * t232 * t174 - t146 * t2 * theta[i] * t161 - t147 * t111 - t146 * t232 * t124;
    double t243 = 0.1e1 / t19;
    double t244 = t21 * t21;
    double t245 = t36 * t244;
    double t248 = t159 * t159;
    double t249 = t56 * t248;
    double t250 = t249 * t3;
    double t255 = t2 * t248 * t3;
    double t260 = t248 * t45;
    double t264 = t248 * t3;
    double t267 = t243 + t245 * t3 - t245 - t243 * t3 + t128 * t250 + 0.2e1 * t163 * t215 + t128 * t255 + 0.2e1 * t217 + t82 * t248 * t84 - 0.3e1 * t61 * t260 - t82 * t260 - t27 * t264 + t32 * t264;
    double t269 = t244 * t133;
    double t283 = t36 * t36;
    double t287 = t37 * t264 + t223 * t269 - t41 * t264 - t6 * t260 + t53 * t264 + t6 * t264 + 0.2e1 * t218 * t2 * t160 + t57 * t264 + t61 * t264 - t68 * t264 + t71 * t264 + t74 * t264 + t50 * t283 * t269 - t77 * t264;
    H(i,0) = (t60 + t104) * t5 - 0.2e1 * t120 * t2 * t124 + 0.2e1 * t139 * t56 * t141 * t15 * t45 - t147 * t149 + t147 * t151 * t123 + t147 * t151 * t148 - t146 * t56 * t149;
    H(i,1) = t242;
    H(i,2)= (t267 + t287) * t5 - 0.2e1 * t228 * t215 + 0.2e1 * t139 * t249 * t45 - t146 * t255 - t146 * t250;
  }
  return(H);
}

