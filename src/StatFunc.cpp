/*
 * Implementations of the statistical functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "StatFunc.hpp"

#include <algorithm>
#include <cmath>

#include "CommFunc.hpp"

extern "C" {
#include "cdflib.h"
}

////////// P-value Calculatiion Functions Start ////////////////

double StatFunc::gammln(double xx) {
  int j;
  double x, y, tmp, ser;
  static const double cof[6] = {76.18009172947146,  -86.50532032941677,    24.01409824083091,
                                -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j < 6; j++) ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}

double StatFunc::chi_prob(double df, double chi_sqr_val) { return 1 - StatFunc::gammp(df * 0.5, chi_sqr_val * 0.5); }

double StatFunc::gammp(const double a, const double x) {
  double gamser, gammcf, gln;

  if (x < 0.0 || a <= 0.0) throw("Invalid arguments in routine gammp");

  if (x < a + 1.0) {
    gser(gamser, a, x, gln);
    return gamser;
  } else {
    gcf(gammcf, a, x, gln);
    return 1.0 - gammcf;
  }
}

void StatFunc::gser(double& gamser, const double a, const double x, double& gln) {
  const int ITMAX = 500;
  const double Eps = 1.0e-08;
  int n;
  double sum, del, ap;

  gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) throw("x less than 0 in routine gser");
    gamser = 0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 0; n < ITMAX; n++) {
      ++ap;
      del *= x / ap;
      sum += del;
      if (CommFunc::Abs(del) < CommFunc::Abs(sum) * Eps) {
        gamser = sum * exp(-x + a * log(x) - gln);
        return;
      }
    }
    throw("a too large, ITMAX too small in routine gser");
    return;
  }
}

void StatFunc::gcf(double& gammcf, const double a, const double x, double& gln) {
  const int ITMAX = 100;
  const double EPS = std::numeric_limits<double>::epsilon();
  const double FPMIN = (std::numeric_limits<double>::min)() / EPS;
  int i;
  double an, b, c, d, del, h;

  gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= ITMAX; i++) {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (CommFunc::Abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (CommFunc::Abs(del - 1.0) <= EPS) break;
  }
  if (i > ITMAX) throw("a too large, ITMAX too small in gcf");
  gammcf = exp(-x + a * log(x) - gln) * h;
}
////////// P-value Calculatiion Functions End ////////////////

// not good, do not use
double StatFunc::chi_val(double df, double prob) {
  double walk = 100.0;
  double chi_val = walk;
  double way = 0.0, preway = 0.0;
  double prob_buf = chi_prob(df, chi_val);

  if (CommFunc::Abs(prob_buf - prob) < 1.0e-08) return chi_val;

  if (prob_buf > prob) {
    preway = way = 1.0;
    chi_val += walk;
  } else {
    preway = way = -1.0;
    chi_val -= walk;
  }

  while (true) {
    prob_buf = chi_prob(df, chi_val);

    if (prob_buf > prob) way = 1.0;
    else way = -1.0;

    if (CommFunc::Abs(preway - way) > 1.0e-08) walk *= 0.5;
    chi_val += walk * way;
    preway = way;

    if (walk < 1.0e-08) break;
  }

  return chi_val;
}

// Default: upper-tail
double StatFunc::pnorm(double x) {
  double z = 0.0;
  if (x > 0) z = -1.0 * x;
  else z = x;
  double sqrt2pi = 2.50662827463;
  double t0, z1, p0;
  t0 = 1 / (1 + 0.2316419 * fabs(z));
  z1 = exp(-0.5 * z * z) / sqrt2pi;
  p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
  return x >= 0 ? p0 : 1.0 - p0;
}

double StatFunc::pchisq(double x, double df) {
  if (x < 0) return -9;

  double p, q;
  int st = 0;      // error variable
  int w = 1;       // function variable
  double bnd = 1;  // boundary function

  // NCP is set to 0
  cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

  // Check status
  if (st != 0) return -9;

  // Return p-value
  return q;
}

// q is not good to be less than 1e-161
double StatFunc::qchisq(double q, double df) {
  if (q < 0) return -9;
  else if (q >= 1) return 0;

  double x;
  double p = 1 - q;
  int st = 0;      // error variable
  int w = 2;       // function variable
  double bnd = 1;  // boundary function

  // NCP is set to 0
  cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

  // Check status
  if (st != 0) return -9;
  if (q < 1e-161) {
    double tmp = pchisq(x, 1);
    while (q / tmp > 10) {
      x -= 10;
      tmp = pchisq(x, 1);
    }
    if (q / tmp > 1) {
      while (q / tmp > 1.01) {
        x -= 0.1;
        tmp = pchisq(x, 1);
      }
    } else if (q / tmp < 1) {
      while (q / tmp < 1.01) {
        x += 0.1;
        tmp = pchisq(x, 1);
      }
    }
  }

  // Return p-value
  return x;
}

// #################
//  functions to calculate pchisqsum

double StatFunc::pchisqsum(double x, VectorXd lambda) {
  double pval = psadd(x, lambda);
  if (pval > 1.0) pval = psatt(x, lambda);
  return pval;
}

double StatFunc::psadd(double x, VectorXd lambda) {
  double d = lambda.maxCoeff();
  if (d <= 0.0) return 2.0;
  lambda = lambda.array() / d;
  x = x / d;

  double lmin = 0.0;
  double m = lambda.minCoeff();
  if (m < 0.0) lmin = 0.499995 / m;
  else if (x > lambda.sum()) lmin = -0.01;
  else lmin = -0.5 * (double)lambda.size() / x;
  double lmax = 0.499995 / lambda.maxCoeff();

  double hatzeta = Brents_Kp_min_x(lambda, x, lmin, lmax, 1e-08);
  if (hatzeta > lmax + 9) return 2.0;
  double sign = (hatzeta < 0.0) ? -1.0 : 1.0;
  double w = sign * sqrt(2 * (hatzeta * x - K(hatzeta, lambda)));
  double v = hatzeta * sqrt(Kpp(hatzeta, lambda));

  // debug
  // std::cout<<"hatzeta = "<<hatzeta<<std::endl;
  // std::cout<<"w = "<<w<<std::endl;
  // std::cout<<"v = "<<v<<std::endl;

  if (fabs(hatzeta) < 1e-04) return 2.0;
  else return pnorm(w + log(v / w) / w);
}

double StatFunc::psatt(double x, VectorXd lambda) {
  double sum = lambda.sum();
  if (CommFunc::FloatEqual(sum, 0.0)) return 2.0;

  double sq_sum = lambda.dot(lambda);
  double sum_sq = sum * sum;
  double a = sq_sum / sum;
  double b = sum_sq / sq_sum;

  return pchisq(x / a, b);
}

double StatFunc::K(double zeta, VectorXd& lambda) { return -0.5 * (1.0 - 2.0 * zeta * lambda.array()).log().sum(); }

double StatFunc::Kp(double zeta, VectorXd& lambda) {
  return (lambda.array() / (1.0 - 2.0 * zeta * lambda.array())).sum();
}

double StatFunc::Kpp(double zeta, VectorXd& lambda) {
  return 2.0 * (lambda.array().square() / (1.0 - 2.0 * zeta * lambda.array()).array().square()).sum();
}

double StatFunc::Kp_min_x(double zeta, VectorXd& lambda, double x) { return Kp(zeta, lambda) - x; }

double StatFunc::Brents_Kp_min_x(VectorXd& lambda, double x, double lowerLimit, double upperLimit, double errorTol) {
  double a = lowerLimit;
  double b = upperLimit;
  double c = 0;
  double d = 1.7976931348623157E+308;

  double fa = Kp_min_x(a, lambda, x);
  double fb = Kp_min_x(b, lambda, x);

  double fc = 0;
  double s = 0;
  double fs = 0;

  // if f(a) f(b) >= 0 then error-exit
  if (fa * fb >= 0) {
    if (fa < fb) return a;
    else return b;
  }

  // if |f(a)| < |f(b)| then std::swap (a,b) end if
  if (fabs(fa) < fabs(fb)) {
    double tmp = a;
    a = b;
    b = tmp;
    tmp = fa;
    fa = fb;
    fb = tmp;
  }

  c = a;
  fc = fa;
  bool mflag = true;
  int i = 0;

  while (!(fb == 0) && (fabs(a - b) > errorTol)) {
    if ((fa != fc) && (fb != fc))
      // Inverse quadratic interpolation
      s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) +
          c * fa * fb / (fc - fa) / (fc - fb);
    else
      // Secant Rule
      s = b - fb * (b - a) / (fb - fa);

    double tmp2 = (3 * a + b) / 4;
    if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) ||
        (!mflag && (fabs(s - b) >= (fabs(c - d) / 2)))) {
      s = (a + b) / 2;
      mflag = true;
    } else {
      if ((mflag && (fabs(b - c) < errorTol)) || (!mflag && (fabs(c - d) < errorTol))) {
        s = (a + b) / 2;
        mflag = true;
      } else mflag = false;
    }
    fs = Kp_min_x(s, lambda, x);
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    // if |f(a)| < |f(b)| then std::swap (a,b) end if
    if (fabs(fa) < fabs(fb)) {
      double tmp = a;
      a = b;
      b = tmp;
      tmp = fa;
      fa = fb;
      fb = tmp;
    }
    i++;
    if (i > 1000) return upperLimit + 10;
  }
  return b;
}
