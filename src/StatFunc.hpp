/*
 * Interface to the statistical functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _STATFUNC_H
#define _STATFUNC_H

#include <vector>

#include <Eigen/Dense>

using namespace Eigen;

namespace StatFunc {
////////// P-value Calculatiion Functions Start ////////////////
double gammln(double xx);

double chi_prob(double df, double chi_sqr_val);

double gammp(const double a, const double x);

void gser(double& gamser, const double a, const double x, double& gln);

void gcf(double& gammcf, const double a, const double x, double& gln);
////////// P-value Calculatiion Functions End ////////////////

// Function for get a Chi value of right tail when given a df and prob.
double chi_val(double df, double prob);

// normal distribution
double pnorm(double x);

// chisq distribution
double pchisq(double x, double df);
double qchisq(double q, double df);

// sum of chisq distribution
double pchisqsum(double x, VectorXd lambda);
double psadd(double x, VectorXd lambda);
double psatt(double x, VectorXd lambda);
double K(double zeta, VectorXd& lambda);
double Kp(double zeta, VectorXd& lambda);
double Kpp(double zeta, VectorXd& lambda);
double Kp_min_x(double zeta, VectorXd& lambda, double x);
double Brents_Kp_min_x(VectorXd& lambda, double x, double lowerLimit, double upperLimit, double errorTol);
}  // namespace StatFunc

#endif
