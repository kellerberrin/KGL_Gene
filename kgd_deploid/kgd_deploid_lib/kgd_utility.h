/*
 * kgd_deploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgd_deploid.
 *
 * kgd_deploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef KGD_UTILITY_H
#define KGD_UTILITY_H


#include <vector>
#include <iostream>
#include <algorithm>    /* min_element, max_element */

#include "kgd_mersenne_twister.h"
#include "kgd_global.h"

#include "kgd_variantIndex.h"
#include "kgd_exceptions.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace


template<typename T>
std::vector<T> vecDiff(std::vector<T> &vecA, std::vector<T> &vecB) {

  assert(vecA.size() == vecB.size());

  std::vector<T> difference(vecA.size(), (T) 0);

  for (size_t i = 0; i < vecA.size(); i++) {

    difference[i] = vecA[i] - vecB[i];

  }

  return difference;

}


template<typename T>
std::vector<T> vecSum(std::vector<T> &vecA, std::vector<T> &vecB) {

  assert(vecA.size() == vecB.size());

  std::vector<T> tmpSum(vecA.size(), (T) 0);

  for (size_t i = 0; i < vecA.size(); i++) {

    tmpSum[i] = vecA[i] + vecB[i];

  }

  return tmpSum;
}


template<typename T>
std::vector<T> vecProd(std::vector<T> &vecA, std::vector<T> &vecB) {

  assert(vecA.size() == vecB.size());

  std::vector<T> tmpProd(vecA.size(), (T) 0);

  for (size_t i = 0; i < vecA.size(); i++) {

    tmpProd[i] = vecA[i] * vecB[i];

  }

  return tmpProd;

}


template<typename T>
T sumOfVec(std::vector<T> &array) {

  T tmp = 0;

  for (auto const &value: array) {

    tmp += value;

  }

  return tmp;

}

/*! \brief Compute factorial of a \return double a! */
template<class T>
T factorial(T a) {

  if (a > 1) return (a * factorial(a - 1));

  else return (1);

}

/*! \brief Compute a permutations of n \return double */
template<class T>
T n_permu_a(T n, T a) {

  if (a > 1) return (n * n_permu_a(n - 1, a - 1));

  else if (a == 1) return (n);

  else return (1);

}

/*! \brief Compute n choose k \return double */
template<class T>
T n_choose_k(T n, T k) {

  if (k < (n / 2)) return (n_choose_k(n, n - k));

  else return (n_permu_a(n, k) / factorial(k));

}


double normal_pdf(double x, double m, double s);

double min_value(std::vector<double> x);

double max_value(std::vector<double> x);

std::vector<double> computeCdf(std::vector<double> &dist);

double sumOfMat(std::vector<std::vector<double> > &matrix);

void normalizeBySum(std::vector<double> &array);

void normalizeByMax(std::vector<double> &array);

void normalizeBySumMat(std::vector<std::vector<double> > &matrix);

std::vector<double>
calcLLKs(std::vector<double> &refCount, std::vector<double> &altCount, std::vector<double> &expectedWsaf,
         size_t firstIndex, size_t length, double fac, double err = 0.01);

double calcLLK(double ref, double alt, double unadjustedWsaf, double err, double fac);

size_t sampleIndexGivenProp(RandomGenerator *rg, std::vector<double> proportion);

std::vector<double> reshapeMatToVec(std::vector<std::vector<double> > &Mat);

double betaPdf(double x, double a, double b);

double logBetaPdf(double x, double a, double b);

double binomialPdf(int s, int n, double p);

//double betaDistConst( double a , double b);
double rBeta(double alpha, double beta, RandomGenerator *rg);


}   // organization level namespace
}   // project level namespace


#endif
