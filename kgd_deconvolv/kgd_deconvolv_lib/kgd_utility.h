/*
 * kgd_deconvolv is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgd_deconvolv.
 *
 * kgd_deconvolv is free software: you can redistribute it and/or modify
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

#include "kgd_global.h"

#include "kgd_variant_index.h"
#include "kgd_exceptions.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace

/// A static singleton class to control visibility and access to misc. utility functions.
class Utility {

public:

  Utility() = delete;
  ~Utility() = delete;


  template<typename T>
  static std::vector<T> vecDiff(std::vector<T> &vecA, std::vector<T> &vecB) {

    assert(vecA.size() == vecB.size());

    std::vector<T> difference(vecA.size(), (T) 0);

    for (size_t i = 0; i < vecA.size(); i++) {

      difference[i] = vecA[i] - vecB[i];

    }

    return difference;

  }


  template<typename T>
  static std::vector<T> vecSum(std::vector<T> &vecA, std::vector<T> &vecB) {

    assert(vecA.size() == vecB.size());

    std::vector<T> tmpSum(vecA.size(), (T) 0);

    for (size_t i = 0; i < vecA.size(); i++) {

      tmpSum[i] = vecA[i] + vecB[i];

    }

    return tmpSum;
  }


  template<typename T>
  static std::vector<T> vecProd(const std::vector<T> &vecA, const std::vector<T> &vecB) {

    assert(vecA.size() == vecB.size());

    std::vector<T> tmpProd(vecA.size(), (T) 0);

    for (size_t i = 0; i < vecA.size(); i++) {

      tmpProd[i] = vecA[i] * vecB[i];

    }

    return tmpProd;

  }


  template<typename T>
  static T sumOfVec(const std::vector<T> &array) {

    T tmp = 0;

    for (auto const &value: array) {

      tmp += value;

    }

    return tmp;

  }

  static std::vector<std::vector<size_t> > uniqueMatrixColumns(std::vector<std::vector<size_t> > &matrix);

  static bool vectorEquivalence(std::vector<size_t> vec1, std::vector<size_t> vec2);

  static double max_value(std::vector<double> x);

  static double sumOfMat(std::vector<std::vector<double> > &matrix);

  static void normalizeBySum(std::vector<double> &array);

  static void normalizeBySumMat(std::vector<std::vector<double> > &matrix);

  static std::vector<double> calcLLKs(const std::vector<double> &refCount,
                                      const std::vector<double> &altCount,
                                      const std::vector<double> &expectedWsaf,
                                      size_t firstIndex,
                                      size_t length,
                                      double fac,
                                      double err = 0.01);

  static double calcLLK(double ref, double alt, double unadjustedWsaf, double err, double fac);

  static size_t sampleIndexGivenProp(const std::vector<double> &proportion);

  static std::vector<double> reshapeMatToVec(std::vector<std::vector<double> > &Mat);
// Returns trimmed string.
  static std::string trimWhiteSpace(const std::string& s);


private:


};



}   // organization level namespace
}   // project level namespace


#endif
