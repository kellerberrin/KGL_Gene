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

#include "logbeta.h"
#include <iterator>     // std::distance
#include <algorithm> // find
#include "kgd_utility.h"
#include "kgd_deconvolv_app.h"
#include "kgd_distribution.h"
#include "loggammasum.h" // which includes log_gamma.h
#include "gamma.h"

#include <boost/math/special_functions/gamma.hpp>

namespace bm = boost::math;
namespace kgd = kellerberrin::deconvolv;



std::vector<std::vector<size_t> > kgd::Utility::uniqueMatrixColumns(std::vector<std::vector<size_t> > &matrix) {

  std::vector<std::vector<size_t> > ret;

  ret.push_back(matrix[0]);

  for (size_t i = 1; i < matrix.size(); i++) {

    bool is_new_column = true;

    for (auto state : ret) {

      if (vectorEquivalence(state, matrix[i])) {

        is_new_column = false;
        break;

      }

    }

    if (is_new_column) {

      ret.push_back(matrix[i]);

    }

  }

  return ret;

}


bool kgd::Utility::vectorEquivalence(std::vector<size_t> vec1, std::vector<size_t> vec2) {

  if (vec1.size() != vec2.size()) {

    return false;

  }

  bool ret = true;

  for (size_t i = 0; i < vec1.size(); i++) {

    if (vec1[i] != vec2[i]) {

      ret = false;
      break;

    }

  }

  return ret;

}


double kgd::Utility::max_value(std::vector<double> x) {

  assert(x.size() > 0);

  auto tmpMaxIt = std::max_element(std::begin(x), std::end(x));

  return *tmpMaxIt;

}


double kgd::Utility::sumOfMat(std::vector<std::vector<double> > &matrix) {

  double tmp = 0.0;

  for (auto const &array: matrix) {

    for (auto const &value: array) {

      tmp += value;

    }

  }

  return tmp;
}


void kgd::Utility::normalizeBySum(std::vector<double> &array) {

  double inverse_sum = 1.0 / sumOfVec(array);

  for (size_t idx = 0; idx < array.size(); ++idx) {

    array[idx] = array[idx] * inverse_sum;

  }

}


void kgd::Utility::normalizeBySumMat(std::vector<std::vector<double> > &matrix) {

  double tmpsum = sumOfMat(matrix);

  for (size_t i = 0; i < matrix.size(); i++) {

    for (std::vector<double>::iterator it = matrix[i].begin(); it != matrix[i].end(); ++it) {

      *it /= tmpsum;

    }

  }

}


std::vector<double> kgd::Utility::calcLLKs(const std::vector<double> &refCount,
                                           const std::vector<double> &altCount,
                                           const std::vector<double> &expectedWsaf,
                                           size_t firstIndex,
                                           size_t length,
                                           double fac,
                                           double err) {

  assert (expectedWsaf.size() == length);

  std::vector<double> tmpLLKs(length, 0.0);

  size_t index = firstIndex;

  for (size_t i = 0; i < length; i++) {

    assert (expectedWsaf[i] >= 0);
    //assert (expectedWsaf[i] <= 1);

    tmpLLKs[i] = calcLLK(refCount[index],
                         altCount[index],
                         expectedWsaf[i],
                         err,
                         fac);

    index++;

  }

  return tmpLLKs;

}


double kgd::Utility::calcLLK(double ref, double alt, double unadjustedWsaf, double err, double fac) {

  double adjustedWsaf = unadjustedWsaf + ((1 - (2 * unadjustedWsaf)) * err);

  double a2_arg = adjustedWsaf * fac;
  double b2_arg = (1 - adjustedWsaf) * fac;

  if (b2_arg <= 0 or a2_arg <= 0) {

    ExecEnv::log().error("calcLLK(); BetaBinomialDistribution Invalid Alpha:{} or Beta: {}, unadjustedWsaf: {}, fac: {}, adjustedWsaf: {}",
                         a2_arg, b2_arg, unadjustedWsaf, fac, adjustedWsaf);

  }

#define BOOST_BBD_  1
#ifdef BOOST_BBD_

  double prob = BetaBinomialDistribution::logPartialPdf(ref+alt, alt, a2_arg, b2_arg);

  return prob;

#else

  double a1_arg = alt + a2_arg;
  double b1_arg = ref + b2_arg;

  //Binomial coefficients are expressed in terms of the (log) beta function
  double inverse_log_binonimal_coeff = std::log(alt + ref + 1.0) + Maths::Special::Gamma::logBeta(ref+1.0, alt+1.0);
  double llk = Maths::Special::Gamma::logBeta(a1_arg, b1_arg) - Maths::Special::Gamma::logBeta(a2_arg, b2_arg) - inverse_log_binonimal_coeff;

  return llk;

#endif

}




size_t kgd::Utility::sampleIndexGivenProp(const std::vector<double>& proportion) {

// Linear search is costly.

  auto biggest = std::max_element(std::begin(proportion), std::end(proportion));

  return std::distance(proportion.begin(), biggest);

#ifdef LEGACY_CODE_

  std::vector <size_t> tmpIndex;

  for ( size_t i = 0; i < proportion.size(); i++ ){

      tmpIndex.push_back(i);

  }

  std::vector <double> tmpCdf = computeCdf(proportion);

  double u = randomGenerator->sample();



  for (size_t i = 0; ; i < tmpCdf.size() ; i++){

      if ( u < tmpCdf[i] ){
          break;

      }

  }

  return i;

#endif

}


std::vector<double> kgd::Utility::reshapeMatToVec(std::vector<std::vector<double> > &Mat) {

  std::vector<double> tmp;

  for (auto const &array: Mat) {

    for (auto const &value: array) {

      tmp.push_back(value);

    }

  }

  return tmp;

}



// Returns trimmed string.
std::string kgd::Utility::trimWhiteSpace(const std::string& s) {

  std::string clean_string;
  auto lambda_not_whitespace = [](unsigned char c){ return not std::isspace(c); };
  std::copy_if(s.begin(), s.end(), back_inserter(clean_string), lambda_not_whitespace);

  return clean_string;

}
