//
// Created by kellerberrin on 9/12/17.
//

#include <map>
#include <iostream>
#include <fstream>
#include <utility>

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "kgl_exec_env.h"
#include "kgl_statistics_phylo.h"


namespace kgl = kellerberrin::genome;
namespace bnu = boost::numeric::ublas;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of a strict lower triangular distance matrix using the Boost library.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using DistanceMatrixImplType = bnu::triangular_matrix<kgl::DistanceType_t, bnu::strict_lower>;

class kgl::DistanceMatrix::DistanceMatrixImpl {

public:

  explicit DistanceMatrixImpl(size_t matrix_size) : lower_triangular_(matrix_size, matrix_size) {}
  explicit DistanceMatrixImpl(DistanceMatrixImpl& distance_matrix) = default;
  ~DistanceMatrixImpl() = default;

  DistanceMatrixImpl& operator=(const DistanceMatrixImpl& distance_matrix) = default;

  size_t size() const { return lower_triangular_.size1(); }
  void reduce(size_t i, size_t j);
  DistanceType_t findMin(size_t& i, size_t& j) const;
  DistanceType_t getDistance(size_t i, size_t j) const;
  void setDistance(size_t i, size_t j, DistanceType_t distance);

private:

  DistanceMatrixImplType lower_triangular_;

};


kgl::DistanceType_t kgl::DistanceMatrix::DistanceMatrixImpl::getDistance(size_t i, size_t j) const {

  if (i >= size() or i >= size()) {

    ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
    return 0;

  }

  if (i == j) {

    return 0;

  }

  if (j > i) {

    std::swap(i, j);

  }

  return lower_triangular_(i,j);

}

void kgl::DistanceMatrix::DistanceMatrixImpl::setDistance(size_t i, size_t j, DistanceType_t distance) {

  if (i >= size() or i >= size()) {

    ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
    return;

  }

  if (i == j) {

    return;

  }

  if (j > i) {

    std::swap(i, j);

  }

  lower_triangular_(i,j) = distance;

}

kgl::DistanceType_t kgl::DistanceMatrix::DistanceMatrixImpl::findMin(size_t& i, size_t& j) const {

  bool first_pass = true;
  DistanceType_t min_distance = 0;
//  i = 0;
//  j = 0;
  for (size_t row = 0; row < size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      if (first_pass) {

        min_distance = getDistance( row, column);
        i = row;
        j = column;
        first_pass = false;

      } else {

        DistanceType_t check_min = getDistance( row, column);
        if (check_min < min_distance) {

          min_distance = check_min;
          i = row;
          j = column;

        }

      }

    }

  }

  return min_distance;

}

// Reduces the distance matrix.
// The reduced column is the left most column (column, j index = 0)
void kgl::DistanceMatrix::DistanceMatrixImpl::reduce(size_t i, size_t j) {

  std::cout << "merge i:" << i << " j:" << j << " initial matrix" << lower_triangular_ << "\n";

  // Save and resize
  DistanceMatrixImpl temp_distance(*this);

  size_t reduce_size = size() - 1;
  lower_triangular_.resize(reduce_size, reduce_size, false);

  // re-populate merged distances.
  size_t idx_row = 1;
  for(size_t row = 0; row < temp_distance.size(); ++row) {

    if (row != i and row != j) {

      DistanceType_t calc_dist = 0.5 * (temp_distance.getDistance(row, i) + temp_distance.getDistance(row, j));
      setDistance(idx_row, 0, calc_dist);
      ++idx_row;

    }

  }

  std::cout << "merged distances, reduced matrix:" << lower_triangular_ << "\n";

  // re-populate other distances.
  if (size() <= 2) {

    return;

  }

  idx_row = 2;  // shift down 2 rows.
  bool update = false;
  for (size_t row = 0; row < temp_distance.size(); ++row) {

    if (row != i and row != j) {

      size_t idx_column = 1;  // shift to column = 1
      for (size_t column = 0; column < row; ++column) {

        if (column != i and column != j) {

          setDistance(idx_row, idx_column, temp_distance.getDistance(row, column));
          idx_column++;
          update = true;

        }

      }

      if (update) {

        idx_row++;
        update = false;

      }

    }

  }

  std::cout << "final reduced matrix:" << lower_triangular_ << "\n";

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public class of a strict lower triangular distance matrix.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DistanceMatrix::DistanceMatrix(size_t matrix_size) : diagonal_impl_ptr_(std::make_unique<kgl::DistanceMatrix::DistanceMatrixImpl>(matrix_size)) {}
kgl::DistanceMatrix::~DistanceMatrix() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete PIMPL type.

size_t kgl::DistanceMatrix::size() const {

  return diagonal_impl_ptr_->size();

}


kgl::DistanceType_t kgl::DistanceMatrix::getDistance(size_t i, size_t j) const {

  return diagonal_impl_ptr_->getDistance( i, j);

}

void kgl::DistanceMatrix::setDistance(size_t i, size_t j, DistanceType_t value) {

  diagonal_impl_ptr_->setDistance( i, j, value);

}

kgl::DistanceType_t kgl::DistanceMatrix::reduceMinimum(size_t& i, size_t& j) const {

  DistanceType_t minimum;

  minimum = diagonal_impl_ptr_->findMin( i, j);
  ExecEnv::log().info("findMin() i: {}, j: {} size: {}, minimum: {}", i, j, size(), minimum);
  diagonal_impl_ptr_->reduce(i, j);

  return minimum;

}

