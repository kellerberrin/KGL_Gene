//
// Created by kellerberrin on 16/12/17.
//


#include <boost/numeric/ublas/triangular.hpp>
// #include <boost/numeric/ublas/io.hpp>

#include "kgl_exec_env.h"
#include "kgl_statistics_upgma.h"


namespace kgl = kellerberrin::genome;
namespace bnu = boost::numeric::ublas;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of a strict lower triangular distance matrix using the Boost library.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using DistanceMatrixImplType = bnu::triangular_matrix<kgl::DistanceType_t, bnu::strict_lower>;

class kgl::DistanceMatrix::BoostDistanceMatrix {

public:

  explicit BoostDistanceMatrix(size_t matrix_size) : lower_triangular_(matrix_size, matrix_size) {}
  explicit BoostDistanceMatrix(const BoostDistanceMatrix&) = default;
  ~BoostDistanceMatrix() = default;

  size_t size() const { return lower_triangular_.size1(); }
  void resize(size_t new_size) { lower_triangular_.resize(new_size, new_size, false); }
  kgl::DistanceType_t getDistance(size_t i, size_t j) const;
  void setDistance(size_t i, size_t j, kgl::DistanceType_t distance);

private:

  DistanceMatrixImplType lower_triangular_;

};


kgl::DistanceType_t kgl::DistanceMatrix::BoostDistanceMatrix::getDistance(size_t i, size_t j) const {

  if (i >= size() or i >= size()) {

    kgl::ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
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

void kgl::DistanceMatrix::BoostDistanceMatrix::setDistance(size_t i, size_t j, kgl::DistanceType_t distance) {

  if (i >= size() or i >= size()) {

    kgl::ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public class of the UPGMA distance matrix.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::DistanceMatrix::DistanceMatrix(size_t matrix_size) : diagonal_impl_ptr_(std::make_unique<kgl::DistanceMatrix::BoostDistanceMatrix>(matrix_size)) {}

kgl::DistanceMatrix::DistanceMatrix(const DistanceMatrix& copy) : diagonal_impl_ptr_(std::make_unique<kgl::DistanceMatrix::BoostDistanceMatrix>(*copy.diagonal_impl_ptr_)) {}

kgl::DistanceMatrix::~DistanceMatrix() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete PIMPL type.


size_t kgl::DistanceMatrix::size() const {

  return diagonal_impl_ptr_->size();

}


void kgl::DistanceMatrix::resize(size_t new_size) {

  diagonal_impl_ptr_->resize(new_size);

}


kgl::DistanceType_t kgl::DistanceMatrix::getDistance(size_t i, size_t j) const {

  return diagonal_impl_ptr_->getDistance( i, j);

}

void kgl::DistanceMatrix::setDistance(size_t i, size_t j, DistanceType_t value) {

  diagonal_impl_ptr_->setDistance( i, j, value);

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DistanceType_t kgl::DistanceMatrix::minimum(size_t& i, size_t& j) const {

  bool first_pass = true;
  kgl::DistanceType_t min_distance = 0;
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

        kgl::DistanceType_t check_min = getDistance( row, column);
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
void kgl::DistanceMatrix::reduce(size_t i, size_t j) {

  // Save and resize
  DistanceMatrix temp_distance(*this);

  size_t reduce_size = size() - 1;

  resize(reduce_size);

  // re-populate merged distances.
  size_t idx_row = 1;
  for(size_t row = 0; row < temp_distance.size(); ++row) {

    if (row != i and row != j) {

      kgl::DistanceType_t i_leaf = getLeafCount(i);
      kgl::DistanceType_t j_leaf = getLeafCount(j);
      kgl::DistanceType_t calc_dist = (i_leaf * temp_distance.getDistance(row, i)) + (j_leaf * temp_distance.getDistance(row, j));
      calc_dist = calc_dist / (i_leaf + j_leaf);
      setDistance(idx_row, 0, calc_dist);
      ++idx_row;

    }

  }

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

}
