//
// kgl_distance_matrix.h — symmetric distance matrix, PIMPL pattern.
//
// Created by kellerberrin on 30/10/23.
//

#ifndef KGL_DISTANCE_MATRIX_H
#define KGL_DISTANCE_MATRIX_H

#include <memory>
#include <tuple>

namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance matrix. Implements the PIMPL pattern; the packed triangular
// implementation is defined in kgl_distance_matrix.cpp.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


using DistanceType_t = double;
class DistanceMatrixImpl;       // Forward declaration of the distance matrix implementation class.

class DistanceMatrix {

public:

  DistanceMatrix();
  explicit DistanceMatrix(size_t matrix_size);
  DistanceMatrix(DistanceMatrix&& matrix) noexcept;
  DistanceMatrix& operator=(DistanceMatrix&& matrix) noexcept;
  ~DistanceMatrix();  // Do not use the default destructor, see PIMPL fwd decl.

  void addMatrix(const DistanceMatrix& add_matrix);
  void copyMatrix(const DistanceMatrix& copy_matrix);

  // Tuple returns value, row index, column index in that order.
  [[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> minimum() const;
  [[nodiscard]] std::tuple<DistanceType_t, size_t, size_t> maximum() const;

  [[nodiscard]] DistanceType_t getDistance(size_t row, size_t column) const;
  void setDistance(size_t row, size_t column, DistanceType_t distance);

  [[nodiscard]] size_t size() const;
  void resize(size_t new_size);

  // Rescale elements to the interval [0, 1].
  void normalizeDistance();

private:

  std::unique_ptr<DistanceMatrixImpl> impl_ptr_;    // PIMPL

  [[nodiscard]] const DistanceMatrixImpl& impl() const noexcept { return *impl_ptr_; }
  [[nodiscard]] DistanceMatrixImpl& impl() noexcept { return *impl_ptr_; }

};


} // Namespace.


#endif //KGL_DISTANCE_MATRIX_H