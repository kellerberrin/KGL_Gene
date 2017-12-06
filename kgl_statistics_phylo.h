//
// Created by kellerberrin on 6/12/17.
//

#ifndef KGL_STATISTICS_PHYLO_H
#define KGL_STATISTICS_PHYLO_H

#include "kgl_statistics.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class T> class PhyloNode;  // fwd.
using DistanceType_t = uint64_t;
template <class T> using OutNodes = std::multimap<DistanceType_t , std::shared_ptr<PhyloNode<T>>>;

template<class T> class PhyloNode {

public:

  explicit PhyloNode(std::shared_ptr<const T> leaf, DistanceType_t distance) : leaf_(leaf), distance_(distance) {}
  ~PhyloNode() = default;

  DistanceType_t distance() const { return distance_; }
  void addNode(std::shared_ptr<PhyloNode<T>> node) { out_nodes_.insert(std::pair<DistanceType_t , std::shared_ptr<PhyloNode<T>>>(node->distance(), node)); }
  std::shared_ptr<T> leaf() const { return leaf_; }
  const OutNodes<T>& getMap() const { return out_nodes_; }

private:

  DistanceType_t distance_;
  OutNodes<T> out_nodes_;
  std::shared_ptr<T> leaf_;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class T> using NodeVector = std::vector<std::shared_ptr<PhyloNode<T>>>;

template<class T> class UPGMAMatrix {

public:

  explicit UPGMAMatrix(const NodeVector<T>& node_vector) : node_vector_(node_vector) {}
  ~UPGMAMatrix() = default;

  void calculateReduce();
  DistanceType_t distance(std::shared_ptr<PhyloNode<T>> row_node, std::shared_ptr<PhyloNode<T>> column_node) const;

private:

  NodeVector<T> node_vector_;

};


template<class T> void UPGMAMatrix<T>::calculateReduce() {

  while (node_vector_.size() > 1) {

    struct {
      size_t row_index;
      size_t column_index;
      DistanceType_t distance;
    } minimum_distance;

    bool first_pass = true;
    for (size_t row = 0; row < node_vector_.size(); ++row) {

      for (size_t column = 0; column < row; ++column) {

        if (first_pass) {

          minimum_distance.row_index = row;
          minimum_distance.column_index = column;
          minimum_distance.distance = distance(node_vector_[row], node_vector_[column]);
          first_pass = false;

        } else {

          DistanceType_t node_distance = distance(node_vector_[row], node_vector_[column]);
          if (node_distance < minimum_distance.distance) {

            minimum_distance.row_index = row;
            minimum_distance.column_index = column;
            minimum_distance.distance = node_distance;

          } // if new minimum

        } // else first pass

      } // for column

    } // for row

    // we have calculated the minimum distance pair, so now reduce them.



  } // while reduce.

}



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_STATISTICS_PHYLO_H
