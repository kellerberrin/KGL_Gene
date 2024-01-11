//
// Created by kellerberrin on 25/12/23.
//

#ifndef KGL_SEQUENCE_VIEW_TREE_H
#define KGL_SEQUENCE_VIEW_TREE_H


#include "kgl_sequence_amino.h"
#include "kgl_sequence_base.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_distance_tree_node.h"

#include <concepts>
#include <cmath>


namespace kellerberrin::genome {   //  organization level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Important - for performance reasons; only sequence views are stored.
// The underlying sequences must remain on memory while this tree node exists or a seg-fault will surely occur.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T>
concept NodeSequenceView = ( std::same_as<T, AminoSequenceView>
                          || std::same_as<T, DNA5SequenceLinearView>
                          || std::same_as<T, DNA5SequenceCodingView>);


template<NodeSequenceView SequenceViewType>
class SequenceViewTreeNode : public TreeNodeDistance {

public:

  SequenceViewTreeNode(const SequenceViewType& sequence_view, std::string sequence_tag, SequenceDistanceMetric<SequenceViewType> sequence_distance)
  : sequence_view_(sequence_view), sequence_tag_(std::move(sequence_tag)), sequence_distance_(sequence_distance) {}
  ~SequenceViewTreeNode() override = default;

  [[nodiscard]] std::string nodeText() const override { return sequence_tag_; }
  [[nodiscard]] DistanceType_t distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const override;

private:

  SequenceViewType sequence_view_;
  std::string sequence_tag_;
  SequenceDistanceMetric<SequenceViewType> sequence_distance_;

};


template<NodeSequenceView SequenceViewType>
DistanceType_t SequenceViewTreeNode<SequenceViewType>::distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const {

  auto sequence_node_ptr = std::dynamic_pointer_cast<const SequenceViewTreeNode>(distance_node);
  if (not sequence_node_ptr) {

    ExecEnv::log().error("Mismatched sequence types");
    return 0.0;

  }

  DistanceType_t distance = sequence_distance_(sequence_view_, sequence_node_ptr->sequence_view_);
  if (std::isnan(distance)) {

    ExecEnv::log().warn("Invalid distance (NaN) between sequence node: {} and node: {}", nodeText(), sequence_node_ptr->nodeText());
    distance = 0.0;

  }

  return distance;

}


using AminoSequenceViewNode = SequenceViewTreeNode<AminoSequenceView>;
using CodingSequenceViewNode = SequenceViewTreeNode<DNA5SequenceCodingView>;
using LinearSequenceViewNode = SequenceViewTreeNode<DNA5SequenceLinearView>;



} // namespace.


#endif //KGL_SEQUENCE_VIEW_TREE_H
