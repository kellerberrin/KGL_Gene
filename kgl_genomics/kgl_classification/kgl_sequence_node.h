//
// Created by kellerberrin on 25/12/23.
//

#ifndef KGL_SEQUENCE_TREE_H
#define KGL_SEQUENCE_TREE_H


#include "kgl_sequence_amino.h"
#include "kgl_sequence_base.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_distance_tree_node.h"

#include <concepts>
#include <cmath>

namespace kellerberrin::genome {   //  organization level namespace


template<typename T>
concept NodeSequence = (std::same_as<T, AminoSequence> || std::same_as<T, DNA5SequenceLinear> || std::same_as<T, DNA5SequenceCoding>);


template<NodeSequence SequenceType>
class SequenceNode : public TreeNodeDistance {

public:

  SequenceNode(SequenceType&& sequence, std::string sequence_tag, SequenceDistanceMetric<SequenceType> sequence_distance)
  : sequence_(std::move(sequence)), sequence_tag_(std::move(sequence_tag)), sequence_distance_(sequence_distance) {}
  ~SequenceNode() override = default;

  [[nodiscard]] std::string nodeText() const override { return sequence_tag_; }
  [[nodiscard]] DistanceType_t distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const override;

private:

  SequenceType sequence_;
  std::string sequence_tag_;
  SequenceDistanceMetric<SequenceType> sequence_distance_;

};


template<NodeSequence SequenceType>
DistanceType_t SequenceNode<SequenceType>::distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const {

  auto sequence_node_ptr = std::dynamic_pointer_cast<const SequenceNode>(distance_node);
  if (not sequence_node_ptr) {

    ExecEnv::log().error("Mismatched sequence types");
    return 0.0;

  }

  DistanceType_t distance = sequence_distance_(sequence_, sequence_node_ptr->sequence_);
  if (std::isnan(distance)) {

    ExecEnv::log().warn("Invalid distance (NaN) between sequence node: {} and node: {}", nodeText(), sequence_node_ptr->nodeText());
    distance = 0.0;

  }

  return distance;

}


using AminoSequenceNode = SequenceNode<AminoSequence>;
using CodingSequenceNode = SequenceNode<DNA5SequenceCoding>;
using LinearSequenceNode = SequenceNode<DNA5SequenceLinear>;



} // namespace.


#endif //KGL_SEQUENCE_TREE_H
