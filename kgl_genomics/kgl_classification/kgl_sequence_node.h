//
// Created by kellerberrin on 11/01/24.
//

#ifndef KGL_SEQUENCE_NODE_H
#define KGL_SEQUENCE_NODE_H



#include "kgl_sequence_amino.h"
#include "kgl_sequence_base.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_distance_tree_node.h"

#include <concepts>
#include <cmath>


namespace kellerberrin::genome {   //  organization level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Important - for actual sequences (not sequence views)
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T>
concept NodeSequence = ( std::same_as<T, AminoSequence> || std::same_as<T, DNA5SequenceLinear> || std::same_as<T, DNA5SequenceCoding>);


template<NodeSequence SequenceType>
class SequenceTreeNode : public TreeNodeDistance {

public:

  SequenceTreeNode(SequenceType&& sequence, std::string sequence_tag, SequenceDistanceMetric<SequenceType> sequence_distance)
      : sequence_(std::move(sequence)), sequence_tag_(std::move(sequence_tag)), sequence_distance_(sequence_distance) {}
  ~SequenceTreeNode() override = default;

  [[nodiscard]] std::string nodeText() const override { return sequence_tag_; }
  [[nodiscard]] DistanceType_t distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const override;

private:

  SequenceType sequence_;
  std::string sequence_tag_;
  SequenceDistanceMetric<SequenceType> sequence_distance_;

};


template<NodeSequence SequenceType>
DistanceType_t SequenceTreeNode<SequenceType>::distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const {

  auto sequence_node_ptr = std::dynamic_pointer_cast<const SequenceTreeNode>(distance_node);
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


using AminoSequenceNode = SequenceTreeNode<AminoSequence>;
using CodingSequenceNode = SequenceTreeNode<DNA5SequenceCoding>;
using LinearSequenceNode = SequenceTreeNode<DNA5SequenceLinear>;


} // namespace.






#endif //KGL_SEQUENCE_NODE_H
