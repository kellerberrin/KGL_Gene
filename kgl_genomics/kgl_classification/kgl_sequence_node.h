//
// kgl_sequence_node.h — sequence tree nodes for UPGMA classification.
//
// Created by kellerberrin on 11/01/24.
//

#ifndef KGL_SEQUENCE_NODE_H
#define KGL_SEQUENCE_NODE_H

#include "kgl_sequence.h"
#include "kgl_genetic_code.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_distance_tree_node.h"

#include <cmath>
#include <concepts>


namespace kellerberrin::genome {   //  organization level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Sequence classification nodes.
//
// Important - for performance reasons, the view variants only store sequence views.
// The underlying sequences must remain in memory while the view tree nodes exist
// or a seg-fault will surely occur.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T>
concept NodeSequence = ( std::same_as<T, AminoSequence> || std::same_as<T, DNA5SequenceLinear> || std::same_as<T, DNA5SequenceCoding>);

template<typename T>
concept NodeSequenceView = ( std::same_as<T, AminoSequenceView>
                          || std::same_as<T, DNA5SequenceLinearView>
                          || std::same_as<T, DNA5SequenceCodingView>);

// Unifies the reference's SequenceTreeNode and SequenceViewTreeNode into a single template.
template<NodeSequence SequenceType>
class SequenceClassificationNode : public TreeNodeDistance {

public:

  SequenceClassificationNode(SequenceType&& sequence, std::string sequence_tag, SequenceDistanceMetric<SequenceType> sequence_distance)
      : sequence_(std::move(sequence)), sequence_tag_(std::move(sequence_tag)), sequence_distance_(sequence_distance) {}
  ~SequenceClassificationNode() override = default;

  [[nodiscard]] std::string nodeText() const override { return sequence_tag_; }
  [[nodiscard]] DistanceType_t distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const override {

    auto sequence_node_ptr = std::dynamic_pointer_cast<const SequenceClassificationNode>(distance_node);
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

private:

  SequenceType sequence_;
  std::string sequence_tag_;
  SequenceDistanceMetric<SequenceType> sequence_distance_;

};


// Unifies the reference's SequenceViewTreeNode in the same template shape.
template<NodeSequenceView SequenceViewType>
class SequenceViewClassificationNode : public TreeNodeDistance {

public:

  SequenceViewClassificationNode(const SequenceViewType& sequence_view, std::string sequence_tag, SequenceDistanceMetric<SequenceViewType> sequence_distance)
  : sequence_view_(sequence_view), sequence_tag_(std::move(sequence_tag)), sequence_distance_(sequence_distance) {}
  ~SequenceViewClassificationNode() override = default;

  [[nodiscard]] std::string nodeText() const override { return sequence_tag_; }
  [[nodiscard]] DistanceType_t distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const override {

    auto sequence_node_ptr = std::dynamic_pointer_cast<const SequenceViewClassificationNode>(distance_node);
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

private:

  SequenceViewType sequence_view_;
  std::string sequence_tag_;
  SequenceDistanceMetric<SequenceViewType> sequence_distance_;

};


using AminoSequenceNode = SequenceClassificationNode<AminoSequence>;
using CodingSequenceNode = SequenceClassificationNode<DNA5SequenceCoding>;
using LinearSequenceNode = SequenceClassificationNode<DNA5SequenceLinear>;

using AminoSequenceViewNode = SequenceViewClassificationNode<AminoSequenceView>;
using CodingSequenceViewNode = SequenceViewClassificationNode<DNA5SequenceCodingView>;
using LinearSequenceViewNode = SequenceViewClassificationNode<DNA5SequenceLinearView>;


} // namespace.


#endif //KGL_SEQUENCE_NODE_H