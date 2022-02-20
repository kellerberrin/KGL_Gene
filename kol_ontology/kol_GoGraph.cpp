//
// Created by kellerberrin on 20/2/22.
//

#include "kol_GoGraph.h"
#include "kol_GoGraphImpl.h"

namespace kol = kellerberrin::ontology;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implements simple PIMPL encapsulation for the boost acyclic graph object.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kol::GoGraph::GoGraph(std::unique_ptr<GoGraphImpl>&& go_graph_impl_ptr) : go_graph_impl_ptr_(std::move(go_graph_impl_ptr)) {}

kol::GoGraph::~GoGraph() {

  go_graph_impl_ptr_ = nullptr;

}


const kol::GoGraphImpl& kol::GoGraph::getGoGraphImpl() const {

  return *go_graph_impl_ptr_;


}

