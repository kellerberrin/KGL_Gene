//
// Created by kellerberrin on 20/2/22.
//

#ifndef KOL_GO_GRAPH
#define KOL_GO_GRAPH

#include <memory>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// PIMPL class to hide and localize the boost acyclic graph implementation details.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::ontology {


// Forward for the implementation class.
class GoGraphImpl;


class GoGraph {

public:


  explicit GoGraph(std::unique_ptr<GoGraphImpl>&& go_graph_impl_ptr);
  ~GoGraph();

  [[nodiscard]] const GoGraphImpl& getGoGraphImpl() const;


private:

  std::unique_ptr<GoGraphImpl> go_graph_impl_ptr_;

};



} // namespace


#endif //KGL_KOL_GOGRAPH_H
