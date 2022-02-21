//
// Created by kellerberrin on 16/4/21.
//

#include "contrib/kol_GoGraphImpl.h"
#include "kel_exec_env.h"

namespace kol = kellerberrin::ontology;


//! A Breath first search visitor that creates the induced subgraph of a graph
/*!
  This class extends breadth first search and adds a vertex to a subgraph
    of the graph being visited.
*/

class SubgraphBFSVisitor : public boost::default_bfs_visitor {

public:
  SubgraphBFSVisitor(kol::GoGraphImpl::Graph &sub) : subgraph(sub) {}

  template<typename Vertex, typename Graph>
  void discover_vertex(Vertex u, const Graph &) {
    boost::add_vertex(u, subgraph);
  }

  kol::GoGraphImpl::Graph &subgraph;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kol::GoTermRecord::clearRecord() {

  term_id_.clear();
  name_.clear();
  definition_.clear();
  ontology_ = GO::Ontology::ONTO_ERROR;
  relations_.clear();
  alt_id_.clear();
  attributes_.clear();

}

// Check the record for integrity
bool kol::GoTermRecord::validRecord() const {

  bool valid;
  valid = not term_id_.empty();
  valid = valid and not name_.empty();
  valid = valid and not definition_.empty();
  valid = valid and ontology_ != kol::GO::Ontology::ONTO_ERROR;
  // All nodes except root nodes must have relations to another node.
  valid = valid and (not relations().empty() or GO::isRootTerm(term_id_));

  return valid;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kol::GoGraphImpl::GoGraphImpl(const GoTermMap& go_term_map) {

  size_t relation_count{0};

  for (auto const& [term_id, term_record_ptr] : go_term_map) {

    if (not term_record_ptr->validRecord()) {

      ExecEnv::log().error("GoGraphImpl::GoGraphImpl; invalid term record at term: {}",  term_record_ptr->termId());
      ExecEnv::log().error("GoGraphImpl::GoGraphImpl; id: {}, name: {}, def: {}, ontology: {}",
                           term_record_ptr->termId(), term_record_ptr->name(), term_record_ptr->definition(), GO::ontologyToString(term_record_ptr->ontology()));

    } else {

      insertTerm(term_record_ptr->termId(), term_record_ptr->name(), term_record_ptr->definition(), GO::ontologyToString(term_record_ptr->ontology()));
      for (auto const& [related_term, relation] : term_record_ptr->relations()) {
        //insert related terms, there are just stubs to be overwritten later on
        insertTerm(related_term, "name", "description", "ontology");

        //insert edge
        insertRelationship(term_record_ptr->termId(), related_term, GO::relationshipToString(relation));
        ++relation_count;

      }

    }

  }

//  ExecEnv::log().info("GoGraphImpl::GoGraphImpl: {}, relationships: {}", go_term_map.size(), relation_count);
  //call to initialize the graph's vertex to index maps
  initMaps();

}


//! Method to insert terms into the graph
/*!
  This method takes a go term, description, and ontology information (MOLECULAR_FUNCTION,BIOLOGICAL_PROCESS,CELLULAR_COMPONENT).
   The method will check if the term already exists in the graph then add the vertex
   or update the meta data accordingly. The parser can call this method without having
   to consider if terms have already been added or not.
*/

void kol::GoGraphImpl::insertTerm(const std::string &termId, const std::string &name, const std::string &description, const std::string &ontology) {

  //term already exists, update its information,
  if (_nameToIndex.find(termId) != _nameToIndex.end()) {

    std::size_t index = _nameToIndex[termId];

    //Term needs to be updated
    //If name is "name", this is a stub, no need to update
    if (name != "name") {

      _names.at(index) = name;
      _descriptions.at(index) = description;
      GoVertex V = boost::vertex(index, _goGraph);
      _goGraph[V].ontology = GO::ontologyStringToCode(ontology);

    }

  } else {
    //term is new and must be added
    //map termId to index
    _nameToIndex[termId] = boost::num_vertices(_goGraph);

    //add the vertex to the graph, get back a reference to the vertex, V
    GoVertex V = boost::add_vertex(_goGraph);

    //set the termId
    _goGraph[V].termId = termId;
    //set the ontology
    _goGraph[V].ontology = GO::ontologyStringToCode(ontology);

    //add name to list
    _names.push_back(name);
    //add description to list
    _descriptions.push_back(description);

  }

}


//! Method to insert relationship edges into the graph
/*!
  This method takes a parent term, child term, and relationshp type as arguments.
   The method will insert the edge into the graph, setting the relationship type based
   on the data provided.
*/
void kol::GoGraphImpl::insertRelationship(const std::string &termParent, const std::string &termChild, const std::string &relationship) {

  //get the vertices by name, they should already exit in the graph
  GoVertex v = boost::vertex(_nameToIndex[termParent], _goGraph);
  GoVertex u = boost::vertex(_nameToIndex[termChild], _goGraph);

  //get the relationship type as its enum value
  GO::Relationship relType = GO::relationshipStringToCode(relationship);

  //add the edge to the graph, get a reference to the edge
  std::pair<GoEdge, bool> myPair = boost::add_edge(v, u, _goGraph);

  //set that edge's internal value for relationship type
  GoEdge e = myPair.first;
  _goGraph[e].relType = relType;

}


//! A method to initialize internal index maps
/*!
  This method sets the private map variables by call calling boost get on the property maps.
*/
void kol::GoGraphImpl::initMaps() {

  _vMap = boost::get(boost::vertex_index, _goGraph);
  _eMap = boost::get(boost::edge_index, _goGraph);

}//end method initMaps


//! A helper method to return the index of the term
/*!
  This method returns the index of the given term.
*/
size_t kol::GoGraphImpl::getTermIndex(const std::string &term) const {

  if (not hasTerm(term)) {

    std::string error_message = "GoGraphImpl::getTermIndex; term: " + term + " not defined";
    throw std::runtime_error(error_message);

  }

  auto const&[term_key, index] = *(_nameToIndex.find(term));

  return index;

}

//! A helper method to return the string id based on the index
/*!
  This method returns the term's string id using its index.
    Used mainly for testing.
*/
std::string kol::GoGraphImpl::getTermStringIdByIndex(std::size_t index) const {

  GoVertex v = getVertexByIndex(index);
  return _goGraph[v].termId;

}


//! A helper method to return the string name based on the go term
/*!
  This method returns the term's string name using the go term.
*/
std::string kol::GoGraphImpl::getTermName(const std::string &term) const {

  if (hasTerm(term)) {

    size_t index = getTermIndex(term);
    return _names.at(index);

  } else {

    return "";

  }

}


//! A helper method to return the string description based on the go term
/*!
  This method returns the term's description string using the go term.
*/
std::string kol::GoGraphImpl::getTermDescription(const std::string &term) const {

  if (hasTerm(term)) {

    size_t index = getTermIndex(term);
    return _descriptions.at(index);

  } else {

    return "";

  }

}


//! A helper method to return the root of the graph
/*!
  This method returns the root vertex or first root vertex of a graph.
*/
kol::GoGraphImpl::GoVertex kol::GoGraphImpl::getRoot() const {

  //creat a vertex variable
  GoVertex root;
  auto [vertex_iter, vertex_end] = boost::vertices(_goGraph);
  for (; vertex_iter != vertex_end; ++vertex_iter) {
    //if it has no out edges it is a root
    if (boost::out_degree(*vertex_iter, _goGraph) == 0) {
      //set the variable and break the loop
      root = *vertex_iter;
      break;
    }
  }
  //return the root
  return root;

}

//! A helper method to return the ontology of a term by term string
/*!
  This method returns the term's ontology taking a string term as an argument
*/
kol::GO::Ontology kol::GoGraphImpl::getTermOntology(const std::string &term) const {

  if (not hasTerm(term)) {

    return GO::Ontology::ONTO_ERROR;

  } else {

    size_t index = getTermIndex(term);
    return _goGraph[index].ontology;

  }

}

//! A helper method to return the ontology of a term by index
/*!
  This method returns the term's ontoogy taking an index as an argument
*/
kol::GO::Ontology kol::GoGraphImpl::getTermOntologyByIndex(std::size_t index) const {

  GoVertex vertex = getVertexByIndex(index);
  return _goGraph[vertex].ontology;

}


//! A helper method to get the desendant terms for a given term.
/*!
  This method takes a term and returns a list of desendant terms.
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getDescendantTerms(const std::string &term) const {
  //return empty set, if term is not found
  if (!hasTerm(term)) {

    return OntologySetType<std::string>();

  }

  //get the correct index from the term string
  std::size_t vIndex = getTermIndex(term);

  //get the vertex from the term index
  GoVertex vertex = getVertexByIndex(vIndex);

  //create a map set
  OntologyMapType<std::size_t, bool> desendantMap;
  //call the recursive helper method.
  getDescendantTermsHelper(vertex, desendantMap);

  //create output container
  OntologySetType<std::string> desendantTerms;
  //create an iterator for the map
  for (auto const&[index, value] : desendantMap) {

    std::string index_term = getTermStringIdByIndex(index);
    desendantTerms.insert(index_term);

  }

  return desendantTerms;

}

kol::OntologySetType<std::string> kol::GoGraphImpl::getSelfDescendantTerms(const std::string &term) const {

  OntologySetType<std::string> self_ancestor_set = getDescendantTerms(term);
  self_ancestor_set.insert(term);
  return self_ancestor_set;

}


//! A helper method to get the ancestor terms for a given term.
/*!
  This method takes a term and returns a list of ancestor terms.
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getAncestorTerms(const std::string &term) const {
  //return empty set, if term is not found
  if (!hasTerm(term)) {

    return OntologySetType<std::string>();

  }

  //get the correct index from the term string
  std::size_t vIndex = getTermIndex(term);
  GoVertex vertex = getVertexByIndex(vIndex);

  //create a map set
  OntologyMapType<std::size_t, bool> ancestorMap;

  //call the recursive helper method.
  getAncestorTermsHelper(vertex, ancestorMap);

  //create output container
  OntologySetType<std::string> ancestorTerms;
  for (auto const&[term, value] : ancestorMap) {

    std::string index_term = getTermStringIdByIndex(term);
    ancestorTerms.insert(index_term);

  }

  return ancestorTerms;

}


kol::OntologySetType<std::string> kol::GoGraphImpl::getSelfAncestorTerms(const std::string &term) const {

  OntologySetType<std::string> self_ancestor_set = getAncestorTerms(term);
  self_ancestor_set.insert(term);
  return self_ancestor_set;

}

//! A method for calculating the extended term set. The set of all terms in the induced subgraph of the ontology
/*!
  This method returns the extended term set of a set of terms. Basically the set of terms and all their ancestors.
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getExtendedTermSet(const OntologySetType<std::string> &terms) const {


  // For each term create set of ancestor terms plus the term itself and push onto a vector
  std::vector<OntologySetType<std::string>> ancestor_set_vector;
  for (auto const &term : terms) {

    auto term_ancestors = getAncestorTerms(term);
    term_ancestors.insert(term);
    ancestor_set_vector.push_back(std::move(term_ancestors));

  }

  // Generate the union of all the ancestor sets.
  OntologySetType<std::string> inducedSet;
  for (auto const &term_ancestor_set : ancestor_set_vector) {

    // add the new terms to the set using union and the ancestors from the go graph
    inducedSet.insert(term_ancestor_set.begin(), term_ancestor_set.end());

  }

  return inducedSet;

}


//! A helper method to get the parent terms for a given term.
/*!
  This method takes a term and returns a list of parent terms (immediate ancestors).
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getParentTerms(const std::string &term) const {
  //return empty set, if term is not found
  if (!hasTerm(term)) {

    return OntologySetType<std::string>();

  } else {

    GoVertex vertex = getVertexByName(term);
    OntologySetType<std::string> parents;

    auto [parent_vertex_iter, parent_vertex_end] = boost::out_edges(vertex, _goGraph);
    for (; parent_vertex_iter != parent_vertex_end; ++parent_vertex_iter) {

      GoVertex parent_vertex = boost::target(*parent_vertex_iter, _goGraph);
      parents.insert(_goGraph[parent_vertex].termId);

    }

    return parents;

  }

}

//! A helper method to get the child terms for a given term.
/*!
  This method takes a term and returns a list of child terms.
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getChildTerms(const std::string &term) const {
  //return empty set, if term is not found
  if (!hasTerm(term)) {

    return OntologySetType<std::string>();

  } else {

    GoVertex vertex = getVertexByName(term);
    OntologySetType<std::string> children;

    auto [child_vertex_iter, child_vertex_end] = boost::in_edges(vertex, _goGraph);
    for (; child_vertex_iter != child_vertex_end; ++child_vertex_iter) {

      GoVertex child_vertex = boost::source(*child_vertex_iter, _goGraph);
      children.insert(_goGraph[child_vertex].termId);

    }

    return children;

  }

}

//! A helper method determines if a term is a leaf on the graph.
/*!
  This method takes a term and returns is true boolean if it is a leaf.
*/
bool kol::GoGraphImpl::isLeaf(const std::string &term) const {
  //return empty set, if term is not found
  if (!hasTerm(term)) {

    return false;

  } else {

    GoVertex vertex = getVertexByName(term);
    auto [vertex_begin, vertex_end] = boost::in_edges(vertex, _goGraph);
    return vertex_begin == vertex_end;

  }

}


kol::GoGraphImpl::GoVertex kol::GoGraphImpl::getRightLeaf(GoVertex vertex) const {

  do {

    auto [edge_begin, edge_end] = boost::in_edges(vertex, _goGraph);
    if (edge_begin == edge_end) {

      return vertex;

    }
    vertex = boost::target(*edge_begin, _goGraph);

  } while(true);

}


//! A helper method to retrieve all terms in the GoGraphImpl
/*!
  This method returns a set of term strings
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getAllTerms() const {

  //create a collection to return
  OntologySetType<std::string> outSet;
  for (std::size_t i = 0; i < getNumVertices(); ++i) {

    GoVertex vertex = getVertexByIndex(i);
    outSet.insert(_goGraph[vertex].termId);

  }
  return outSet;

}


//! A helper method to retrieve all terms in the GoGraphImpl
/*!
  This method returns a set of term strings
*/
kol::OntologyMapType<std::string, kol::GO::Ontology> kol::GoGraphImpl::getAllOntTerms() const {

  //create a collection to return
  OntologyMapType<std::string, GO::Ontology> outMap;
  for (std::size_t i = 0; i < getNumVertices(); ++i) {

    GoVertex vertex = getVertexByIndex(i);
    auto [iter, result] = outMap.try_emplace(_goGraph[vertex].termId, _goGraph[vertex].ontology);
    if (not result) {

      std::string term_id = _goGraph[vertex].termId;
      ExecEnv::log().error("GoGraphImpl::getAllOntTerms; duplicate GO term: {}", term_id);

    }

  }

  return outMap;

}


//! A helper method to retrieve all terms in the GoGraphImpl belonging to the BIOLOGICAL_PROCESS ontology
/*!
  This method returns a set of BIOLOGICAL_PROCESS terms in the graph
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getAllTermsBP() const {

  OntologySetType<std::string> outSet = getDescendantTerms(GO::getRootTermBP());
  outSet.insert(GO::getRootTermBP());
  return outSet;

}

//! A helper method to retrieve all terms in the GoGraphImpl belonging to the MOLECULAR_FUNCTION ontology
/*!
  This method returns a set of MOLECULAR_FUNCTION terms in the graph
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getAllTermsMF() const {

  OntologySetType<std::string> outSet = getDescendantTerms(GO::getRootTermMF());
  outSet.insert(GO::getRootTermMF());
  return outSet;

}

//! A helper method to retrieve all terms in the GoGraphImpl belonging to the CELLULAR_COMPONENT ontology
/*!
  This method returns a set of CELLULAR_COMPONENT terms in the graph
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getAllTermsCC() const {

  OntologySetType<std::string> outSet = getDescendantTerms(GO::getRootTermCC());
  outSet.insert(GO::getRootTermCC());
  return outSet;

}

//! A helper method to filter out all terms not belonging to a particular ontology
/*!
  This method returns a filtered set of ontology terms matching the given ontology
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::filterSetForOntology(const OntologySetType<std::string> &inSet, GO::Ontology onto) const {
  //create a collection to return
  OntologySetType<std::string> outSet;
  //iterate over the collection
  for (auto const &term : inSet) {

    if (getTermOntology(term) == onto) {

      outSet.insert(term);

    }

  }

  return outSet;

}

//! A helper method to filter out all terms not belonging to a particular ontology from a vector
/*!
  This method returns a filtered set of ontology terms matching the given ontology
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::filterSetForOntology(const std::vector<std::string> &inSet, GO::Ontology onto) const {
  //create a collection to return
  OntologySetType<std::string> outSet;

  //iterate over the collection
  std::vector<std::string>::const_iterator iter;
  for (iter = inSet.begin(); iter != inSet.end(); ++iter) {
    std::string term = *iter;

    if (getTermOntology(term) == onto) {
      outSet.insert(term);
    }
  }

  return outSet;

}

//! Get the root term for a particular term
/*!
  Return the root node for a term's ontology
*/
std::string kol::GoGraphImpl::getTermRoot(const std::string &term) const {

  GO::Ontology ontology = getTermOntology(term);

  switch (ontology) {
    case GO::Ontology::BIOLOGICAL_PROCESS:
      return GO::getRootTermBP();
    case GO::Ontology::MOLECULAR_FUNCTION:
      return GO::getRootTermMF();
    case GO::Ontology::CELLULAR_COMPONENT:
      return GO::getRootTermCC();
    default:
      return "";

  }

}


//!	A helper method to return only the terms of the give ontology.
/*!
  Returns only those terms used that occur for the given ontology.
*/
kol::OntologySetType<std::string> kol::GoGraphImpl::getOntologyTerms(GO::Ontology ontology) const {
  //Use only terms in the annotation database, this will save on space and computation time.
  std::string rootId;
  switch (ontology) {
    case GO::Ontology::BIOLOGICAL_PROCESS:
      rootId = GO::getRootTermBP();
      break;
    case GO::Ontology::MOLECULAR_FUNCTION:
      rootId = GO::getRootTermMF();
      break;
    case GO::Ontology::CELLULAR_COMPONENT:
      rootId = GO::getRootTermCC();
      break;
    default:
      rootId = "";
      break;
  }
  //All Ontology terms are descendants of the root.
  OntologySetType<std::string> ontologyTerms = getDescendantTerms(rootId);
  // Add the root.
  ontologyTerms.insert(rootId);

  return ontologyTerms;

}

//! A method to return the induced subgraph of a given term, ancestor graph
/*!
  This method returns a subgraph of the graph induced by traversing the ancestors of
    the given vertex.
*/
kol::GoGraphImpl::Graph& kol::GoGraphImpl::getInducedSubgraph2(const std::string &termId) {

  Graph &subgraph = _goGraph.create_subgraph();

  SubgraphBFSVisitor subgraphVisitor(subgraph);
  boost::breadth_first_search(_goGraph, getVertexByName(termId), boost::visitor(subgraphVisitor));
  return subgraph;

}

//! A method to return the induced subgraph of a given term, ancestor graph
/*!
  This method returns a subgraph of the graph induced by traversing the ancestors of
    the given vertex.
*/
kol::GoGraphImpl::Graph& kol::GoGraphImpl::getInducedSubgraph(const std::string &termId) {

  Graph &subgraph = _goGraph.create_subgraph();
  boost::add_vertex(getVertexByName(termId), subgraph);

  OntologySetType<std::string> ancestors = getAncestorTerms(termId);
  OntologySetType<std::string>::iterator iter;
  for (iter = ancestors.begin(); iter != ancestors.end(); ++iter) {

    boost::add_vertex(getVertexByName(*iter), subgraph);

  }

  return subgraph;
}

//! A method to calculate the number of connected components of the graph
/*!
  This method calculates the number of connected components in the graph.
    This is used to check if the GO graph contains only the 3 sub-ontologies.
*/
std::size_t kol::GoGraphImpl::getNumComponents() const {

  //Define undirected graph type
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> undirected_graph_t;
  undirected_graph_t undirected_g;

  //Make an undirected copy of the graph
  for (std::size_t i = 0; i < getNumVertices(); ++i) {
    //add graph vertices
    boost::add_vertex(undirected_g);
  }
  boost::graph_traits<Graph>::edge_iterator iter, end;
  boost::tie(iter, end) = boost::edges(_goGraph);
  for (; iter != end; ++iter) {
    GoVertex s = boost::source(*iter, _goGraph);
    GoVertex t = boost::target(*iter, _goGraph);
    //add edges, undirected
    boost::add_edge(_vMap[s], _vMap[t], undirected_g);
  }

  //calculate the connected components
  std::vector<std::size_t> componentAssignment(getNumVertices());
  return boost::connected_components(boost::make_reverse_graph(undirected_g), &componentAssignment[0]);

}


//! A private recursive helper method to get the descendant terms for a given term.
/*!
  This method is wrapped by a public method. It traverses the children of a node,
    populating the map with node indices of descendant terms.
*/
void kol::GoGraphImpl::getDescendantTermsHelper(GoVertex vertex, OntologyMapType<std::size_t, bool> &desendantMap) const {
  //create edge iterators
  InEdgeIterator ei, end;
  //loop over each edge
  for (boost::tie(ei, end) = boost::in_edges(vertex, _goGraph); ei != end; ++ei) {
    //get the soruce vertex ( specific --is_a--> general )
    GoVertex v = boost::source(*ei, _goGraph);
    //add the vertex index to the desendant map, addressed in the method call
    //  redundancies are handled by the map
    desendantMap[_vMap[v]] = true;
    //make the recursive call
    getDescendantTermsHelper(v, desendantMap);
  }
}

//! A private recursive helper method to get the ancestor terms for a given term.
/*!
  This method is wrapped by a public method. It traverses the parents of a node,
    populating the map with node indices of ancestor terms.
*/
void kol::GoGraphImpl::getAncestorTermsHelper(GoVertex vertex, OntologyMapType<std::size_t, bool> &ancestorMap) const {

  //create edge iterators
  OutEdgeIterator ei, end;

  //loop over each edge
  for (boost::tie(ei, end) = boost::out_edges(vertex, _goGraph); ei != end; ++ei) {
    //get the source vertex ( specific --is_a--> general )
    GoVertex v = boost::target(*ei, _goGraph);

    //add the vertex index to the descendant map, addressed in the method call
    //  redundancies are handled by the map
    ancestorMap[_vMap[v]] = true;

    //make the recursive call
    getAncestorTermsHelper(v, ancestorMap);
  }

}

