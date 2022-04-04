
#ifndef KOL_GO_GRAPH_IMPL
#define KOL_GO_GRAPH_IMPL

#include "kol_GoEnums.h"

#include "kol_OntologyTypes.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/reverse_graph.hpp>


namespace kellerberrin::ontology {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GoTermRecord {

public:

  GoTermRecord() = default;
  GoTermRecord(const GoTermRecord&) = default;
  ~GoTermRecord() = default;

  // Getters
  [[nodiscard]] const std::string& termId() const { return term_id_; }
  [[nodiscard]] const std::string& name() const { return name_; }
  [[nodiscard]] const std::string& definition() const { return definition_; }
  [[nodiscard]] GO::Ontology ontology() const { return ontology_; }
  [[nodiscard]] const std::vector<std::pair<std::string, GO::Relationship>>& relations() const { return relations_; }
  [[nodiscard]] const std::vector<std::string>& altId() const { return alt_id_; }
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& attributes() const { return attributes_; }

  // Setters.
  void termId(const std::string& term_id) { term_id_ = term_id; }
  void name(const std::string& name) { name_ = name; }
  void definition(const std::string& definition) { definition_ = definition; }
  void ontology(GO::Ontology ontology) { ontology_ = ontology; }
  void relations(std::pair<std::string, GO::Relationship> relation) { relations_.push_back(std::move(relation)); }
  void relations(std::vector<std::pair<std::string, GO::Relationship>>&& relations) { relations_ = relations; }
  void altId(const std::string& alt_id) { alt_id_.push_back(alt_id); }
  void attributes(std::pair<std::string, std::string> attribute) { attributes_.push_back(std::move(attribute)); }

  void clearRecord();
  [[nodiscard]] bool validRecord() const;

private:

  std::string term_id_;
  std::string name_;
  std::string definition_;
  GO::Ontology ontology_{GO::Ontology::ONTO_ERROR};
  std::vector<std::pair<std::string, GO::Relationship>> relations_;
  std::vector<std::string> alt_id_;
  std::vector<std::pair<std::string, std::string>> attributes_;  // key:value pairs is all misc attributes.

};

using GoTermMap = OntologyMapType<std::string, std::shared_ptr<const GoTermRecord>>;


//! A Vertex Property object
/*!
This struct represent the data needed by each vertex. Boost provides
 constant time access to these members by querying them using the vertex
 and graph objects (graphVar[vertex].termId ... etc.).
*/
struct VertexProps {

  VertexProps() = default;

  ~VertexProps() = default;

  /*!
    The term id of the go term, the GO accession, GO:########.
  */
  std::string termId;

  /*!
    The ontology of the GO term, BIOLOGICAL_PROCESS, MOLECULAR_FUNCTION, CELLULAR_COMPONENT
  */
  GO::Ontology ontology;

};

//! An Edge Property object
/*!
This struct represent the data needed by each edge. Boost provides
 constant time access to these members by querying them using the vertex
 and graph objects (graphVar[edge].relType).
*/
struct EdgeProps {

  EdgeProps() = default;

  ~EdgeProps() = default;

  /*!
    The type of relationship between the terms, is_a, part_of etc.
  */
  GO::Relationship relType;


};


/*! \class GoGraphImpl
	\brief This class holds the Gene Ontology directed acyclic graph

	This class holds the Gene Ontology as a boost graph. It provides the graph data,
	 as well as other structures which make working with the graph easier.

*/
class GoGraphImpl {

private:

  //! A Graph type representing Go
  /*!
  This typedef defines a graph type used as the basic go graph. This typedef
   takes VertexProps and EdgeProps as templete arguments.
*/
  using Graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
  boost::property<boost::vertex_index_t, size_t, VertexProps>,
  boost::property<boost::edge_index_t, size_t, EdgeProps> >;

public:



  //! The main Graph type representing Go
  /*!
  This typedef defines the main type as a subgraph of Graph_t. This allows the
   the graph to be divided into subgraphs if needed. Virtually not differnece but
   can cause problems with some boost constructors such as random graph generators.
*/
  using Graph = boost::subgraph<Graph_t>;

  ////////////////////////////
  //  Vertex and Edge typedefs
  ////////////////////////////
  //! A vertex object
  /*!
  A typedef of the boost vertex_descriptor. Saves typing by using GoVertex.
*/
  using GoVertex = boost::graph_traits<Graph>::vertex_descriptor;

  //! An edge object
  /*!
  A typedef of the boost edge_descriptor. Saves typing by using GoEdge.
*/
  using GoEdge = boost::graph_traits<Graph>::edge_descriptor;


  /////////////////////////////
  //  Vertex and Edge iterators
  /////////////////////////////
  //! A vertex iterator
  /*!
  A typedef of the boost vertex_iterator. Saves typing by using GoVertexIterator.
*/
  using GoVertexIterator = boost::graph_traits<Graph>::vertex_iterator;

  //! An in edge iterator
  /*!
  A typedef of the boost in_edge_iterator. Saves typing by using InEdgeIterator.
*/
  using InEdgeIterator = boost::graph_traits<Graph>::in_edge_iterator;

  //! An out edge iterator
  /*!
  A typedef of the boost out_edge_iterator. Saves typing by using OutEdgeIterator.
*/
  using OutEdgeIterator = boost::graph_traits<Graph>::out_edge_iterator;

  //////////////////////////////////
  //  Index maps for edge and vertex
  //////////////////////////////////
  //! A vertex to index map
  /*!
  A typedef of the boost property_map. Saves typing by using VertexIndexMap.
*/
  using VertexIndexMap = boost::property_map<Graph, boost::vertex_index_t>::type;

  //! An edge to index map
  /*!
  A typedef of the boost property_map. Saves typing by using EdgeIndexMap.
*/
  using EdgeIndexMap = boost::property_map<Graph, boost::edge_index_t>::type;


  GoGraphImpl() = default;
  explicit GoGraphImpl(const GoTermMap& go_term_map);
  ~GoGraphImpl() = default;


  //! A method to return the boost graph
  /*!
    This method is needed to return the graph to other boost algorithms. As stated by boost,
     subclasses of adjacency_list are not recommended.
     http://www.boost.org/doc/libs/1_54_0/libs/graph/doc/graph_concepts.html
  */
  [[nodiscard]] const Graph &getGraph() const { return _goGraph; }

  //! Vertex map
  /*!
    This maps a vertex object (GoVertex) to its index.
  */
  [[nodiscard]] const VertexIndexMap &vertexIndexMap() const { return _vMap; }

  //! Edge map
  /*!
    This maps an edge object (GoEdge) to its index.
  */
  [[nodiscard]] const EdgeIndexMap &edgeIndexMap() const { return _eMap; }


  //! Method to insert terms into the graph
  /*!
    This method takes a go term, description, and ontology information (MOLECULAR_FUNCTION,BIOLOGICAL_PROCESS,CELLULAR_COMPONENT).
     The method will check if the term already exists in the graph then add the vertex
     or update the meta data accordingly. The parser can call this method without having
     to consider if terms have already been added or not.
  */

  void insertTerm(const std::string &termId, const std::string &name, const std::string &description, const std::string &ontology);
  //! Method to insert relationship edges into the graph
  /*!
    This method takes a parent term, child term, and relationshp type as arguments.
     The method will insert the edge into the graph, setting the relationship type based
     on the data provided.
  */
  void insertRelationship(const std::string &termParent, const std::string &termChild, const std::string &relationship);


  //! Helper method to get number of vertices
  /*!
    This method calls boost num_vertices on the go graph
  */
  [[nodiscard]] size_t getNumVertices() const { return boost::num_vertices(_goGraph); }

  //! Helper method to get number of edges
  /*!
    This method calls boost num_edges on the go graph
  */
  [[nodiscard]] size_t getNumEdges() const { return boost::num_edges(_goGraph); }


  //! A method to initialize internal index maps
  /*!
    This method sets the private map variables by call calling boost get on the property maps.
  */
  void initMaps();



  //! A helper method to test term existence
  /*!
    Tests the map for existence of the term.
  */
  [[nodiscard]] bool hasTerm(const std::string &term) const { return _nameToIndex.find(term) != _nameToIndex.end(); }


  //! A helper method to return the index of the term
  /*!
    This method returns the index of the given term.
  */
  [[nodiscard]] size_t getTermIndex(const std::string &term) const;

  //! A helper method to return the string id based on the index
  /*!
    This method returns the term's string id using its index.
      Used mainly for testing.
  */
  [[nodiscard]] std::string getTermStringIdByIndex(std::size_t index) const;

  //! A helper method to return the string name based on the index
  /*!
    This method returns the term's string name using its index.
      Used mainly for testing.
  */
  [[nodiscard]] std::string getTermNameByIndex(std::size_t index) const { return _names.at(index); }

  //! A helper method to return the string name based on the go term
  /*!
    This method returns the term's string name using the go term.
  */
  [[nodiscard]] std::string getTermName(const std::string &term) const;

  //! A helper method to return the string description based on the index
  /*!
    This method returns the term's description string using its index.
      Used mainly for testing.
  */
  [[nodiscard]] std::string getTermDescriptionByIndex(std::size_t index) const {

    return _descriptions.at(index);

  }

  //! A helper method to return the string description based on the go term
  /*!
    This method returns the term's description string using the go term.
  */
  [[nodiscard]] std::string getTermDescription(const std::string &term) const;

  //! A helper method to return the root of the graph
  /*!
    This method returns the root vertex or first root vertex of a graph.
  */
  [[nodiscard]] GoVertex getRoot() const;

  //! A helper method to return the ontology of a term by term string
  /*!
    This method returns the term's ontoogy taking a string term as an argument
  */
  [[nodiscard]] GO::Ontology getTermOntology(const std::string &term) const;

  //! A helper method to return the ontology of a term by index
  /*!
    This method returns the term's ontoogy taking an index as an argument
  */
  [[nodiscard]] GO::Ontology getTermOntologyByIndex(std::size_t index) const;

  //! A helper method to return the ontology of a term by GoVertex
  /*!
    This method returns the term's ontoogy taking GoVertex as an argument
  */
  [[nodiscard]] GO::Ontology getTermOntologyByVertex(GoVertex vertex) const {

    return _goGraph[vertex].ontology;

  }

  //! A helper method to return the index of a GoVertex
  /*!
    This method returns the index of a GoVertex
  */
  [[nodiscard]] std::size_t getVertexIndex(GoVertex vertex) const {

    return _vMap[vertex];

  }

  //! A helper method to return the GoVertex for the given index
  /*!
    This method returns the GoVertex based on the given index
  */
  [[nodiscard]] GoVertex getVertexByIndex(std::size_t index) const {

    return boost::vertex(index, _goGraph);

  }

  //! A helper method to return the GoVertex for the given term
  /*!
    This method returns the GoVertex based on the given term
  */
  [[nodiscard]] GoVertex getVertexByName(const std::string &term) const {

    return boost::vertex(getTermIndex(term), _goGraph);

  }

  //! A helper method to get the desendant terms for a given term.
  /*!
    This method takes a term and returns a list of desendant terms.
  */
  [[nodiscard]] OntologySetType<std::string> getDescendantTerms(const std::string &term) const;

  [[nodiscard]] OntologySetType<std::string> getSelfDescendantTerms(const std::string &term) const;


  //! A helper method to get the ancestor terms for a given term.
  /*!
    This method takes a term and returns a list of ancestor terms.
  */
  [[nodiscard]] OntologySetType<std::string> getAncestorTerms(const std::string &term) const;

  [[nodiscard]] OntologySetType<std::string> getSelfAncestorTerms(const std::string &term) const;

  //! A method for calculating the extended term set. The set of all terms in the induced subgraph of the ontology
  /*!
    This method returns the extended term set of a set of terms. Basically the set of terms and all their ancestors.
  */
  [[nodiscard]] OntologySetType<std::string> getExtendedTermSet(const OntologySetType<std::string> &terms) const;

  //! A helper method to get the parent terms for a given term.
  /*!
    This method takes a term and returns a list of parent terms (immediate ancestors).
  */
  [[nodiscard]] OntologySetType<std::string> getParentTerms(const std::string &term) const;

  //! A helper method to get the child terms for a given term.
  /*!
    This method takes a term and returns a list of child terms.
  */
  [[nodiscard]] OntologySetType<std::string> getChildTerms(const std::string &term) const;

  //! A helper method determines if a term is a leaf on the graph.
  /*!
    This method takes a term and returns is true boolean if it is a leaf.
  */
  [[nodiscard]] bool isLeaf(const std::string &term) const;
  [[nodiscard]] GoVertex getRightLeaf(GoVertex vertex) const;

  //! A helper method to retrieve all terms in the GoGraphImpl
  /*!
    This method returns a set of term strings
  */
  [[nodiscard]] OntologySetType<std::string> getAllTerms() const;
  // Returns all terms and term ontology.
  [[nodiscard]] OntologyMapType<std::string, GO::Ontology> getAllOntTerms() const;

  //! A helper method to retrieve all terms in the GoGraphImpl belonging to the BIOLOGICAL_PROCESS ontology
  /*!
    This method returns a set of BIOLOGICAL_PROCESS terms in the graph
  */
  [[nodiscard]] OntologySetType<std::string> getAllTermsBP() const;


  //! A helper method to retrieve all terms in the GoGraphImpl belonging to the MOLECULAR_FUNCTION ontology
  /*!
    This method returns a set of MOLECULAR_FUNCTION terms in the graph
  */
  [[nodiscard]] OntologySetType<std::string> getAllTermsMF() const;

  //! A helper method to retrieve all terms in the GoGraphImpl belonging to the CELLULAR_COMPONENT ontology
  /*!
    This method returns a set of CELLULAR_COMPONENT terms in the graph
  */
  [[nodiscard]] OntologySetType<std::string> getAllTermsCC() const;

  //! A helper method to filter out all terms not belonging to a particular ontology
  /*!
    This method returns a filtered set of ontology terms matching the given ontology
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForOntology(const OntologySetType<std::string> &inSet, GO::Ontology onto) const;

  //! A helper method to filter out all terms not belonging to a particular ontology from a vector
  /*!
    This method returns a filtered set of ontology terms matching the given ontology
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForOntology(const std::vector<std::string> &inSet, GO::Ontology onto) const;

  //! Get the root term for a particular term
  /*!
    Return the root node for a term's ontology
  */
  [[nodiscard]] std::string getTermRoot(const std::string &term) const;

  //! Get the root vertex for a particular term
  /*!
    Get the root vertex for a particular term
  */
  [[nodiscard]] GoVertex getTermRootVertex(const std::string &term) const {

    return getVertexByName(getTermRoot(term));

  }


  //! A helper method to retrun only BIOLOGICAL_PROCESS terms from a vector
  /*!
    This method returns a filtered set containing only BIOLOGICAL_PROCESS terms
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForBP(const std::vector<std::string> &inSet) const {

    return filterSetForOntology(inSet, GO::Ontology::BIOLOGICAL_PROCESS);

  }

  //! A helper method to retrun only BIOLOGICAL_PROCESS terms from a set
  /*!
    This method returns a filtered set containing only BIOLOGICAL_PROCESS terms
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForBP(const OntologySetType<std::string> &inSet) const {

    return filterSetForOntology(inSet, GO::Ontology::BIOLOGICAL_PROCESS);

  }

  //! A helper method to retrun only MOLECULAR_FUNCTION terms from a vector
  /*!
    This method returns a filtered set containing only MOLECULAR_FUNCTION terms
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForMF(const std::vector<std::string> &inSet) const {

    return filterSetForOntology(inSet, GO::Ontology::MOLECULAR_FUNCTION);

  }

  //! A helper method to retrun only MOLECULAR_FUNCTION terms from a set
  /*!
    This method returns a filtered set containing only MOLECULAR_FUNCTION terms
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForMF(const OntologySetType<std::string> &inSet) const {

    return filterSetForOntology(inSet, GO::Ontology::MOLECULAR_FUNCTION);

  }

  //! A helper method to retrun only CELLULAR_COMPONENT terms from a vector
  /*!
    This method returns a filtered set containing only CELLULAR_COMPONENT terms
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForCC(const std::vector<std::string> &inSet) const {

    return filterSetForOntology(inSet, GO::Ontology::CELLULAR_COMPONENT);

  }

  //! A helper method to retrun only CELLULAR_COMPONENT terms from a set
  /*!
    This method returns a filtered set containing only CELLULAR_COMPONENT terms
  */
  [[nodiscard]] OntologySetType<std::string> filterSetForCC(const OntologySetType<std::string> &inSet) const {

    return filterSetForOntology(inSet, GO::Ontology::CELLULAR_COMPONENT);

  }


  //!	A helper method to return only the terms of the give ontology.
  /*!
    Returns only those terms used that occur for the given ontology.
  */
  [[nodiscard]] OntologySetType<std::string> getOntologyTerms(GO::Ontology ontology) const;

  //! A method to return the induced subgraph of a given term, ancestor graph
  /*!
    This method returns a subgraph of the graph induced by traversing the ancestors of
      the given vertex.
  */
  [[nodiscard]] Graph &getInducedSubgraph2(const std::string &termId);

  //! A method to return the induced subgraph of a given term, ancestor graph
  /*!
    This method returns a subgraph of the graph induced by traversing the ancestors of
      the given vertex.
  */
  [[nodiscard]] Graph &getInducedSubgraph(const std::string &termId);

  //! A method to calculate the number of connected components of the graph
  /*!
    This method calculates the number of connected components in the graph.
      This is used to check if the GO graph conatains only the 3 sub-ontologies.
  */
  [[nodiscard]] std::size_t getNumComponents() const;

private:

  //! Private graph member
  /*!
    The go graph defined as subgraph<Graph_t>.
  */
  Graph _goGraph;

  /////////////////////////////
  //  maps from vertex to index
  /////////////////////////////
  //! Private vertex map
  /*!
    This maps a vertex object (GoVertex) to its index.
  */
  VertexIndexMap _vMap;

  //! Private edge map
  /*!
    This maps an edge object (GoEdge) to its index.
  */
  EdgeIndexMap _eMap;

  //! A map from term name to index
  /*!
    This maps a term string to its index. Boost unordered_map has O(1) find like a hash map.
  */
  OntologyMapType<std::string, std::size_t> _nameToIndex;


  /*
    list of go short names and long descriptions (termId --> index --> name/description)
  */

  //! A list of term names, titles
  /*!
    A list of go names, the title of the term such as "positive regulation of cell cycle".
  */
  std::vector<std::string> _names;

  //! A list of term descriptions.
  /*!
    A list of go term descriptions. Long explanation and detailed description of the term.
  */
  std::vector<std::string> _descriptions;



  //! A private recursive helper method to get the desendant terms for a given term.
  /*!
    This method is wrapped by a public method. It traverses the children of a node,
      populating the map with node indices of desendant terms.
  */
  void getDescendantTermsHelper(GoVertex vertex, OntologyMapType<std::size_t, bool> &descendantMap) const;

  //! A private recursive helper method to get the ancestor terms for a given term.
  /*!
    This method is wrapped by a public method. It traverses the parents of a node,
      populating the map with node indices of ancestor terms.
  */
  void getAncestorTermsHelper(GoVertex vertex, OntologyMapType<std::size_t, bool> &ancestorMap) const;

  //! A Breath first search visitor that creates the induced subgraph of a graph
  /*!
    This class extends breadth first search and adds a vertex to a subgraph
      of the graph being visited.
  */


};


} // namespace


#endif
