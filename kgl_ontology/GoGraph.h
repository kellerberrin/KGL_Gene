/*=============================================================================
    Copyright (c) 2016 Paul W. Bible

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef GO_GRAPH
#define GO_GRAPH

#include <GoEnums.h>
#include <OntologyTypes.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/reverse_graph.hpp>


//! A Vertex Property object
  /*!
  This struct represent the data needed by each vertex. Boost provides
   constant time access to these members by querying them using the vertex
   and graph objects (graphVar[vertex].termId ... etc.).
*/
struct VertexProps{

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
struct EdgeProps{

  EdgeProps() = default;
  ~EdgeProps() = default;

  /*!
    The type of relationship between the terms, is_a, part_of etc.
  */
  GO::Relationship relType;


};


/*! \class GoGraph
	\brief This class holds the Gene Ontology directed acyclic graph

	This class holds the Gene Ontology as a boost graph. It provides the graph data,
	 as well as other strucutres which make working with the graph easier.

*/
class GoGraph{

public:


	//! A Graph type representing Go
    /*!
		This typedef defines a graph type used as the basic go graph. This typedef
		 takes VertexProps and EdgeProps as templete arguments.
	*/
	using Graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
	                                      boost::property< boost::vertex_index_t, size_t, VertexProps>,
			                                  boost::property< boost::edge_index_t, size_t, EdgeProps> >;

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
	using VertexIndexMap = boost::property_map<Graph, boost::vertex_index_t >::type;

	//! An edge to index map
    /*!
		A typedef of the boost property_map. Saves typing by using EdgeIndexMap.
	*/
	using EdgeIndexMap = boost::property_map<Graph, boost::edge_index_t >::type;


  GoGraph() = default;
	~GoGraph() = default;

  //! A method to return the boost graph
  /*!
    This method is needed to return the graph to other boost algorithms. As stated by boost,
     subclasses of adjacency_list are not recommended.
     http://www.boost.org/doc/libs/1_54_0/libs/graph/doc/graph_concepts.html
  */
  [[nodiscard]] const Graph& getGraph() const { return _goGraph; }

  //! Vertex map
  /*!
    This maps a vertex object (GoVertex) to its index.
  */
  [[nodiscard]] const VertexIndexMap& vertexIndexMap() const { return _vMap; }

  //! Edge map
  /*!
    This maps an edge object (GoEdge) to its index.
  */
  [[nodiscard]] const EdgeIndexMap& edgeIndexMap() const { return _eMap; }


  //! Method to insert terms into the graph
  /*!
    This method takes a go term, description, and ontology information (MOLECULAR_FUNCTION,BIOLOGICAL_PROCESS,CELLULAR_COMPONENT).
     The method will check if the term already exists in the graph then add the vertex
     or update the meta data accordingly. The parser can call this method without having
     to consider if terms have already been added or not.
  */

	void insertTerm(const std::string &termId, const std::string &name, const std::string &description, const std::string &ontology) {

		//term already exists, update its information,
		if(_nameToIndex.find(termId) != _nameToIndex.end()){

			std::size_t index = _nameToIndex[termId];

			//Term needs to be updated
			//If name is "name", this is a stub, no need to update
			if(name != "name") {

				_names.at(index) = name;
				_descriptions.at(index) = description;
				GoVertex V = boost::vertex(index,_goGraph);
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
	inline void insertRelationship(const std::string &termParent, const std::string &termChild, const std::string &relationship){

		//get the vertices by name, they should already exit in the graph 
		GoVertex v = boost::vertex(_nameToIndex[termParent],_goGraph);
		GoVertex u = boost::vertex(_nameToIndex[termChild],_goGraph);

		//get the relationship type as its enum value
		GO::Relationship relType = GO::relationshipStringToCode(relationship);

		//add the edge to the graph, get a reference to the edge
		std::pair<GoEdge,bool> myPair = boost::add_edge(v,u,_goGraph);

		//set that edge's internal value for relationship type
		GoEdge e = myPair.first;
		_goGraph[e].relType = relType;

	}


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
	void initMaps() {

		_vMap = boost::get(boost::vertex_index,_goGraph);
		_eMap = boost::get(boost::edge_index,_goGraph);

	}//end method initMaps



	//! A helper method to test term existence
	/*!
		Tests the map for existence of the term.
	*/
	[[nodiscard]] bool hasTerm(const std::string &term) const { return _nameToIndex.find(term) != _nameToIndex.end(); }


	//! A helper method to return the index of the term
	/*!
		This method returns the index of the given term.
	*/
	[[nodiscard]] size_t getTermIndex(const std::string &term) const
	{

    if (not hasTerm(term)) {

      std::string error_message = "GoGraph::getTermIndex; term: " + term + " not defined";
      throw std::runtime_error(error_message);

    }

    auto const& [term_key, index] = *(_nameToIndex.find(term));

	  return index;

	}

	//! A helper method to return the string id based on the index
	/*!
		This method returns the term's string id using its index.
		  Used mainly for testing.
	*/
	[[nodiscard]] std::string getTermStringIdByIndex(std::size_t index) const {

		GoVertex v = getVertexByIndex(index);
		return _goGraph[v].termId;

	}

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
	[[nodiscard]] std::string getTermName(const std::string& term) const {

		if (hasTerm(term)){

			size_t index = getTermIndex(term);
			return _names.at(index);

		}else{

			return "";

		}

	}

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
  [[nodiscard]] std::string getTermDescription(const std::string& term) const {

		if (hasTerm(term)){

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
  [[nodiscard]] GoVertex getRoot() const {

		//Create vertex iterators
		GoVertexIterator vi,vend;

		//creat a vertex variable
		GoVertex root;

		for(boost::tie(vi,vend) = boost::vertices(_goGraph); vi != vend; ++vi){
			//if it has no out edges it is a root
			if(boost::out_degree(*vi,_goGraph) == 0){
				//set the variable and break the loop
				root = *vi;
				break;
			}
		}
		//return the root
		return root;

	}

	//! A helper method to return the ontology of a term by term string
	/*!
		This method returns the term's ontoogy taking a string term as an argument
	*/
  [[nodiscard]] GO::Ontology getTermOntology(const std::string &term) const {

		if (not hasTerm(term)){

			return GO::Ontology::ONTO_ERROR;

		}else{

			return _goGraph[getTermIndex(term)].ontology;

		}

	}

	//! A helper method to return the ontology of a term by index
	/*!
		This method returns the term's ontoogy taking an index as an argument
	*/
  [[nodiscard]] GO::Ontology getTermOntologyByIndex(std::size_t index) const {

		GoVertex v = getVertexByIndex(index);
		return _goGraph[v].ontology;

	}

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

		return boost::vertex(index,_goGraph);

	}

	//! A helper method to return the GoVertex for the given term
	/*!
		This method returns the GoVertex based on the given term
	*/
  [[nodiscard]] GoVertex getVertexByName(const std::string &term) const {

		return boost::vertex(getTermIndex(term),_goGraph);

	}

	//! A helper method to get the desendant terms for a given term.
	/*!
		This method takes a term and returns a list of desendant terms.
	*/
  [[nodiscard]] OntologySetType<std::string> getDescendantTerms(const std::string &term) const {
		//return empty set, if term is not found
		if (!hasTerm(term)){

			return OntologySetType<std::string>();

		}

		//get the correct index from the term string
		std::size_t vIndex = getTermIndex(term);

		//get the vertex from the term index
		GoVertex vertex = getVertexByIndex(vIndex);

		//create a map set
		OntologyMapType<std::size_t,bool> desendantMap;
		//call the recursive helper method.
		getDescendantTermsHelper(vertex, desendantMap);

		//create output container
		OntologySetType<std::string> desendantTerms;
		//create an iterator for the map
		for(auto const& [index, value] : desendantMap) {

			std::string index_term = getTermStringIdByIndex(index);
			desendantTerms.insert(index_term);

		}

		return desendantTerms;

	}



	//! A helper method to get the ancestor terms for a given term.
	/*!
		This method takes a term and returns a list of ancestor terms.
	*/
  [[nodiscard]] OntologySetType<std::string> getAncestorTerms(const std::string &term) const {
		//return empty set, if term is not found
		if (!hasTerm(term)){

			return OntologySetType<std::string>();

		}

		//get the correct index from the term string
		std::size_t vIndex = getTermIndex(term);
		GoVertex vertex = getVertexByIndex(vIndex);

		//create a map set
		OntologyMapType<std::size_t,bool> ancestorMap;

		//call the recursive helper method.
		getAncestorTermsHelper(vertex, ancestorMap);

		//create output container
		OntologySetType<std::string> ancestorTerms;
		for(auto const& [term, value] : ancestorMap) {

			std::string index_term = getTermStringIdByIndex(term);
			ancestorTerms.insert(index_term);

		}

		return ancestorTerms;

	}


  [[nodiscard]] OntologySetType<std::string> getSelfAncestorTerms(const std::string &term) const {

    OntologySetType<std::string> self_ancestor_set = getAncestorTerms(term);
    self_ancestor_set.insert(term);
    return self_ancestor_set;

  }

    //! A method for calculating the extended term set. The set of all terms in the induced subgraph of the ontology
  /*!
    This method returns the extended term set of a set of terms. Basically the set of terms and all their ancestors.
  */
  [[nodiscard]] OntologySetType<std::string> getExtendedTermSet(const OntologySetType<std::string> &terms) const {


    // For each term create set of ancestor terms plus the term itself and push onto a vector
    std::vector<OntologySetType<std::string>> ancestor_set_vector;
    for (auto const& term : terms) {

      auto term_ancestors = getAncestorTerms(term);
      term_ancestors.insert(term);
      ancestor_set_vector.push_back(std::move(term_ancestors));

    }

    // Generate the union of all the ancestor sets.
    OntologySetType<std::string> inducedSet;
    for (auto const& term_ancestor_set : ancestor_set_vector) {

      // add the new terms to the set using union and the ancestors from the go graph
      inducedSet.insert(term_ancestor_set.begin(), term_ancestor_set.end());

    }

    return inducedSet;

  }


	//! A helper method to get the parent terms for a given term.
	/*!
		This method takes a term and returns a list of parent terms (immediate ancestors).
	*/
  [[nodiscard]] OntologySetType<std::string> getParentTerms(const std::string &term) const {
		//return empty set, if term is not found
		if (!hasTerm(term)){

			return OntologySetType<std::string>();

		} else{

			GoVertex vertex = getVertexByName(term);
			OntologySetType<std::string> parents;

			OutEdgeIterator ei, end;
			for (boost::tie(ei, end) = boost::out_edges(vertex, _goGraph); ei != end; ++ei){
				GoVertex v = boost::target(*ei, _goGraph);
				parents.insert(_goGraph[v].termId);
			}

			return parents;

		}

	}

	//! A helper method to get the child terms for a given term.
	/*!
		This method takes a term and returns a list of child terms.
	*/
  [[nodiscard]] OntologySetType<std::string> getChildTerms(const std::string &term) const {
		//return empty set, if term is not found
		if (!hasTerm(term)){

			return OntologySetType<std::string>();

		}else{

			GoVertex vertex = getVertexByName(term);
			OntologySetType<std::string> children;

			InEdgeIterator ei, end;
			for (boost::tie(ei, end) = boost::in_edges(vertex, _goGraph); ei != end; ++ei){
				GoVertex v = boost::source(*ei, _goGraph);
				children.insert(_goGraph[v].termId);
			}

			return children;

		}

	}

  //! A helper method determines if a term is a leaf on the graph.
  /*!
    This method takes a term and returns is true boolean if it is a leaf.
  */
  [[nodiscard]] bool isLeaf(const std::string &term) const {
    //return empty set, if term is not found
    if (!hasTerm(term)){

      return false;

    }else{

      GoVertex vertex = getVertexByName(term);

      InEdgeIterator ei, end;
      boost::tie(ei, end) = boost::in_edges(vertex, _goGraph);

      return ei == end;

    }

  }

	//! A helper method to retrieve all terms in the GoGraph
	/*!
		This method returns a set of term strings
	*/
  [[nodiscard]] OntologySetType<std::string> getAllTerms() const {

		//create a collection to return
		OntologySetType<std::string> outSet;
		for (std::size_t i = 0; i < getNumVertices(); ++i){
			GoVertex v = getVertexByIndex(i);
			outSet.insert(_goGraph[v].termId);
		}
		return outSet;

	}

	//! A helper method to retrieve all terms in the GoGraph belonging to the BIOLOGICAL_PROCESS ontology
	/*!
		This method returns a set of BIOLOGICAL_PROCESS terms in the graph
	*/
  [[nodiscard]] OntologySetType<std::string> getAllTermsBP() const {

		OntologySetType<std::string> outSet = getDescendantTerms(GO::getRootTermBP());
		outSet.insert(GO::getRootTermBP());
		return outSet;

	}

	//! A helper method to retrieve all terms in the GoGraph belonging to the MOLECULAR_FUNCTION ontology
	/*!
		This method returns a set of MOLECULAR_FUNCTION terms in the graph
	*/
  [[nodiscard]] OntologySetType<std::string> getAllTermsMF() const {

		OntologySetType<std::string> outSet = getDescendantTerms(GO::getRootTermMF());
		outSet.insert(GO::getRootTermMF());
		return outSet;

	}

	//! A helper method to retrieve all terms in the GoGraph belonging to the CELLULAR_COMPONENT ontology
	/*!
		This method returns a set of CELLULAR_COMPONENT terms in the graph
	*/
  [[nodiscard]] OntologySetType<std::string> getAllTermsCC() const {

		OntologySetType<std::string> outSet = getDescendantTerms(GO::getRootTermCC());
		outSet.insert(GO::getRootTermCC());
		return outSet;

	}

	//! A helper method to filter out all terms not belonging to a particular ontology
	/*!
		This method returns a filtered set of ontology terms matching the given ontology
	*/
  [[nodiscard]] OntologySetType<std::string> filterSetForOntology(const OntologySetType<std::string> &inSet, GO::Ontology onto) const {
		//create a collection to return
		OntologySetType<std::string> outSet;
		//iterate over the collection
		for(auto const& term : inSet) {

			if(getTermOntology(term) == onto){

				outSet.insert(term);

			}

		}

		return outSet;

	}

	//! A helper method to filter out all terms not belonging to a particular ontology from a vector
	/*!
		This method returns a filtered set of ontology terms matching the given ontology
	*/
  [[nodiscard]] OntologySetType<std::string> filterSetForOntology(const std::vector<std::string> &inSet, GO::Ontology onto) const {
		//create a collection to return
		OntologySetType<std::string> outSet;

		//iterate over the collection
		std::vector<std::string>::const_iterator iter;
		for (iter = inSet.begin(); iter != inSet.end(); ++iter){
			std::string term = *iter;

			if (getTermOntology(term) == onto){
				outSet.insert(term);
			}
		}

		return outSet;

	}

	//! Get the root term for a particular term
	/*!
		Return the root node for a term's ontology
	*/
  [[nodiscard]] std::string getTermRoot(const std::string &term) const {

		GO::Ontology ontology = getTermOntology(term);

		switch (ontology){
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
  [[nodiscard]] OntologySetType<std::string> getOntologyTerms(GO::Ontology ontology) const {
		//Use only terms in the annotation database, this will save on space and computation time.
		std::string rootId;
		switch (ontology){
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
  [[nodiscard]] Graph& getInducedSubgraph2(const std::string &termId) {
		
		Graph& subgraph = _goGraph.create_subgraph();

		SubgraphBFSVisitor subgraphVisitor(subgraph);
		boost::breadth_first_search(_goGraph,getVertexByName(termId),boost::visitor(subgraphVisitor));
		return subgraph;

	}

	//! A method to return the induced subgraph of a given term, ancestor graph
	/*!
		This method returns a subgraph of the graph induced by traversing the ancestors of 
		  the given vertex.
	*/
  [[nodiscard]] Graph& getInducedSubgraph(const std::string &termId) {
		
		Graph& subgraph = _goGraph.create_subgraph();
		boost::add_vertex(getVertexByName(termId), subgraph);

		OntologySetType<std::string> ancestors = getAncestorTerms(termId);
		OntologySetType<std::string>::iterator iter;
		for(iter = ancestors.begin(); iter != ancestors.end(); ++iter){

			boost::add_vertex(getVertexByName(*iter), subgraph);

		}

		return subgraph;
	}

	//! A method to calculate the number of connected components of the graph
	/*!
		This method calculates the number of connected components in the graph.
		  This is used to check if the GO graph conatains only the 3 sub-ontologies.
	*/
	[[nodiscard]] std::size_t getNumComponents() const {

		//Define undirected graph type
		typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS> undirected_graph_t;
		undirected_graph_t undirected_g;

		//Make an undirected copy of the graph
		for (std::size_t i = 0; i < getNumVertices(); ++i){
			//add graph vertices
			boost::add_vertex(undirected_g);
		}
		boost::graph_traits<Graph>::edge_iterator iter, end;
		boost::tie(iter, end) = boost::edges(_goGraph);
		for (; iter != end; ++iter){
			GoVertex s = boost::source(*iter, _goGraph);
			GoVertex t = boost::target(*iter, _goGraph);
			//add edges, undirected
			boost::add_edge(_vMap[s], _vMap[t], undirected_g);
		}

		//calculate the connected components
		std::vector<std::size_t> componentAssignment(getNumVertices());
		return boost::connected_components(boost::make_reverse_graph(undirected_g), &componentAssignment[0]);

	}




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
	EdgeIndexMap   _eMap;

	//! A map from term name to index
	/*!
		This maps a term string to its index. Boost unordered_map has O(1) find like a hash map.
	*/
	OntologyMapType<std::string,std::size_t> _nameToIndex;


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
	void getDescendantTermsHelper(GoVertex vertex, OntologyMapType<std::size_t,bool> &desendantMap) const {
		//create edge iterators
		InEdgeIterator ei,end;
		//loop over each edge
		for(boost::tie(ei,end) = boost::in_edges(vertex, _goGraph); ei != end; ++ei){
			//get the soruce vertex ( specific --is_a--> general )
			GoVertex v = boost::source(*ei,_goGraph);
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
	void getAncestorTermsHelper(GoVertex vertex, OntologyMapType<std::size_t,bool> &ancestorMap) const {

		//create edge iterators
		OutEdgeIterator ei,end;

		//loop over each edge
		for(boost::tie(ei,end) = boost::out_edges(vertex,_goGraph); ei != end; ++ei){
			//get the soruce vertex ( specific --is_a--> general )
			GoVertex v = boost::target(*ei,_goGraph);

			//add the vertex index to the desendant map, addressed in the method call
			//  redundancies are handled by the map
			ancestorMap[_vMap[v]] = true;

			//make the recursive call
			getAncestorTermsHelper(v,ancestorMap);
		}

	}

  //! A Breath first search visitor that creates the induced subgraph of a graph
  /*!
    This class extends breadth first search and adds a vertex to a subgraph
      of the graph being visited.
  */

  class SubgraphBFSVisitor:public boost::default_bfs_visitor{

  public:
    SubgraphBFSVisitor(GoGraph::Graph& sub): subgraph(sub) {}

    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph&)
    {
      boost::add_vertex(u, subgraph);
    }

    GoGraph::Graph& subgraph;

  };

};


#endif
