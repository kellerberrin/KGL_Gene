/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef TERM_DEPTH_MAP
#define TERM_DEPTH_MAP

#include <GoGraph.hpp>
#include <AnnotationData.hpp>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>

/*! \class TermDepthMap
	\brief A class to calculate the depth of a GO term in the ontology.

	This class provides a map that returns the depth of a GO term. This method
	 is used in graph and edge based similarity methods to calculate a node's depth
*/
using TermDepthType = size_t;
class TermDepthMap{
public:
	//! A parameterized constructor
	/*!
		This constructor takes pointers to GoGraph and AnnotationData objects.
		  Only the parameterized construtor is allowed to ensure these objects are
		  created with valid parameters.
	*/
	explicit TermDepthMap(const GoGraph& graph) { initializeDepthMap(graph); }
	//! A default constructor
	/*!
		This constructor initialized the storage structures. Should not be used.
	*/
	TermDepthMap() = delete;
	//! A default desctructor
	/*!
		This desctructor clears the containters
	*/
	~TermDepthMap() = default;

	//! Accessor for probablities vector
	/*!
		Get the vector of values
	*/
	[[nodiscard]] const std::vector<TermDepthType>& getValues() const { return _depths; }

	//! Function to return all the keys in the map
	/*!
		Returns all valid keys in the map.
	*/
	[[nodiscard]] std::vector<std::string> getKeys() const {

		std::vector<std::string> keys;
		for (auto const& [key, value] : _nameToIndex)
		{

			keys.push_back(key);

		}

		return keys;

	}

	//! Method to test if the id exists in the map
	/*!
		Return true the id is found, false if not
	*/
	[[nodiscard]] bool hasTerm(const std::string& testTerm) const { return _nameToIndex.find(testTerm) != _nameToIndex.end(); }

	//! Overloaded [] bracket operator to mimic Map
	/*!
		This defines a bracket operator to access the data inside of the map.
		  This is done to mimic the behavior of the map class
	*/
	[[nodiscard]] TermDepthType operator[](const std::string& termId) const {
	//inline size_t& operator[](std::string termId){

		//return the depth
		return getValue(termId);

	}

	//! Mapping function to return the value mapped by key
	/*!
		Get the value mapped by the given key. A specified function for the [] operator
	*/
	[[nodiscard]] TermDepthType getValue(const std::string& termId) const {

    if (not hasTerm(termId) or termId.empty()) {

      return 0;

    }
		//get index
		auto const& [term, index] = *(_nameToIndex.find(termId));
		//return the depth
		return _depths.at(index);

	}

  //! A method for calculating the least common ancestor
  /*!
    This method searches the sets to determine the deepest common ancestor
  */
  [[nodiscard]] std::string getLCA(const OntologySetType<std::string> &ancestorsA, const OntologySetType<std::string> &ancestorsB) const {
    //get the first term as a start
    if(ancestorsA.size() < ancestorsB.size()){

      return getEfficientLCA( ancestorsA, ancestorsB);

    } else {

      return getEfficientLCA( ancestorsB, ancestorsA);

    }

  }


private:
	//! A private map that returns the index of a term.
	/*!
		This map takes string term ids and returns the index for annotation count access.
	*/
	OntologyMapType<std::string,std::size_t> _nameToIndex;

	//! A private list of term depths
	/*!
		This vector of doubles holds the depth for each term
	*/
	std::vector<TermDepthType> _depths;

  void initializeDepthMap(const GoGraph& graph) {

    //Initialize an annotation list the size of verticies in go, each value is 0
    //_depths = std::vector<std::size_t>(graph->getNumVertices(),0);
    _depths = std::vector<TermDepthType>(graph.getNumVertices(), 0);

    //get the (first) root of the ontology.
    //GoGraph::GoVertex root = graph->getRoot();
    //TESTING
    //std::cout << root << std::endl;
    //std::cout << graph->getTermStringIdByIndex(root) << std::endl;

    //get the boost graph from the GoGraph object. Must be done to utilize boost algorithms
    const GoGraph::Graph& go_graph = graph.getGraph();

    //wrap _depth with a vertex map
    const GoGraph::VertexIndexMap& vMap = graph.vertexIndexMap();

    /* // Temparary fix until I can get SWIG to recognize std::size_t
    boost::iterator_property_map< std::vector<std::size_t>::iterator,
                                GoGraph::VertexIndexMap >
                    d_map(_depths.begin(), vMap);
    */
    boost::iterator_property_map< std::vector<TermDepthType>::iterator, GoGraph::VertexIndexMap > d_map(_depths.begin(), vMap);

    //call the boost depth first search using our custom visitor
    // revering the graph is necessary otherwise the root vertex would have no edges.
    //boost::depth_first_search(boost::make_reverse_graph(*go_graph),boost::visitor(vis).root_vertex(root));
    GoGraph::GoVertex bpRoot = graph.getVertexByName(GO::getRootTermBP());
    GoGraph::GoVertex mfRoot = graph.getVertexByName(GO::getRootTermMF());
    GoGraph::GoVertex ccRoot = graph.getVertexByName(GO::getRootTermCC());

    //Start at bproot, record depths
    // must reverse graph due to edge relationship direction
    boost::breadth_first_search(boost::make_reverse_graph(go_graph),
                                bpRoot,
                                boost::visitor(boost::make_bfs_visitor(boost::record_distances(d_map, boost::on_tree_edge()))));

    //Start at bproot, record depths
    // must reverse graph due to edge relationship direction
    boost::breadth_first_search(boost::make_reverse_graph(go_graph),
                                mfRoot,
                                boost::visitor(boost::make_bfs_visitor(boost::record_distances(d_map, boost::on_tree_edge()))));

    //Start at bproot, record depths
    // must reverse graph due to edge relationship direction
    boost::breadth_first_search(boost::make_reverse_graph(go_graph),
                                ccRoot,
                                boost::visitor(boost::make_bfs_visitor(boost::record_distances(d_map, boost::on_tree_edge()))));

    //initialize the term to index map
    _nameToIndex = OntologyMapType<std::string,std::size_t>(boost::num_vertices(go_graph));

    // Vertex Iterators
    GoGraph::GoVertexIterator vi,vend;
    for(boost::tie(vi,vend) = boost::vertices(go_graph); vi != vend; ++vi){

      GoGraph::GoVertex v = *vi;
      _nameToIndex[graph.getTermStringIdByIndex(vMap[v])] = vMap[v];
      //std::cout << vMap[v] << " " << _depths.at(vMap[v]) << " " << graph->getTermNameByIndex(vMap[v]) << std::endl;
      //std::cin.get();
    }

  }

  //! A method for calculating the least common ancestor
  /*!
    This method searches the sets to determine the deepest common ancestor
  */
  [[nodiscard]] std::string getEfficientLCA(const OntologySetType<std::string>& smaller_set, const OntologySetType<std::string>& larger_set) const {
    //get the first term as a start
    std::string lca;
    //max depth
    TermDepthType max{0};

    //loop over shorter list
    for(auto const& currentTerm : smaller_set){

      if(larger_set.find(currentTerm) != smaller_set.end()){

        //if new max, update
        if(getValue(currentTerm) > max){

          lca = currentTerm;
          max = getValue(currentTerm);

        }

      }

    }

    return lca;

  }

};
#endif
