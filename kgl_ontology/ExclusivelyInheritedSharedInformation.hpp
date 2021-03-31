/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef EXCLUSIVELY_INHERITED_SHARED_INFORMATION
#define EXCLUSIVELY_INHERITED_SHARED_INFORMATION

#include <TermInformationContentMap.hpp>
#include <SharedInformationInterface.hpp>
#include <SetUtilities.hpp>
#include <Accumulators.hpp>
#include <GoGraph.hpp>

#include <utility>
#include <algorithm>

#include <boost/unordered_set.hpp>
#include <boost/graph/breadth_first_search.hpp>

/*! \class ExclusivelyInheritedSharedInformation
	\brief A class to calculate shared infromation in linear time after Zhang and Lai.

	Shu-Bo Zhang and Jian-Huang Lai. Semantic Similarity measurement between gene ontology
	  terms based on exclusively inherited shared informaiton. Gene 558 (2015) 108-117.

*/
class ExclusivelyInheritedSharedInformation : public SharedInformationInterface{

public:
	
	//! A constructor
	/*!
		Creates the CoutoGraSMSharedInformation class
	*/
	ExclusivelyInheritedSharedInformation(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const TermInformationContentMap> icMap)
	: _goGraph(std::move(goGraph)), _icMap(std::move(icMap)) {}
  ~ExclusivelyInheritedSharedInformation() override = default;

	//! A method for determining the common disjunctive ancestors
	/*!
		This method returns the common disjunctive ancestors for two terms
	*/
	[[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1, const std::string &termC2) const {

		//Z&L: EICommonSet <- <empty_set>
    OntologySetType<std::string> cda;

		//EICA(t,t) = {t}
		if(termC1 == termC2) {

			cda.insert(termC1);
			return cda;

		}

		//Z&L: CommonAnSet <- GetCommonAnSet(t1,t2,AncestorSet)
    OntologySetType<std::string> ancestorsC1 = _goGraph->getAncestorTerms(termC1);
		ancestorsC1.insert(termC1);
    OntologySetType<std::string> ancestorsC2 = _goGraph->getAncestorTerms(termC2);
		ancestorsC2.insert(termC2);
    OntologySetType<std::string> commonAncestors = SetUtilities::set_intersection(ancestorsC1,ancestorsC2);
		//std::cout << commonAncestors.size() << std::endl;


		//commonDisjointAncestors(c,c) = {c}, by definition
		if(commonAncestors.size() == 1){

			return commonAncestors;

		}


		//Z&L: UnionAnSet <- GetAnSet(t1,AncestorSet) U GetAnSet(t2,AncestorSet)
    OntologySetType<std::string> unionAncestors = SetUtilities::set_union(ancestorsC1,ancestorsC2);

		//Z&L: DiffAnSet <- UnionAnSet - CommonAnSet
    OntologySetType<std::string> diffSet = SetUtilities::set_difference(unionAncestors,commonAncestors);

		//get the boost graph
		const GoGraph::Graph& g = _goGraph->getGraph();
		//Z&L: for each a in CommonAnSet do ...
		for(auto const& term : commonAncestors) {

			bool isDisj = false;

			//Z&L: DirectChildSet <- getDirectDescendant(a,ChildSet)
			//boost::unordered_set<std::string> directChildSet;
			//std::cout << "-----> " << term << std::endl;

			GoGraph::GoVertex v = _goGraph->getVertexByName(term);
			GoGraph::InEdgeIterator ei,end;
			for(boost::tie(ei,end) = boost::in_edges(v,g); ei != end; ++ei){

				GoGraph::GoVertex child = boost::source(*ei,g);
				size_t index = _goGraph->getVertexIndex(child);
				std::string cTerm = _goGraph->getTermStringIdByIndex(index);
				//std::cout << cTerm << std::endl;
				//directChildSet.insert(cTerm);
				
				if(diffSet.find(cTerm) != diffSet.end()){
					//early exit should improve runtime
					isDisj = true;
					break;

				}

			}

			//the below code was dropped to improve runtime, improvement is minor though

			//Z&L: tmpset <- DiffAnSet <intersect> DirectChildSet
			//boost::unordered_set<std::string> tempSet = SetUtilities::set_intersection(diffSet,directChildSet);
			//Z&L: if tmpset != <empty_set>
			//if(tempSet.size() != 0){
			//	isDisj = true;
			//}

			//if the term is a disjoint ancestor add it to the set
			if(isDisj){

				cda.insert(term);
				//std::cout << "eisi " << term << " " << _icMap[term] << std::endl;
			}

		}

		//std::cout << "eisi cda size "<< cda.size() << std::endl;
		return cda;
	
	}

	//! An method for returning the shared information of two terms
	/*!
		This method returns the mean information content of the frontier ancestors
	*/
	[[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override {
		// return 0 for any terms not in the datbase
		if (not _icMap->hasTerm(termA) or not _icMap->hasTerm(termB)){

			return 0.0;

		}
		// return 0 for terms in different ontologies
		if (_goGraph->getTermOntology(termA) != _goGraph->getTermOntology(termB)){

			return 0.0;

		}

		Accumulators::MeanAccumulator meanIC;
		OntologySetType<std::string> cda = getCommonDisjointAncestors( termA, termB);
		//std::cout << "size " << cda.size() << std::endl;

		for(auto const& element : cda){
			//std::cout << _icMap[*iter] << std::endl;
			meanIC(_icMap->getValue(element));

		}

		return Accumulators::extractMean(meanIC);

	}

	//! An interface method for returning the shared information of a single terms,or information content
	/*!
		This method privdes a mechanism for returing a term's infromation content.
	*/
	[[nodiscard]] double sharedInformation(const std::string &term) const override {
		// return 0 for any terms not in the datbase
		if (not _icMap->hasTerm(term)) {

			return 0.0;

		}

		return _icMap->getValue(term);

	}


	//! An interface method for returning the maximum information content for a term
	/*!
		This method provides the absolute max information content within a corpus for normalization purposes.
	*/
	[[nodiscard]] double maxInformationContent(const std::string &term) const override {

		double maxIC{0.0};
		//select the correct ontology normalization factor
		GO::Ontology ontoType = _goGraph->getTermOntology(term);
		if(ontoType == GO::Ontology::BIOLOGICAL_PROCESS){

			maxIC = -std::log(_icMap->getMinBP());

		}else if(ontoType == GO::Ontology::MOLECULAR_FUNCTION){

			maxIC = -std::log(_icMap->getMinMF());

		}else{

			maxIC = -std::log(_icMap->getMinCC());

		}

		return maxIC;
	}

	//! An interface method for determining if a term can be found
	/*!
	  Determines if the term can be found in the current map.
	*/
	[[nodiscard]] bool hasTerm(const std::string &term) const override { return _icMap->hasTerm(term); }

	//! An interface method for determining if the two terms are of like ontologies.
	/*!
		Determine if two terms are of the same ontology.
	*/
	[[nodiscard]] bool isSameOntology(const std::string &termA, const std::string &termB) const override {

		return _goGraph->getTermOntology(termA) == _goGraph->getTermOntology(termB);

	}

private:


	//breadth first search to calculate the visited edges
	class EdgeSetVisitor:public boost::default_bfs_visitor{

	public:
		EdgeSetVisitor( OntologySetType<std::size_t>& inSet,
                    const GoGraph::EdgeIndexMap& inMap,
                    OntologyMapType<std::string, OntologySetType<std::size_t> >& termToEdges)
                    : edgeSet(inSet),eMap(inMap),termEdgesMap(termToEdges){}


		template < typename Edge, typename Graph >
		void examine_edge(Edge e, const Graph & g)
		{
			//add the edge to the set of visited edges
			edgeSet.insert(eMap[e]);

			//get the vertex of the target
			typename Graph::vertex_descriptor v = boost::target(e,g);
			std::string term = g[v].termId;
			
			//create a set for the term if none exists
			if(termEdgesMap.find(term) == termEdgesMap.end()){
				termEdgesMap[term] = OntologySetType<std::size_t>();
			}
			//add the edge to the map
			termEdgesMap[term].insert(eMap[e]);


		}

		OntologySetType<std::size_t>& edgeSet;
		const GoGraph::EdgeIndexMap& eMap;
		OntologyMapType<std::string, OntologySetType<std::size_t> >& termEdgesMap;

	};

	std::shared_ptr<const GoGraph> _goGraph;
	std::shared_ptr<const TermInformationContentMap> _icMap;


};
#endif
