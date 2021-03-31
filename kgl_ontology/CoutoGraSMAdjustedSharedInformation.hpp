/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef GRASM_SHARED_INFORMATION_ADJUSTED
#define GRASM_SHARED_INFORMATION_ADJUSTED

#include <SharedInformationInterface.hpp>
#include <TermInformationContentMap.hpp>
#include <GoGraph.hpp>
#include <SetUtilities.hpp>
#include <Accumulators.hpp>

#include <utility>
#include <algorithm>

#include <boost/unordered_set.hpp>
#include <boost/graph/breadth_first_search.hpp>

/*! \class CoutoGraSMAdjustedSharedInformation
	\brief A class to calculate shared infromation accross disjoint common ancetors using an adjusted algorithm.

	This class calculates shared infromation accross disjoint common ancetors. This is a modificaiton of the
	 original algorithm provided by Couto. The adjustment changes the contrain to path lengths to strictly greater than.
	 See line 150.

    F. M. Couto, M. J. Silva, and P. M. Coutinho, "Measuring semantic similarity
	between Gene Ontology terms," Data & Knowledge Engineering, vol. 61, 
	pp. 137-152, Apr 2007.

	Couto proposing calculating this value a subsituite for the IC of the MICA in calculating
	 Resnik, Lin, and Jiang-Conrath

*/
class CoutoGraSMAdjustedSharedInformation : public SharedInformationInterface {

public:
	
	//! Constructor
	/*!
		Creates the CoutoGraSMAdjustedSharedInformation class
	*/
	CoutoGraSMAdjustedSharedInformation(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const TermInformationContentMap> icMap)
	: _goGraph(std::move(goGraph)), _icMap(std::move(icMap)) {}
  ~CoutoGraSMAdjustedSharedInformation() override = default;

	//! Calculate disjunctive ancestors.
	/*!
		A method for determining common disjunctive ancestors for two terms
	*/
	[[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1, const std::string &termC2) const {

		OntologySetType<std::string> ancestorsC1 = _goGraph->getAncestorTerms(termC1);
		ancestorsC1.insert(termC1);
		//std::cout << ancestorsC1.size() << std::endl;
    OntologySetType<std::string> ancestorsC2 = _goGraph->getAncestorTerms(termC2);
		ancestorsC2.insert(termC2);
		//std::cout << ancestorsC2.size() << std::endl;
		
		//Couto: CommonDisjAnc = {}
    OntologySetType<std::string> cda;

		if(termC1.compare(termC2) == 0){
			cda.insert(termC1);
			return cda;
		}

		//Couto: Anc = CommonAnc(c1,c2)
    OntologySetType<std::string> commonAncestors = SetUtilities::set_intersection(ancestorsC1,ancestorsC2);
		//std::cout << commonAncestors.size() << std::endl;

		std::vector<std::pair<double,std::string> > orderedCommonAncestors;
		//create a pair to associate a term with its information content
		for(auto const& term : commonAncestors) {

			orderedCommonAncestors.emplace_back(_icMap->getValue(term) , term);

		}

		//sort descending
		std::sort(orderedCommonAncestors.begin(),orderedCommonAncestors.end(),std::greater<>());

		
		//start of main algorithm
		//Couto: for all a in sortDescByIC(Anc) do ...
		for(auto const& [value, termA] : orderedCommonAncestors) {

			//Couto: isDisj=true
			bool isDisj = true;

			//std::cout << "testing " << termA << std::endl;

			//Couto: for all cda in CommonDisjAnc do ...
			for(auto const& termCda : cda) {

				//std::cout << "VS " << termCda << std::endl;

				//continue if the terms are the same
				if(termCda == termA){

					continue;

				}

				//Couto: isDisj = isDisj ^ ( DisjAnc(c1,(a,cda)) or DisjAnc(c2,(a,cda)) )
				isDisj = isDisj && (isDisjoint(termC1, termA, termCda) || isDisjoint(termC2, termA, termCda));

			}

			//Couto: if isDisj then...
			if(isDisj){
				//std::cout << myPair.second << " is cda " << std::endl;
				//Couto: addTo(CommonDisjAnc,a)
				cda.insert(termA);

			}

		}

		return cda;

	}


	//! Determine if a terms are disjoint in a concept.
	/*!
		A method for determining if, for a term c, a pair (a1,a2) is disjoint in c
	*/
	[[nodiscard]] bool isDisjoint(const std::string &termC, const std::string &termA1, const std::string &termA2) const {

		//std::cout << "isDisjoint " << termC << " ("  << termA1 << " , " << termA2 << ") "; //<< std::endl;
		//if not from same ontology, return 0;
		if(_goGraph->getTermOntology(termA1) != _goGraph->getTermOntology(termA2) ||
		   _goGraph->getTermOntology(termC) != _goGraph->getTermOntology(termA1) ||
		   _goGraph->getTermOntology(termC) != _goGraph->getTermOntology(termA2))
		{

			return false;

		}

		if(_icMap->getValue(termA1) <= _icMap->getValue(termA2)) {
			//std::cout << "case 1" << std::endl;
			size_t nPaths = getNumPaths(termA1,termA2);
			//std::cout << "nPaths " << termA1 << " to "  << termA2 << " " << nPaths << std::endl << std::endl;
			size_t nPaths1 = getNumPaths(termA1,termC);
			//std::cout << "nPaths " << termA1 << " to "  << termC << " " << nPaths1 << std::endl << std::endl;
			size_t nPaths2 = getNumPaths(termA2,termC);
			//std::cout << "nPaths " << termA2 << " to "  << termC << " " << nPaths2 << std::endl << std::endl;
			if(nPaths1 > nPaths*nPaths2) {
				//std::cout << "true" << std::endl;
				return true;

			} else {
				//std::cout << "false" << std::endl;
				return false;

			}
			//return nPaths1 > nPaths*nPaths2;
		} else {

			return false;

		}

	}

	//! Calculate the number of paths between two concept terms.
	/*!
		A method for calculating the number of paths from one term to another.
	*/
  [[nodiscard]] std::size_t getNumPaths(const std::string &termA, const std::string &termB) const {

		if(_icMap->getValue(termA) > _icMap->getValue(termB)) {

			return 0;

		}

		return pathCount(termA, termB);

	}

	//! Shared infromation between two conecepts.
	/*!
		A method for calculating the shared infromation between two concepts.
	*/
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override {
		// return 0 for any terms not in the datbase
		if (!_icMap->hasTerm(termA) || !_icMap->hasTerm(termB)){

			return 0.0;

		}
		// return 0 for terms in different ontologies
		if (_goGraph->getTermOntology(termA) != _goGraph->getTermOntology(termB)){

			return 0.0;

		}

		Accumulators::MeanAccumulator meanIC;
		OntologySetType<std::string> cda = getCommonDisjointAncestors( termA, termB);
		//std::cout << "size " << cda.size() << std::endl;

		for(auto const& term : cda) {
			//std::cout << _icMap[*iter] << std::endl;
			meanIC(_icMap->getValue(term));

		}

		return Accumulators::extractMean(meanIC);

	}

	//! Term information content.
	/*!
		An interface method to conventiently get information content of a single term
	*/
  [[nodiscard]] double sharedInformation(const std::string &term) const override {
		// return 0 for any terms not in the datbase

		if (not _icMap->hasTerm(term)){

			return 0.0;

		}

		return _icMap->getValue(term);

	}

	//! Maximum Ontology IC for normalization.
	/*!
		An interface method for returning the maximum information content for a term within a corpus for normalization purposes.
	*/
	[[nodiscard]] double maxInformationContent(const std::string &term) const override {

		double maxIC;
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
  [[nodiscard]] bool hasTerm(const std::string &term) const override {

		return _icMap->hasTerm(term);

	}

	//! An interface method for determining if the two terms are of like ontologies.
	/*!
		Determine if two terms are of the same ontology.
	*/
  [[nodiscard]] bool isSameOntology(const std::string &termA, const std::string &termB) const override {

		return _goGraph->getTermOntology(termA) == _goGraph->getTermOntology(termB);

	}


private:

	//! Count paths from B to A
	/*!
		Count paths between B and A
	*/
	[[nodiscard]] std::size_t pathCount(const std::string &termA, const std::string &termB) const {

		if (_icMap->getValue(termA) > _icMap->getValue(termB)) {

			return 0;

		}

		OntologySetType<std::string> ancestors = _goGraph->getAncestorTerms(termB);
		ancestors.insert(termB);

		const GoGraph::Graph& graph = _goGraph->getGraph();
		GoGraph::GoVertex root_vertex = _goGraph->getTermRootVertex(termB);

    OntologySetType<std::string> finished;
    OntologyMapType<std::string, size_t> pathMap;
		visitHelper(root_vertex, graph, ancestors, finished, pathMap);

		return pathMap[termA];

	}

	//! Recursive helper method that performs the DFS topological sort for path counting
	/*!
		A path counting topological sort recursive method.
	*/
	void visitHelper( const GoGraph::GoVertex &v,
                    const GoGraph::Graph& graph,
		                OntologySetType<std::string> &ancestors,
                    OntologySetType<std::string> &finished,
                    OntologyMapType<std::string, size_t> &pathMap) const
	{
		size_t childCount = 0;
		std::string vTerm = graph[v].termId;
		//std::cout << "discover vertex " << vTerm << std::endl;

		//examine children and recurse
		GoGraph::InEdgeIterator it, end;
		for (boost::tie(it, end) = boost::in_edges(v, graph); it != end; ++it){

			GoGraph::GoVertex child = boost::source(*it, graph);
			std::string childTerm = graph[child].termId;

			if (not SetUtilities::set_contains(ancestors, childTerm)){
				continue;
			}
			//recurse if child is not finished
			if (not SetUtilities::set_contains(finished, childTerm)){
				visitHelper(child, graph, ancestors, finished, pathMap);
			}
			++childCount;
		}

		//finish vertex
		finished.insert(vTerm);
		//std::cout << "finish vertex " << vTerm << ", childred " << childCount << std::endl;
		if (childCount == 0){

			pathMap[vTerm] = 1;

		} else{

			pathMap[vTerm] = 0;
			for (boost::tie(it, end) = boost::in_edges(v, graph); it != end; ++it){

				GoGraph::GoVertex child = boost::source(*it, graph);
				std::string childTerm = graph[child].termId;
				if (!SetUtilities::set_contains(ancestors, childTerm)){

					continue;

				}

				pathMap[vTerm] += pathMap[childTerm];

			}
		}
	}

	//! A private function to create a string key from a pair of terms
	/*!
		Creates a string key our of a pair to use in memorizing path counts
	*/
	[[nodiscard]] std::string keyPair(const std::string &termA, const std::string &termB) const {

		if (termA.compare(termB) > 0){

			return termB + "_" + termA;

		} else {

			return termB + "_" + termA;

		}

	}

	//! A private function to test if the key as been seen already
	/*!
	A private function to test if the key as been seen already.
	*/
	[[nodiscard]] bool hasSeenKey(const std::string &key) const {

		if (_pathMemory.find(key) != _pathMemory.end()){

			return true;

		}

		else{

			return false;

		}

	}

	std::shared_ptr<const GoGraph> _goGraph;
	std::shared_ptr<const TermInformationContentMap> _icMap;
	OntologyMapType<std::string, size_t> _pathMemory;

};
#endif
