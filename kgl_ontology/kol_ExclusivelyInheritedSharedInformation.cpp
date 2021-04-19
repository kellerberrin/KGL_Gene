//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_ExclusivelyInheritedSharedInformation.h"

namespace kol = kellerberrin::ontology;


//! A method for determining the common disjunctive ancestors
/*!
  This method returns the common disjunctive ancestors for two terms
*/
kol::OntologySetType<std::string> kol::ExclusivelyInheritedSharedInformation::getCommonDisjointAncestors( const std::string &termC1,
                                                                                                          const std::string &termC2) const {

  //Z&L: EICommonSet <- <empty_set>
  OntologySetType<std::string> cda;

  //EICA(t,t) = {t}
  if (termC1 == termC2) {

    cda.insert(termC1);
    return cda;

  }

  //Z&L: CommonAnSet <- GetCommonAnSet(t1,t2,AncestorSet)
  OntologySetType<std::string> ancestorsC1 = graph_ptr_->getSelfAncestorTerms(termC1);
  OntologySetType<std::string> ancestorsC2 = graph_ptr_->getSelfAncestorTerms(termC2);
  OntologySetType<std::string> commonAncestors = SetUtilities::setIntersection(ancestorsC1, ancestorsC2);
  //std::cout << commonAncestors.size() << std::endl;


  //commonDisjointAncestors(c,c) = {c}, by definition
  if (commonAncestors.size() == 1) {

    return commonAncestors;

  }


  //Z&L: UnionAnSet <- GetAnSet(t1,AncestorSet) U GetAnSet(t2,AncestorSet)
  OntologySetType<std::string> unionAncestors = SetUtilities::setUnion(ancestorsC1, ancestorsC2);

  //Z&L: DiffAnSet <- UnionAnSet - CommonAnSet
  OntologySetType<std::string> diffSet = SetUtilities::setDifference(unionAncestors, commonAncestors);

  //get the boost graph
  const GoGraph::Graph &g = graph_ptr_->getGraph();
  //Z&L: for each a in CommonAnSet do ...
  for (auto const &term : commonAncestors) {

    bool isDisj = false;

    //Z&L: DirectChildSet <- getDirectDescendant(a,ChildSet)
    //boost::unordered_set<std::string> directChildSet;
    //std::cout << "-----> " << term << std::endl;

    GoGraph::GoVertex v = graph_ptr_->getVertexByName(term);
    GoGraph::InEdgeIterator ei, end;
    for (boost::tie(ei, end) = boost::in_edges(v, g); ei != end; ++ei) {

      GoGraph::GoVertex child = boost::source(*ei, g);
      size_t index = graph_ptr_->getVertexIndex(child);
      std::string cTerm = graph_ptr_->getTermStringIdByIndex(index);
      //std::cout << cTerm << std::endl;
      //directChildSet.insert(cTerm);

      if (diffSet.find(cTerm) != diffSet.end()) {
        //early exit should improve runtime
        isDisj = true;
        break;

      }

    }

    //the below code was dropped to improve runtime, improvement is minor though

    //Z&L: tmpset <- DiffAnSet <intersect> DirectChildSet
    //boost::unordered_set<std::string> tempSet = SetUtilities::setIntersection(diffSet,directChildSet);
    //Z&L: if tmpset != <empty_set>
    //if(tempSet.size() != 0){
    //	isDisj = true;
    //}

    //if the term is a disjoint ancestor add it to the set
    if (isDisj) {

      cda.insert(term);
      //std::cout << "eisi " << term << " " << ic_map_ptr_[term] << std::endl;
    }

  }

  //std::cout << "eisi cda size "<< cda.size() << std::endl;
  return cda;

}


//! An method for returning the shared information of two terms
/*!
  This method returns the mean information content of the frontier ancestors
*/
double kol::ExclusivelyInheritedSharedInformation::sharedInformation(const std::string &termA, const std::string &termB) const {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->hasTerm(termA) or not ic_map_ptr_->hasTerm(termB)) {

    return 0.0;

  }
  // return 0 for terms in different ontologies
  if (graph_ptr_->getTermOntology(termA) != graph_ptr_->getTermOntology(termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;
  OntologySetType<std::string> cda = getCommonDisjointAncestors(termA, termB);
  //std::cout << "size " << cda.size() << std::endl;

  for (auto const &element : cda) {
    //std::cout << ic_map_ptr_[*iter] << std::endl;
    meanIC(ic_map_ptr_->getValue(element));

  }

  return Accumulators::extractMean(meanIC);

}


//! An interface method for returning the shared information of a single terms,or information content
/*!
  This method privdes a mechanism for returing a term's infromation content.
*/

double kol::ExclusivelyInheritedSharedInformation::sharedInformation(const std::string &term) const  {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->hasTerm(term)) {

    return 0.0;

  }

  return ic_map_ptr_->getValue(term);

}


//! An interface method for returning the maximum information content for a term
/*!
  This method provides the absolute max information content within a corpus for normalization purposes.
*/
double kol::ExclusivelyInheritedSharedInformation::maxInformationContent(const std::string &term) const  {


  //select the correct ontology normalization factor
  GO::Ontology ontology = graph_ptr_->getTermOntology(term);
  double maxIC;

  switch (ontology) {

    case GO::Ontology::BIOLOGICAL_PROCESS:
      maxIC = ic_map_ptr_->getMinBP();
      break;

    case GO::Ontology::MOLECULAR_FUNCTION:
      maxIC = ic_map_ptr_->getMinMF();
      break;

    case GO::Ontology::CELLULAR_COMPONENT:
      maxIC = ic_map_ptr_->getMinCC();
      break;

    default:
    case GO::Ontology::ONTO_ERROR:
      maxIC = 0.0;
      break;

  }

  if (maxIC <= 0.0) {

    return 0.0;

  }

  return -1.0 * std::log(maxIC);

}

